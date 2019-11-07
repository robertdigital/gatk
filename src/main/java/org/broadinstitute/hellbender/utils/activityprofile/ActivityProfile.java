package org.broadinstitute.hellbender.utils.activityprofile;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.IntStream;

/**
 * Class holding information about per-base activity scores for
 * assembly region traversal
 */
public class ActivityProfile {
    protected final List<ActivityProfileState> stateList;
    private final double activeProbThreshold;
    protected Integer regionStart = null;
    private Integer regionStop = null;
    private SAMFileHeader samHeader;

    /**
     * Create a new empty ActivityProfile
     * @param activeProbThreshold threshold for the probability of an activity profile state being active
     */
    public ActivityProfile(final double activeProbThreshold, final SAMFileHeader header) {
        stateList = new ArrayList<>();
        this.activeProbThreshold = activeProbThreshold;
        this.samHeader = header;
    }

    @Override
    public String toString() {
        return String.format("ActivityProfile{start=%d, stop=%d}", regionStart, regionStop);
    }

    /**
     * Is this profile empty? (ie., does it contain no ActivityProfileStates?)
     * @return true if the profile is empty (ie., contains no ActivityProfileStates)
     */
    public boolean isEmpty() {
        return stateList.isEmpty();
    }

    public int getEnd() {
        return regionStop;
    }

    // --------------------------------------------------------------------------------
    //
    // routines to add states to a profile
    //
    // --------------------------------------------------------------------------------

    /**
     * Add the next ActivityProfileState to this profile.
     *
     * Must be contiguous with the previously added result, or an IllegalArgumentException will be thrown
     *
     * @param state a well-formed ActivityProfileState result to incorporate into this profile
     */
    public void add(final ActivityProfileState state) {
        Utils.nonNull(state);
        final SimpleInterval loc = state.getLoc();

        if ( regionStart == null ) {
            regionStart = loc.getStart();
            regionStop = loc.getStart();
        } else {
            Utils.validateArg( regionStop == loc.getStart() - 1, () ->
                    "Bad add call to ActivityProfile: loc " + loc + " not immediately after last loc " + regionStop);
            regionStop++;
        }

        stateList.add(state);
    }

    // --------------------------------------------------------------------------------
    //
    // routines to get active regions from the profile
    //
    // --------------------------------------------------------------------------------

    /**
     * Get the next completed assembly regions from this profile, and remove all states supporting them from this profile
     *
     * Takes the current profile and finds all of the active / inactive from the start of the profile that are
     * ready.  By ready we mean unable to have their probability modified any longer by future additions to the
     * profile.  The regions that are popped off the profile take their states with them, so the start of this
     * profile will always be after the end of the last region returned here.
     *
     * The regions are returned sorted by genomic position.
     *
     * This function may not return anything in the list, if no regions are ready
     *
     * No returned region will be larger than maxRegionSize.
     *
     * @param assemblyRegionExtension the extension value to provide to the constructed regions
     * @param minRegionSize the minimum region size, in the case where we have to cut up regions that are too large
     * @param maxRegionSize the maximize size of the returned region
     * @param atEndOfInterval if true, we'll return a region whose end isn't sufficiently far from the end of the
     *                        stateList.  Used to close out the active region when we've hit some kind of end (such
     *                        as the end of the contig)
     * @return a non-null list of active regions
     */
    public List<AssemblyRegion> popReadyAssemblyRegions( final int assemblyRegionExtension, final int minRegionSize, final int maxRegionSize, final boolean atEndOfInterval ) {
        Utils.validateArg(assemblyRegionExtension >= 0, () -> "assemblyRegionExtension must be >= 0 but got " + assemblyRegionExtension);
        Utils.validateArg( minRegionSize > 0, () -> "minRegionSize must be >= 1 but got " + minRegionSize);
        Utils.validateArg( maxRegionSize > 0, () -> "maxRegionSize must be >= 1 but got " + maxRegionSize);

        final List<AssemblyRegion> regions = new ArrayList<>();

        // TODO: max region size should be replaced with some max look-ahead distance parameter
        while ( !stateList.isEmpty() && (atEndOfInterval || stateList.size() >= maxRegionSize)) {
            final ActivityProfileState first = stateList.get(0);
            final boolean isActiveRegion = first.getActiveProb() > activeProbThreshold;

            final int offsetOfNextRegionEnd = findEndOfRegion(isActiveRegion, minRegionSize, maxRegionSize);

            // create an active region, discard the states from that region, and advance the start locus
            final List<ActivityProfileState> statesInRegion = stateList.subList(0, offsetOfNextRegionEnd + 1);
            final SimpleInterval regionLoc = new SimpleInterval(first.getLoc().getContig(), first.getLoc().getStart(), first.getLoc().getStart() + offsetOfNextRegionEnd);
            regions.add(new AssemblyRegion(regionLoc, statesInRegion, isActiveRegion, assemblyRegionExtension, samHeader));
            statesInRegion.clear();

            // update the start and stop locations as necessary
            if ( stateList.isEmpty() ) {
                regionStart = null;
                regionStop = null;
            } else {
                regionStart = stateList.get(0).getLoc().getStart();
            }

        }
        return regions;
    }

    /**
     * Find the end of the current region, returning the index into the element isActive element, or -1 if the region isn't done
     *
     * The current region is defined from the start of the stateList, looking for elements that have the same isActiveRegion
     * flag (i.e., if isActiveRegion is true we are looking for states with isActiveProb > threshold, or alternatively
     * for states < threshold).  The maximize size of the returned region is maxRegionSize.  If forceConversion is
     * true, then we'll return the region end even if this isn't safely beyond the max prob propagation distance.
     *
     * Note that if isActiveRegion is true, and we can construct an assembly region > maxRegionSize in bp, we
     * find the further local minimum within that max region, and cut the region there, under the constraint
     * that the resulting region must be at least minRegionSize in bp.
     *
     * @param isActiveRegion is the region we're looking for an active region or inactive region?
     * @param minRegionSize the minimum region size, in the case where we have to cut up regions that are too large
     * @param maxRegionSize the maximize size of the returned region
     * @return the index into stateList of the last element of this region, or -1 if it cannot be found
     */
    private int findEndOfRegion(final boolean isActiveRegion, final int minRegionSize, final int maxRegionSize) {
        final int initialEndOfRegion = findFirstActivityBoundary(isActiveRegion, maxRegionSize);

        return (isActiveRegion && initialEndOfRegion == maxRegionSize) ? findBestCutSite(initialEndOfRegion, minRegionSize) - 1
                : initialEndOfRegion - 1;
    }

    /**
     * If a block of active sites is larger than the maximum region size, determine the best place to cut it.
     *
     * This algorithm finds the global minimum probability state within the region [minRegionSize, endOfActiveRegion)
     * (exclusive of endOfActiveRegion), and returns the state index of that state.
     *
     * @param endOfActiveRegion the last state of the current active region (exclusive)
     * @param minRegionSize the minimum of the left-most region, after cutting
     * @return the index of state after the cut site (just like endOfActiveRegion)
     */
    private int findBestCutSite(final int endOfActiveRegion, final int minRegionSize) {
        Utils.validateArg(endOfActiveRegion >= minRegionSize, "endOfActiveRegion must be >= minRegionSize");
        Utils.validateArg(minRegionSize >= 0, "minRegionSize must be >= 0");

        // TODO: logic for finding best cut site goes here
    }

    /**
     * Find the first index into the state list where the state changes from active to inactive or vice versa
     *
     * Finds the largest contiguous stretch of states starting at the first state (index 0) with the same isActive
     * state as isActiveRegion.  If the entire state list has the same isActive value, then returns
     * maxRegionSize
     *
     * @param isActiveRegion are we looking for a stretch of active states, or inactive ones?
     * @param maxRegionSize don't look for a boundary that would yield a region of size > maxRegionSize
     * @return the index of the first state in the state list with isActive value != isActiveRegion, or maxRegionSize
     *         if no such element exists
     */
    private int findFirstActivityBoundary(final boolean isActiveRegion, final int maxRegionSize) {
        final int max = Math.max(stateList.size(), maxRegionSize);

        // TODO: replace the state probability condition with proximity to active site condition
        return IntStream.range(0, max)
                .filter(n -> getProb(n) > activeProbThreshold != isActiveRegion)
                .findFirst().orElse(max);
    }

    /**
     * Helper function to get the probability of the state at offset index
     * @param index a valid offset into the state list
     * @return the isActiveProb of the state at index
     */
    private double getProb(final int index) {
        Utils.validIndex(index, stateList.size());

        return stateList.get(index).getActiveProb();
    }
}
