package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;

import java.util.*;

/**
 * Region of the genome that gets assembled by the local assembly engine.
 *
 * As AssemblyRegion is defined by two intervals -- a primary interval containing a territory for variant calling and a second,
 * extended, interval for assembly -- as well as the reads overlapping the extended interval.  Although we do not call variants in the extended interval,
 * assembling over a larger territory improves calls in the primary territory.
 *
 * This concept is complicated somewhat by the fact that these intervals are mutable and the fact that the AssemblyRegion onject lives on after
 * assembly during local realignment during PairHMM.  Here is an example of the life cycle of an AssemblyRegion:
 *
 * Suppose that the HaplotypeCaller engine finds an evidence for a het in a pileup at locus 400 -- that is, it produces
 * an {@code ActivityProfileState} with non-zero probability at site 400 and passes it to its {@code ActivityProfile}.
 * The {@code ActivityProfile} eventually produces an AssemblyRegion based on the {@code AssemblyRegionArgumentCollection} parameters.
 * Let's suppose that this initial region has primary span 350-450 and extended span 100 - 700.
 *
 * Next, the assembly engine assembles all reads that overlap the extended interval to find variant haplotypes and the variants
 * they contain.  The AssemblyRegion is then trimmed down to a new primary interval bound by all assembled variants within the original primary interval
 * and a new extended interval.  The amount of padding of the new extended interval around the variants depends on the needs of local realignment
 * and as such need not equal the original padding that was used for assembly.
 */
public final class AssemblyRegion implements Locatable {

    private final SAMFileHeader header;

    /**
     * The reads included in this assembly region.  May be empty upon creation, and expand / contract
     * as reads are added or removed from this region.
     */
    private final List<GATKRead> reads;

    /**
     * The primary of this assembly region over which we call variants
     */
    private final SimpleInterval primaryLoc;

    /**
     * The extended interval over which assembly is performed
     */
    private final SimpleInterval extendedLoc;

    /**
     * Does this region represent an active region (all isActiveProbs above threshold) or
     * an inactive region (all isActiveProbs below threshold)?
     */
    private boolean isActive;

    /**
     * Indicates whether the region has been finalized
     */
    private boolean hasBeenFinalized;

    /**
     * Create a new AssemblyRegion containing no reads
     *  @param primaryLoc the span of this active region
     * @param isActive indicates whether this is an active region, or an inactive one
     * @param padding the active region padding to use for this active region
     */
    public AssemblyRegion(final Locatable primaryLoc, final boolean isActive, final int padding, final SAMFileHeader header) {
        this(primaryLoc, isActive, makeExtendedLoc(primaryLoc, padding, header), header);
    }

    private static SimpleInterval makeExtendedLoc(final Locatable activeRegionLoc, final int padding, final SAMFileHeader header) {
        final String contig = activeRegionLoc.getContig();
        return IntervalUtils.trimIntervalToContig(contig, activeRegionLoc.getStart() - padding, activeRegionLoc.getEnd() + padding, header.getSequence(contig).getSequenceLength());
    }

    private AssemblyRegion(final Locatable primaryLoc, final boolean isActive, final Locatable extendedLoc, final SAMFileHeader header) {
        ParamUtils.isPositive( primaryLoc.getLengthOnReference(), "Active region cannot be of zero size");

        this.header = Utils.nonNull(header);
        reads = new ArrayList<>();
        this.primaryLoc = new SimpleInterval(Utils.nonNull(primaryLoc));
        this.isActive = isActive;
        this.extendedLoc = new SimpleInterval(extendedLoc);
    }

    @Override
    public String getContig() {
        return primaryLoc.getContig();
    }

    @Override
    public int getStart() {
        return primaryLoc.getStart();
    }

    @Override
    public int getEnd() {
        return primaryLoc.getEnd();
    }

    @Override
    public String toString() {
        return "AssemblyRegion "  + primaryLoc.toString() + " active?=" + isActive + " nReads=" + reads.size();
    }

    /**
     * Does this region represent an active region
     */
    public boolean isActive() {
        return isActive;
    }

    /**
     * Override activity state of the region
     *
     * Note: Changing the isActive state after construction is a debug-level operation that only engine classes
     * like AssemblyRegionWalker should be able to do
     *
     * @param value new activity state of this region
     */
    void setIsActive(final boolean value) {
        isActive = value;
    }

    /**
     * Get the raw span of this assembly region (excluding the extension)
     * @return a non-null SimpleInterval
     */
    public SimpleInterval getSpan() { return primaryLoc; }

    /**
     * Get the span of this assembly region including the extension value
     * @return a non-null SimpleInterval
     */
    public SimpleInterval getExtendedSpan() { return extendedLoc; }

    /**
     * Get an unmodifiable copy of the list of reads currently in this assembly region.
     *
     * The reads are sorted by their coordinate position.
     * @return an unmodifiable and inmutable copy of the reads in the assembly region.
    */
    public List<GATKRead> getReads(){
        return Collections.unmodifiableList(new ArrayList<>(reads));
    }

    /**
     * Returns the header for the reads in this region.
     */
    public SAMFileHeader getHeader(){
        return header;
    }

    /**
     * Trim this region to just the span, producing a new assembly region without any reads that has only
     * the extent of newExtend intersected with the current extent
     * @param span the new extend of the active region we want
     * @param extensionSize the extensionSize size we want for the newly trimmed active region
     * @return a non-null, empty assembly region
     */
    public AssemblyRegion trim(final SimpleInterval span, final int extensionSize) {
        Utils.nonNull(span, "Active region extent cannot be null");
        Utils.validateArg( extensionSize >= 0, "the extensionSize size must be 0 or greater");
        final int extendStart = Math.max(1,span.getStart() - extensionSize);
        final int maxStop = header.getSequence(span.getContig()).getSequenceLength();
        final int extendStop = Math.min(span.getEnd() + extensionSize, maxStop);
        final SimpleInterval extendedSpan = new SimpleInterval(span.getContig(), extendStart, extendStop);
        return trim(span, extendedSpan);
    }

    /**
     * Equivalent to trim(span,span).
     */
    public AssemblyRegion trim(final SimpleInterval span) {
        return trim(span, span);
    }

    /**
     * Trim this region to no more than the span and extended span
     *
     * For example, suppose this active region is
     *
     * Active:    100-200 with extended span 50-250
     * NewExtent: 150-225 and extended span 150-275
     *
     * Here we return a region from 150-200 (the intersection of the requested span with the original span)
     * with extended span 150-275 (there is no problem with enlarging the extended span)
     *
     * @param span the new extent of the active region we want
     * @return a non-null, empty active region
     */
    public AssemblyRegion trim(final SimpleInterval span, final SimpleInterval extendedSpan) {
        Utils.nonNull(span, "Active region extent cannot be null");
        Utils.nonNull(extendedSpan, "Active region extended span cannot be null");
        Utils.validateArg(extendedSpan.contains(span), "The requested extended span must fully contain the requested span");

        final SimpleInterval trimmedSpan = primaryLoc.intersect(span);

        final String contig = getContig();
        final SimpleInterval resultExtendedLoc = IntervalUtils.trimIntervalToContig(contig, Math.min(trimmedSpan.getStart(), extendedSpan.getStart()),
                Math.max(trimmedSpan.getEnd(), extendedSpan.getEnd()), header.getSequence(contig).getSequenceLength());

        final AssemblyRegion result = new AssemblyRegion( trimmedSpan, isActive, resultExtendedLoc, header );

        final List<GATKRead> myReads = getReads();
        final int resultExtendedLocStart = resultExtendedLoc.getStart();
        final int resultExtendedLocStop = resultExtendedLoc.getEnd();

        final List<GATKRead> trimmedReads = new ArrayList<>(myReads.size());
        for( final GATKRead read : myReads ) {
            final GATKRead clippedRead = ReadClipper.hardClipToRegion(read, resultExtendedLocStart, resultExtendedLocStop);
            if( result.readOverlapsRegion(clippedRead) && !clippedRead.isEmpty() ) {
                trimmedReads.add(clippedRead);
            }
        }
        result.clearReads();

        trimmedReads.sort(new ReadCoordinateComparator(header));
        result.addAll(trimmedReads);
        return result;
    }

    /**
     * Returns true if read would overlap the extended extent of this region
     * @param read the read we want to test
     * @return true if read can be added to this region, false otherwise
     */
    public boolean readOverlapsRegion(final GATKRead read) {
        if ( read.isEmpty() || read.getStart() > read.getEnd() ) {
            return false;
        }

        return read.overlaps(extendedLoc);
    }

    /**
     * Add read to this region
     *
     * Read must have alignment start >= than the last read currently in this active region.
     *
     * @throws IllegalArgumentException if read doesn't overlap the extended region of this active region
     *
     * @param read a non-null GATKRead
     */
    public void add( final GATKRead read ) {
        Utils.nonNull(read, "Read cannot be null");
        final SimpleInterval readLoc = new SimpleInterval( read );
        Utils.validateArg(readOverlapsRegion(read), () ->
                "Read location " + readLoc + " doesn't overlap with active region extended span " + extendedLoc);

        if ( ! reads.isEmpty() ) {
            final GATKRead lastRead = reads.get(size() - 1);
            Utils.validateArg(lastRead.contigsMatch(read), () ->
                    "Attempting to add a read to ActiveRegion not on the same contig as other reads: lastRead " + lastRead + " attempting to add " + read);
            Utils.validateArg( read.getStart() >= lastRead.getStart(), () ->
                    "Attempting to add a read to ActiveRegion out of order w.r.t. other reads: lastRead " + lastRead + " at " + lastRead.getStart() + " attempting to add " + read + " at " + read.getStart());
        }

        reads.add( read );
    }

    /**
     * Get the number of reads currently in this region
     * @return an integer >= 0
     */
    public int size() { return reads.size(); }

    /**
     * Clear all of the reads currently in this region
     */
    public void clearReads() {
        reads.clear();
    }

    /**
     * Remove all of the reads in readsToRemove from this region
     * @param readsToRemove the set of reads we want to remove
     */
    public void removeAll( final Collection<GATKRead> readsToRemove ) {
        Utils.nonNull(readsToRemove);
        reads.removeAll(readsToRemove);
    }

    /**
     * Add all readsToAdd to this region
     * @param readsToAdd a collection of readsToAdd to add to this active region
     */
    public void addAll(final Collection<GATKRead> readsToAdd){
        Utils.nonNull(readsToAdd).forEach(r -> add(r));
    }

    /**
     * Get the reference bases from referenceReader spanned by the extended location of this region,
     * including additional padding bp on either side.  If this expanded region would exceed the boundaries
     * of the active region's contig, the returned result will be truncated to only include on-genome reference
     * bases.
     *
     * @param referenceReader the source of the reference genome bases
     * @param padding the padding, in BP, we want to add to either side of this active region extended region
     * @param genomeLoc a non-null genome loc indicating the base span of the bp we'd like to get the reference for
     * @return a non-null array of bytes holding the reference bases in referenceReader
     */
    private static byte[] getReference(final ReferenceSequenceFile referenceReader, final int padding, final SimpleInterval genomeLoc) {
        Utils.nonNull(referenceReader, "referenceReader cannot be null");
        Utils.nonNull(genomeLoc, "genomeLoc cannot be null");
        Utils.validateArg( padding >= 0, () -> "padding must be a positive integer but got " + padding);
        Utils.validateArg( genomeLoc.size() > 0, () -> "GenomeLoc must have size > 0 but got " + genomeLoc);

        return referenceReader.getSubsequenceAt( genomeLoc.getContig(),
                Math.max(1, genomeLoc.getStart() - padding),
                Math.min(referenceReader.getSequenceDictionary().getSequence(genomeLoc.getContig()).getSequenceLength(), genomeLoc.getEnd() + padding) ).getBases();
    }

    /**
     * See {@link #getAssemblyRegionReference} with padding == 0
     */
    public byte[] getAssemblyRegionReference( final ReferenceSequenceFile referenceReader ) {
        return getAssemblyRegionReference(referenceReader, 0);
    }

    /**
     * Get the reference bases from referenceReader spanned by the extended location of this active region,
     * including additional padding bp on either side.  If this expanded region would exceed the boundaries
     * of the active region's contig, the returned result will be truncated to only include on-genome reference
     * bases
     *
     * @param referenceReader the source of the reference genome bases
     * @param padding the padding, in BP, we want to add to either side of this active region extended region
     * @return a non-null array of bytes holding the reference bases in referenceReader
     */
    public byte[] getAssemblyRegionReference(final ReferenceSequenceFile referenceReader, final int padding ) {
        return getReference(referenceReader, padding, extendedLoc);
    }

    /**
     * Is this region equal to other, excluding any reads in either region in the comparison
     * @param other the other active region we want to test
     * @return true if this region is equal, excluding any reads and derived values, to other
     */
    public boolean equalsIgnoreReads(final AssemblyRegion other) {
        return other != null && isActive == other.isActive &&
                primaryLoc.equals(other.primaryLoc) && extendedLoc.equals(other.extendedLoc);
    }

    public void setFinalized(final boolean value) {
        hasBeenFinalized = value;
    }

    public boolean isFinalized() {
        return hasBeenFinalized;
    }

}
