package org.broadinstitute.hellbender.utils.activityprofile;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Captures the probability that a specific locus in the genome represents an "active" site containing
 * real variation.
 */
public final class ActivityProfileState {
    private final SimpleInterval loc;
    private double activeProb;
    private final Type type;

    public double getActiveProb() {
        return activeProb;
    }

    public Type getType() {
        return type;
    }

    public enum Type {
        NONE
    }

    /**
     * Create a new ActivityProfileState at loc with probability of being active of activeProb
     *
     * @param loc the position of the result profile (for debugging purposes)
     * @param activeProb the probability of being active (between 0 and 1)
     */
    public ActivityProfileState(final SimpleInterval loc, final double activeProb) {
        this(loc, activeProb, Type.NONE);
    }

    /**
     * Create a new ActivityProfileState at loc with probability of being active of activeProb that maintains some
     * information about the result state and value
     *
     * The only state value in use is HIGH_QUALITY_SOFT_CLIPS, and here the value is interpreted as the number
     * of bp affected by the soft clips.
     *  @param loc the position of the result profile (for debugging purposes)
     * @param activeProb the probability of being active (between 0 and 1)
     */
    public ActivityProfileState(final SimpleInterval loc, final double activeProb, final Type type) {
        Utils.validateArg(loc.size() == 1, () -> "Location for an ActivityProfileState must have to size 1 bp but saw " + loc);
        this.loc = loc;
        this.activeProb = activeProb;
        this.type = type;
    }

    /**
     * Get the locus associated with the ActivityProfileState
     * @return the locus associated with the ActivityProfileState as a SimpleInterval
     */
    public SimpleInterval getLoc() {
        return loc;
    }

    @Override
    public String toString() {
        return String.format("ActivityProfileState{loc=%s, activeProb=%d, type=%s}", loc, activeProb, type) ;
    }
}
