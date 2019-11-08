package org.broadinstitute.hellbender.engine.spark;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.Hidden;

import java.io.Serializable;

public class AssemblyRegionArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final String MIN_ASSEMBLY_LONG_NAME = "min-assembly-region-size";
    public static final String MAX_ASSEMBLY_LONG_NAME = "max-assembly-region-size";
    public static final String ASSEMBLY_PADDING_LONG_NAME = "assembly-region-padding";
    public static final String MAX_STARTS_LONG_NAME = "max-reads-per-alignment-start";
    public static final String THRESHOLD_LONG_NAME = "active-probability-threshold";
    public static final String DONT_TRIM_ACTIVE_REGIONS_LONG_NAME = "dont-trim-active-regions";
    public static final String PADDING_AROUND_VARIANTS_LONG_NAME = "padding-around-variants";

    //NOTE: many of these settings are referenced by HaplotypeCallerSpark
    public static final int DEFAULT_MIN_ASSEMBLY_REGION_SIZE = 50;
    public static final int DEFAULT_MAX_ASSEMBLY_REGION_SIZE = 300;
    public static final int DEFAULT_ASSEMBLY_REGION_PADDING = 100;
    public static final int DEFAULT_MAX_READS_PER_ALIGNMENT = 50;
    public static final double DEFAULT_ACTIVE_PROB_THRESHOLD = 0.002;



    @Argument(fullName = MIN_ASSEMBLY_LONG_NAME, doc = "Minimum size of an assembly region", optional = true)
    public int minAssemblyRegionSize = defaultMinAssemblyRegionSize();

    @Argument(fullName = MAX_ASSEMBLY_LONG_NAME, doc = "Maximum size of an assembly region", optional = true)
    public int maxAssemblyRegionSize = defaultMaxAssemblyRegionSize();

    @Argument(fullName = ASSEMBLY_PADDING_LONG_NAME, doc = "Number of additional bases of context to include around each assembly region", optional = true)
    public int assemblyRegionPadding = defaultAssemblyRegionPadding();

    @Argument(fullName = MAX_STARTS_LONG_NAME, doc = "Maximum number of reads to retain per alignment start position. Reads above this threshold will be downsampled. Set to 0 to disable.", optional = true)
    public int maxReadsPerAlignmentStart = defaultMaxReadsPerAlignmentStart();

    @Advanced
    @Argument(fullName = THRESHOLD_LONG_NAME, doc="Minimum probability for a locus to be considered active.", optional = true)
    public double activeProbThreshold = defaultActiveProbThreshold();

    @Advanced
    @Argument(fullName= DONT_TRIM_ACTIVE_REGIONS_LONG_NAME, doc="If specified, we will not trim down the active region from the full region (active + extension) to just the active interval for genotyping", optional = true)
    public boolean dontTrimActiveRegions = false;

    @Hidden
    @Argument(fullName= PADDING_AROUND_VARIANTS_LONG_NAME, doc = "Number of bases of padding around variants (after assembly and before genotyping)", optional = true)
    public int variantPadding = 20;

    protected int defaultMinAssemblyRegionSize() { return DEFAULT_MIN_ASSEMBLY_REGION_SIZE; }

    protected int defaultMaxAssemblyRegionSize() { return DEFAULT_MAX_ASSEMBLY_REGION_SIZE; }

    protected int defaultAssemblyRegionPadding() { return DEFAULT_ASSEMBLY_REGION_PADDING; }

    protected int defaultMaxReadsPerAlignmentStart() { return DEFAULT_MAX_READS_PER_ALIGNMENT; }

    protected double defaultActiveProbThreshold() { return DEFAULT_ACTIVE_PROB_THRESHOLD; }

    public void validate() {
        if ( minAssemblyRegionSize <= 0 || maxAssemblyRegionSize <= 0 ) {
            throw new CommandLineException.BadArgumentValue("min/max assembly region size must be > 0");
        }

        if ( minAssemblyRegionSize > maxAssemblyRegionSize ) {
            throw new CommandLineException.BadArgumentValue("minAssemblyRegionSize must be <= maxAssemblyRegionSize");
        }

        if ( assemblyRegionPadding < 0 ) {
            throw new CommandLineException.BadArgumentValue("assemblyRegionPadding must be >= 0");
        }

        if ( maxReadsPerAlignmentStart < 0 ) {
            throw new CommandLineException.BadArgumentValue("maxReadsPerAlignmentStart must be >= 0");
        }
    }
}
