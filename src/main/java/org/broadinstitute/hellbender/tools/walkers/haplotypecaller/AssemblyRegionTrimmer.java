package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.SortedSet;
import java.util.stream.Collectors;

/**
 * Helper component to manage active region trimming
 *
 * <p/>
 * It receives the user arguments that controls trimming and also performs the trimming region calculation.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class AssemblyRegionTrimmer {

    private AssemblyRegionArgumentCollection assemblyRegionArgs;

    private SAMSequenceDictionary sequenceDictionary;

    /**
     * Initializes the trimmer.
     *
     * <p/>
     * This method should be called once and only once before any trimming is performed.
     *
     * @param assemblyRegionArgs user arguments for the trimmer
     * @param sequenceDictionary dictionary to determine the bounds of contigs
     * @throws IllegalStateException if this trim calculator has already been initialized.
     * @throws IllegalArgumentException if the input location parser is {@code null}.
     * @throws CommandLineException.BadArgumentValue if any of the user argument values is invalid.
     */
    public void initialize(final AssemblyRegionArgumentCollection assemblyRegionArgs, final SAMSequenceDictionary sequenceDictionary) {
        this.sequenceDictionary = sequenceDictionary;
        this.assemblyRegionArgs = assemblyRegionArgs;
    }

    /**
     * Holds the result of trimming.
     */
    public final class Result {

        /**
         * Holds the input active region.
         */
        protected final AssemblyRegion originalRegion;

        /**
         * Trimmed interval containing all callable variants in the input active region (not considering the extension).
         */
        protected final SimpleInterval variantSpan;

        /**
         * The trimmed variant region span including the extension.
         */
        protected final SimpleInterval extendedSpan;

        /**
         * Holds the flanking spans that do not contain the callable variants.
         * <p/>
         * The first element of the pair is the left (up-stream) non-variant flank, whereas the second element is
         * the right (down-stream) non-variant flank.
         */
        protected final Pair<SimpleInterval, SimpleInterval> nonVariantFlanks;

        /**
         * Holds the collection of callable events within the variant trimming region.
         */
        protected final List<VariantContext> callableEvents;

        /**
         * Creates a trimming result given all its properties.
         * @param originalRegion the original active region.
         * @param overlappingEvents contained callable variation events.
         * @param nonVariantFlanks pair of non-variant flank spans around the variant containing span.
         * @param extendedSpan final trimmed variant span including the extension.
         * @param variantSpan variant containing span without padding.
         */
        protected Result(final AssemblyRegion originalRegion,
                         final List<VariantContext> overlappingEvents,
                         final Pair<SimpleInterval, SimpleInterval> nonVariantFlanks,
                         final SimpleInterval extendedSpan,
                         final SimpleInterval variantSpan) {
            this.originalRegion = originalRegion;
            this.nonVariantFlanks = nonVariantFlanks;
            callableEvents = overlappingEvents;
            this.variantSpan = variantSpan;
            this.extendedSpan = extendedSpan;

            Utils.validateArg(extendedSpan == null || variantSpan == null || extendedSpan.contains(variantSpan), "the extended callable span must include the callable span");
        }

        public boolean isVariationPresent() {
            return !callableEvents.isEmpty();
        }

        /**
         * Returns the trimmed variant containing region
         *
         * @throws IllegalStateException if there is no variation detected.
         *
         * @return never {@code null}.
         */
        public AssemblyRegion getCallableRegion() {
            Utils.validate(isVariationPresent(), "there is no variation thus no variant region");
            return originalRegion.trim(variantSpan, extendedSpan);
        }

        /**
         *  Returns the trimmed out left non-variant region.
         *
         *  Notice that in case of no variation, the whole original region is considered the left flanking region.
         */
        public Optional<AssemblyRegion> nonVariantLeftFlankRegion() {
            return nonVariantFlanks.getLeft() == null ? Optional.empty() : Optional.of(originalRegion.trim(nonVariantFlanks.getLeft(), assemblyRegionArgs.assemblyRegionPadding));
        }

        /**
         *  Returns the trimmed out right non-variant region.
         */
        public Optional<AssemblyRegion> nonVariantRightFlankRegion() {
            return nonVariantFlanks.getRight() == null ? Optional.empty() : Optional.of(originalRegion.trim(nonVariantFlanks.getRight(), assemblyRegionArgs.assemblyRegionPadding));
        }

        public AssemblyResultSet trimAssemblyResult(final AssemblyResultSet untrimmedAssemblyResult) {
            final boolean needsTrimming = nonVariantFlanks.getLeft() != null || nonVariantFlanks.getRight() != null;
            return needsTrimming ? untrimmedAssemblyResult.trimTo(getCallableRegion()) : untrimmedAssemblyResult;
        }
    }

    /**
     * Creates a result indicating that there was no trimming to be done.
     */
    protected Result noTrimming(final AssemblyRegion region, final List<VariantContext> events) {
        final SimpleInterval targetRegionLoc = region.getSpan();
        final Result result = new Result(region, events,Pair.of(null, null), targetRegionLoc, targetRegionLoc);
        return result;
    }

    /**
     * Creates a result indicating that no variation was found.
     */
    protected Result noVariation(final AssemblyRegion region) {
        final Result result = new Result(region, Collections.emptyList(), Pair.of(region.getSpan(), null), null, null);
        return result;
    }

    /**
     * Returns a trimming result from which the variant-trimmed region and non-variant flanks can be recovered.
     *
     * @param region the genome location range to trim.
     * @param variants list of variants to trim to.  Those not overlapping with {@code region} are ignored.
     * @return never {@code null}.
     */
    public Result trim(final AssemblyRegion region, final SortedSet<VariantContext> variants) {
        final List<VariantContext> variantsInRegion = variants.stream().filter(region::overlaps).collect(Collectors.toList());

        if ( variantsInRegion.isEmpty() ) {
            return noVariation(region);
        } else if ( assemblyRegionArgs.dontTrimActiveRegions) {
            return noTrimming(region, variantsInRegion);
        }

        final int minStart = variantsInRegion.stream().mapToInt(VariantContext::getStart).min().getAsInt();
        final int maxEnd = variantsInRegion.stream().mapToInt(VariantContext::getEnd).max().getAsInt();
        final SimpleInterval variantSpan = new SimpleInterval(region.getContig(), minStart, maxEnd);

        // TODO: there is no reason why the padding used here, which adds space around variants for PairHMM, must equal
        // TODO: the padding around the active region originally used for assembly.  We should make a new argument.
        // TODO: furthermore, should it depend on, for example, whether any of the variants are long indels, or where STRs are?
        final SimpleInterval paddedVariantSpan = variantSpan.expandWithinContig(assemblyRegionArgs.variantPadding, sequenceDictionary);
        final Pair<SimpleInterval, SimpleInterval> nonVariantFlanks = getFlanks(region, paddedVariantSpan);

        return new Result(region, variantsInRegion, nonVariantFlanks, paddedVariantSpan, variantSpan);
    }

    /**
     * Calculates the left and right non-variant flanks to trim away.
     *
     * @param core the span, including any padding, of the core region containing variation
     */
    private Pair<SimpleInterval, SimpleInterval> getFlanks(final AssemblyRegion region, final SimpleInterval core) {
        final String contig = region.getContig();
        final SimpleInterval leftFlank = region.getStart() < core.getStart() ? new SimpleInterval(contig, region.getStart(), core.getStart() - 1) : null;
        final SimpleInterval rightFlank = region.getEnd() > core.getEnd() ? new SimpleInterval(contig, core.getEnd() + 1, region.getEnd()) : null;
        return Pair.of(leftFlank, rightFlank);
    }
}