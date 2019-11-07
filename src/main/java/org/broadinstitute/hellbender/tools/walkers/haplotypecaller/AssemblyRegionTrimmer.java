package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
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

    /**
     * Holds the debug flag. If {@code true} the trimmer will output debugging messages into the log.
     */
    private boolean debug;

    /**
     * Holds the extension to be used
     */
    private int usableExtension;

    private ReadThreadingAssemblerArgumentCollection assemblyArgs;

    private SAMSequenceDictionary sequenceDictionary;

    /**
     * Holds a reference the trimmer logger.
     */
    private static final Logger logger = LogManager.getLogger(AssemblyRegionTrimmer.class);

    /**
     * Initializes the trimmer.
     *
     * <p/>
     * This method should be called once and only once before any trimming is performed.
     *
     * @param assemblyArgs user arguments for the trimmer
     * @param sequenceDictionary dictionary to determine the bounds of contigs
     * @throws IllegalStateException if this trim calculator has already been initialized.
     * @throws IllegalArgumentException if the input location parser is {@code null}.
     * @throws CommandLineException.BadArgumentValue if any of the user argument values is invalid.
     */
    public void initialize(final ReadThreadingAssemblerArgumentCollection assemblyArgs, final SAMSequenceDictionary sequenceDictionary) {
        Utils.validate(this.assemblyArgs == null, () -> getClass().getSimpleName() + " instance initialized twice");

        this.assemblyArgs = Utils.nonNull(assemblyArgs);;
        this.sequenceDictionary = sequenceDictionary;

        checkUserArguments();
        this.debug = assemblyArgs.debugAssembly;
        usableExtension = this.assemblyArgs.extension;
    }

    /**
     * Checks user trimming argument values
     *
     * @throws CommandLineException.BadArgumentValue if there is some problem with any of the arguments values.
     */
    private void checkUserArguments() {
        if ( assemblyArgs.snpPadding < 0 ) {
            throw new CommandLineException.BadArgumentValue("paddingAroundSNPs", "" + assemblyArgs.snpPadding + "< 0");
        }
        if ( assemblyArgs.indelPadding < 0 ) {
            throw new CommandLineException.BadArgumentValue("paddingAroundIndels", "" + assemblyArgs.indelPadding + "< 0");
        }
        if ( assemblyArgs.extension < 0) {
            throw new CommandLineException.BadArgumentValue("maxDiscARExtension", "" + assemblyArgs.extension + "< 0");
        }
    }

    /**
     * Holds the result of trimming.
     */
    public static final class Result {

        /**
         * Holds the input active region.
         */
        protected final AssemblyRegion originalRegion;

        /**
         * Holds the smaller range that contain all relevant callable variants in the
         * input active region (not considering the extension).
         *
         */
        protected final SimpleInterval callableSpan;

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
         * Holds variant-containing callable region.
         * <p/>
         * This is lazy-initialized using {@link #callableSpan}.
         */
        protected AssemblyRegion callableRegion;


        /**
         * Non-variant left flank region.
         * <p/>
         * This is lazy-initialized using
         * {@link #nonVariantFlanks}.{@link Pair#getLeft()} () getFirst()}.
         */
        private AssemblyRegion leftFlankRegion;

        /**
         * Non-variant right flank region.
         * <p/>
         * This is lazy-initialized using
         * {@link #nonVariantFlanks}.{@link Pair#getLeft()} () getSecond()}.
         */
        private AssemblyRegion rightFlankRegion;

        /**
         * Creates a trimming result given all its properties.
         * @param originalRegion the original active region.
         * @param overlappingEvents contained callable variation events.
         * @param nonVariantFlanks pair of non-variant flank spans around the variant containing span.
         * @param extendedSpan final trimmed variant span including the extension.
         * @param callableSpan variant containing span without padding.
         */
        protected Result(final AssemblyRegion originalRegion,
                         final List<VariantContext> overlappingEvents,
                         final Pair<SimpleInterval, SimpleInterval> nonVariantFlanks,
                         final SimpleInterval extendedSpan,
                         final SimpleInterval callableSpan) {
            this.originalRegion = originalRegion;
            this.nonVariantFlanks = nonVariantFlanks;
            callableEvents = overlappingEvents;
            this.callableSpan = callableSpan;
            this.extendedSpan = extendedSpan;

            Utils.validateArg(extendedSpan == null || callableSpan == null || extendedSpan.contains(callableSpan), "the extended callable span must include the callable span");
        }


        /**
         * Checks whether there is any variation present in the target region.
         *
         * @return {@code true} if there is any variant, {@code false} otherwise.
         */
        public boolean isVariationPresent() {
            return ! callableEvents.isEmpty();
        }

        /**
         * Checks whether the active region needs trimming.
         */
        public boolean needsTrimming() {
            return hasLeftFlankingRegion() || hasRightFlankingRegion();
        }

        /**
         * Returns the trimmed variant containing region
         *
         * @throws IllegalStateException if there is no variation detected.
         *
         * @return never {@code null}.
         */
        public AssemblyRegion getCallableRegion() {
            if (callableRegion == null && extendedSpan != null) {
                callableRegion = originalRegion.trim(callableSpan, extendedSpan);
            } else if (extendedSpan == null) {
                throw new IllegalStateException("there is no variation thus no variant region");
            }
            return callableRegion;
        }

        /**
         * Checks whether there is a non-empty left flanking non-variant trimmed out region.
         * @return {@code true} if there is a non-trivial left flank region, {@code false} otherwise.
         */
        public boolean hasLeftFlankingRegion() {
            return nonVariantFlanks.getLeft() != null;
        }

        /**
         * Checks whether there is a non-empty right flanking non-variant trimmed out region.
         * @return {@code true} if there is a non-trivial right flank region, {@code false} otherwise.
         */
        public boolean hasRightFlankingRegion() {
            return nonVariantFlanks.getRight() != null;
        }

        /**
         *  Returns the trimmed out left non-variant region.
         *  <p/>
         *  Notice that in case of no variation, the whole original region is considered the left flanking region.
         *
         *  @throws IllegalStateException if there is not such as left flanking region.
         */
        public AssemblyRegion nonVariantLeftFlankRegion() {
            if (leftFlankRegion == null && nonVariantFlanks.getLeft() != null) {
                leftFlankRegion = originalRegion.trim(nonVariantFlanks.getLeft(), originalRegion.getExtension());
            } else if (nonVariantFlanks.getLeft() == null) {
                throw new IllegalStateException("there is no left flank non-variant trimmed out region");
            }
            return leftFlankRegion;
        }

        /**
         *  Returns the trimmed out right non-variant region.
         */
        public AssemblyRegion nonVariantRightFlankRegion() {
            if (rightFlankRegion == null && nonVariantFlanks.getRight() != null) {
                rightFlankRegion = originalRegion.trim(nonVariantFlanks.getRight(), originalRegion.getExtension());
            } else if (nonVariantFlanks.getRight() == null) {
                throw new IllegalStateException("there is no right flank non-variant trimmed out region");
            }
            return rightFlankRegion;
        }

        /**
         * Creates a result indicating that there was no trimming to be done.
         */
        protected static Result noTrimming(final AssemblyRegion targetRegion,
                                           final List<VariantContext> events) {
            final SimpleInterval targetRegionLoc = targetRegion.getSpan();
            final Result result = new Result(targetRegion, events,Pair.of(null, null), targetRegionLoc, targetRegionLoc);
            result.callableRegion = targetRegion;
            return result;
        }

        /**
         * Creates a result indicating that no variation was found.
         */
        protected static Result noVariation(final AssemblyRegion targetRegion) {
            final Result result = new Result(targetRegion, Collections.emptyList(), Pair.of(targetRegion.getSpan(), null), null, null);
            result.leftFlankRegion = targetRegion;
            return result;
        }
    }

    /**
     * Returns a trimming result from which the variant-trimmed region and non-variant flanks can be recovered.
     *
     * @param region the genome location range to trim.
     * @param variants list of variants to trim to.  Those not overlapping with {@code region} are ignored.
     * @return never {@code null}.
     */
    public Result trim(final AssemblyRegion region, final SortedSet<VariantContext> variants) {
        if ( variants.isEmpty() ) {
            return Result.noVariation(region);
        }

        boolean foundNonSnp = false;
        SimpleInterval variantSpan = null;

        final List<VariantContext> variantsInRegion = variants.stream().filter(region::overlaps).collect(Collectors.toList());

        if ( variantsInRegion.isEmpty() ) {
            return Result.noVariation(region);
        } else if ( assemblyArgs.dontTrimActiveRegions) {
            return Result.noTrimming(region, variantsInRegion);
        }

        for ( final VariantContext vc : variantsInRegion ) {
                final SimpleInterval vcLoc = new SimpleInterval(vc);
                foundNonSnp = foundNonSnp || ! vc.isSNP();
                variantSpan = variantSpan == null ? vcLoc : variantSpan.spanWith(vcLoc);
        }

        final int padding = foundNonSnp ? assemblyArgs.indelPadding : assemblyArgs.snpPadding;

        final SimpleInterval maximumSpan = new SimpleInterval(region).expandWithinContig(usableExtension, sequenceDictionary);
        final SimpleInterval idealSpan = variantSpan.expandWithinContig(padding, sequenceDictionary);
        final SimpleInterval finalSpan = maximumSpan.intersect(idealSpan).mergeWithContiguous(variantSpan);
        final SimpleInterval callableSpan = variantSpan.intersect(region);

        final Pair<SimpleInterval, SimpleInterval> nonVariantFlanks = getFlanks(region, callableSpan);

        if ( debug ) {
            logger.info("events       : " + variantsInRegion);
            logger.info("region       : " + region);
            logger.info("callableSpan : " + callableSpan);
            logger.info("padding      : " + padding);
            logger.info("finalSpan    : " + finalSpan);
        }

        return new Result(region, variantsInRegion, nonVariantFlanks,finalSpan, variantSpan);
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