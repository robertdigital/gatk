package org.broadinstitute.hellbender.tools.spark.longread;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.LongReadAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;

/**
 * Extracts alignment records that have strange CIGARs.
 *
 * <p>
 *     Currently we define "strange" as:
 *     <ul>
 *         <li>CIGAR having clipping neighboring non-alignment blocks (alignments will have a tag "XS:i:1")</li>
 *         <li>CIGAR having an insertion/deletion operation neighboring a deletion/insertion operation
 *             (alignments will have a tag "XG:Z:section_of_strange_cigar")</li>
 *         <li>CIGAR having a small alignment block sandwiched between two large gaps,
 *             i.e. [0-9]+(I|D)[0-9]+M[0-9]+(I|D) where the middle alignment block is much smaller than the neighboring gaps
 *             (alignments will be annotated with a tag "XC:Z:section_of_strange_cigar")</li>
 *     </ul>
 * </p>
 *
 * <p>
 *     Note that we currently don't extract alignment records that share the same read name as the strange ones,
 *     i.e. their mates (if applicable), their SA/XA records, unless they themselves are strange too.
 *     If you do want those, we recommend you taking a look at
 *     {@link org.broadinstitute.hellbender.tools.spark.sv.utils.ExtractOriginalAlignmentRecordsByNameSpark}.
 * </p>
 */
@DocumentedFeature
@CommandLineProgramProperties(summary = "Extracts alignment records that seem strange. ",
        oneLineSummary = "Extracts alignment records that seem strange",
        programGroup = LongReadAnalysisProgramGroup.class)
@BetaFeature
public class ExtractAlignmentsWithStrangeCigarSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(ExtractAlignmentsWithStrangeCigarSpark.class);

    @Override
    public boolean requiresReads()
    {
        return true;
    }

    public static final String NEIGHBOR_INSDEL_OP = "XG";
    public static final String NEIGHBOR_CLIPGAP_OP = "XS";
    public static final String SEMI_NEIGHBOR_INSDEL_OP = "XC";

    @Argument(doc = "for CIGARs that have closely neighbored gaps separated only by a small alignment block, " +
            "the triggering threshold below which fraction of the alignment size to the size of the smaller gaps",
            shortName = "f", fullName = "fraction",
            optional = true)
    private Double fraction = 0.1;

    @Argument(doc = "uri for the output file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String output;


    @Override
    protected void runTool( final JavaSparkContext ctx ) {

        final DelNeighboringIns catchDelNeighboringIns = new DelNeighboringIns();
        final NonAlnNeighboringClip catchNonAlnNeighboringClip = new NonAlnNeighboringClip();
        final CloselyNeighboringGaps closelyNeighboringGaps = new CloselyNeighboringGaps(fraction);

        final JavaRDD<GATKRead> strangeCigaredReads =
                getReads()
                        .filter(r -> !r.isUnmapped())
                        .map(read -> {
                            GATKRead copy = read.copy();
                           if (catchNonAlnNeighboringClip.test(read)) {
                                copy.setAttribute(NEIGHBOR_CLIPGAP_OP, 1);
                            } else if (catchDelNeighboringIns.test(read)) {
                               copy.setAttribute(NEIGHBOR_INSDEL_OP,
                                                 catchDelNeighboringIns.getOffendingBlocks(read));
                           } else if (closelyNeighboringGaps.test(read)) {
                                copy.setAttribute(SEMI_NEIGHBOR_INSDEL_OP,
                                                  closelyNeighboringGaps.getOffendingBlocks(read));
                            }
                            return copy;
                        })
                        .filter(read ->
                            read.hasAttribute(NEIGHBOR_INSDEL_OP)
                                    || read.hasAttribute(NEIGHBOR_CLIPGAP_OP)
                                    || read.hasAttribute(SEMI_NEIGHBOR_INSDEL_OP)
                        );
        writeReads(ctx, output, strangeCigaredReads, getHeaderForReads(), true);

        final long countOfNeighboringIndel = strangeCigaredReads.filter(read -> read.hasAttribute(NEIGHBOR_INSDEL_OP)).count();
        final long countOfClippingNeighboringIndel = strangeCigaredReads.filter(read -> read.hasAttribute(NEIGHBOR_CLIPGAP_OP)).count();
        final long countOfCloselyNeighboredIndel = strangeCigaredReads.filter(read -> read.hasAttribute(SEMI_NEIGHBOR_INSDEL_OP)).count();
        String message  = String.format("The input alignments has %d records with neighboring indel operations\n", countOfNeighboringIndel);
               message += String.format("                         %d records with clipping neighboring an indel\n", countOfClippingNeighboringIndel);
               message += String.format("                         %d records with small alignment sandwiched between large indels", countOfCloselyNeighboredIndel);
        localLogger.info(message);
    }

    /**
     * A read filter designed to catch (NOT filter out) reads that have CIGAR of the form
     * "...[0-9]+I[0-9]+D..." or "...[0-9]+D[0-9]+I...".
     *
     * Note: read is assumed to be mapped. Will throw if not.
     */
    private static final class DelNeighboringIns extends ReadFilter {
        private static final long serialVersionUID = 1L;

        @Override
        public boolean test( GATKRead read ) {
            final List<CigarElement> cigarElements = read.getCigarElements();
            CigarOperator precedingOp = cigarElements.get(0).getOperator();
            for (int i = 1; i < cigarElements.size(); ++i) {
                final CigarOperator currentOp = cigarElements.get(i).getOperator();
                if (precedingOp.isIndel() && currentOp.isIndel() &&
                        precedingOp != currentOp) {
                    return true;
                }
                precedingOp = currentOp;
            }
            return false;
        }


        public String getOffendingBlocks(final GATKRead read) {
            final List<CigarElement> cigarElements = read.getCigarElements();

            final StringBuilder attribute = new StringBuilder();

            CigarElement precedingElement = cigarElements.get(0);
            for (int i = 1; i < cigarElements.size(); ++i) {
                final CigarElement currentElement = cigarElements.get(i);
                if (precedingElement.getOperator().isIndel() && currentElement.getOperator().isIndel() &&
                        precedingElement.getOperator() != currentElement.getOperator()) {
                    attribute.append(precedingElement.toString());
                    attribute.append(currentElement.toString());
                    attribute.append(",");
                }
                precedingElement = currentElement;
            }


            return StringUtils.removeEnd(attribute.toString(), ",");
        }
    }

    /**
     * A read filter designed to catch (NOT filter out) reads that have CIGAR of the form
     * "([0-9]+H)?[0-9]+S[0-9]+(I|D)..." or "...[0-9]+(I|D)[0-9]+S([0-9]+H)?".
     *
     * Note: read is assumed to be mapped. Will throw if not.
     */
    private static final class NonAlnNeighboringClip extends ReadFilter {
        private static final long serialVersionUID = 1L;

        @Override
        public boolean test( GATKRead read ) {
            final List<CigarElement> cigarElements = read.getCigarElements();

            // front
            CigarOperator terminalOP = cigarElements.get(0).getOperator();
            if (terminalOP == CigarOperator.H) {
                CigarOperator operator = cigarElements.get(1).getOperator();
                if (operator == CigarOperator.S) {
                    if (cigarElements.get(2).getOperator().isIndel())
                        return true;
                } else if (operator.isIndel()) {
                    return true;
                }
            } else if (terminalOP == CigarOperator.S) {
                if (cigarElements.get(1).getOperator().isIndel())
                    return true;
            }

            // end
            final int n = cigarElements.size();
            terminalOP = cigarElements.get(n-1).getOperator();
            if (terminalOP == CigarOperator.H) {
                CigarOperator operator = cigarElements.get(n-2).getOperator();
                if (operator == CigarOperator.S) {
                    if (cigarElements.get(n-3).getOperator().isIndel())
                        return true;
                } else if (operator.isIndel()) {
                    return true;
                }
            } else if (terminalOP == CigarOperator.S) {
                if (cigarElements.get(n-2).getOperator().isIndel())
                    return true;
            }

            return false;
        }
    }

    /**
     * A read filter designed to catch (NOT filter out) reads that have CIGAR of the form
     * "...[0-9]+I[0-9]+M[0-9]+D..." or "...[0-9]+D[0-9]+M[0-9]+I...",
     * where the alignment block is shorter than both gap lengths.
     *
     * See {@link #fraction} for controlling how small the alignment needs to be for the filter to be triggered.
     *
     * Note: read is assumed to be mapped. Will throw if not.
     */
    private static final class CloselyNeighboringGaps extends ReadFilter {
        private static final long serialVersionUID = 1L;

        private final Double fraction;

        CloselyNeighboringGaps(final Double fraction) {
            this.fraction = fraction;
        }

        @Override
        public boolean test( GATKRead read ) {
            final List<CigarElement> cigarElements = read.getCigarElements();


            for (int i = 1; i < cigarElements.size() - 1; ++i) {
                final CigarElement precedingElement = cigarElements.get(i - 1);
                final CigarElement currentElement   = cigarElements.get(i);
                final CigarElement nextElement      = cigarElements.get(i + 1);

                if (precedingElement.getOperator().isIndel()
                        && nextElement.getOperator().isIndel()
                        && currentElement.getOperator().isAlignment()) {
                    final double threshold = Math.min(precedingElement.getLength(), nextElement.getLength()) * fraction;

                    if( currentElement.getLength() < 3 || currentElement.getLength() < threshold ) {
                        return true;
                    }
                }
            }
            return false;
        }

        public String getOffendingBlocks(final GATKRead read) {

            final List<CigarElement> cigarElements = read.getCigarElements();

            final StringBuilder attribute = new StringBuilder();
            for (int i = 1; i < cigarElements.size() - 1; ++i) {
                final CigarElement precedingElement = cigarElements.get(i - 1);
                final CigarElement currentElement   = cigarElements.get(i);
                final CigarElement nextElement      = cigarElements.get(i + 1);

                if (precedingElement.getOperator().isIndel()
                        && nextElement.getOperator().isIndel()
                        && currentElement.getOperator().isAlignment()) {
                    final double threshold = Math.min(precedingElement.getLength(), nextElement.getLength()) * fraction;

                    if( currentElement.getLength() < 3 || currentElement.getLength() < threshold ) {
                        attribute.append(precedingElement.toString());
                        attribute.append(currentElement.toString());
                        attribute.append(nextElement.toString());
                        attribute.append(",");
                    }
                }
            }

            return StringUtils.removeEnd(attribute.toString(), ",");
        }
    }
}
