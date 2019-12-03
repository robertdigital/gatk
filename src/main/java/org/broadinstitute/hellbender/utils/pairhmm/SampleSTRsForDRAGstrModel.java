package org.broadinstitute.hellbender.utils.pairhmm;


import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.RuntimeIOException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UncheckedIOException;
import java.nio.file.Paths;
import java.util.*;

@CommandLineProgramProperties(
        programGroup = ReferenceProgramGroup.class,
        summary = "Determine the presence of STR in a reference sequence",
        oneLineSummary = "Determines the presence of STR in a reference sequence"
)
public class SampleSTRsForDRAGstrModel extends GATKTool {

    private static final Logger logger = LogManager.getLogger(SampleSTRsForDRAGstrModel.class);

    public static class DecimationTable {

        private static final int[][] DEFAULT_DECIMATION_MATRIX = new int[][] {
                {0}, // 0, 0, 0, 0, 0, 0, 0, 0 ...
                {0, 10, 10, 9, 8, 7, 5, 3, 1, 0},
                {0, 0, 9, 6, 3, 0}, // 0, 0, 0 ...
                {0, 0, 8, 4, 1, 0},
                {0, 0, 6, 0},
                {0, 0, 5, 0},
                {0, 0, 4, 0},
                {0}};

        public static final String NO_DECIMATION_STR = "NONE";

        public static final String DEFAULT_DECIMATION_STR = "DEFAULT";

        public static final DecimationTable DEFAULT = new DecimationTable(DEFAULT_DECIMATION_STR);

        public static final DecimationTable NONE = new DecimationTable(NO_DECIMATION_STR);



        private final long[][] decimationMask;

        private final long[][] counts;


        public DecimationTable(final String spec) {
            Utils.nonNull(spec);
            final int[][] decimation;
            if (spec.equalsIgnoreCase(NO_DECIMATION_STR)) {
                decimation = new int[][] {{0}};
            } else if (spec.equalsIgnoreCase(DEFAULT_DECIMATION_STR)) {
                decimation = DEFAULT_DECIMATION_MATRIX;
            } else {
                decimation = parseDecimationMatrixFromPath(spec);
            }
            decimationMask = calculateDecimationMask(decimation);
            counts = composeDecimationCounts(decimationMask);
        }

        private long[][] composeDecimationCounts(final long[][] decimationMask) {
            final long[][] result = new long[decimationMask.length][];
            for (int i = 0; i < result.length; i++) {
                result[i] = new long[decimationMask[i].length];
            }
            return result;
        }

        private static int[][] parseDecimationMatrixFromPath(String spec) {
            try (final BufferedReader reader = new BufferedReader(IOUtils.makeReaderMaybeGzipped(Paths.get(spec)))) {
                final String[][] values = reader.lines()
                        .filter(str -> !str.startsWith("#") && !str.trim().isEmpty())
                        .map(str -> Arrays.stream(str.split("\\s+"))
                                   .mapToDouble(Double::parseDouble)
                                   .toArray())
                        .toArray(String[][]::new);
                return parseDecimationMatrixValues(values, spec);
            } catch (final IOException ex) {
                throw new UserException.CouldNotReadInputFile(spec, ex);
            } catch (final NumberFormatException ex){
                throw new UserException.BadInput(String.format("input decimation file %s contains non-valid values: %s", spec, ex.getMessage()));
            }
        }

        private static int[][] parseDecimationMatrixValues(final String[][] values, final String path) {
            Utils.nonNull(values);
            if (values.length == 0) {
                logger.warn("Decimation matrix path provided does not seem to contain any values, we will proceed without any decimation");
                return new int[0][];
            } else {
                int totalValues = 0;
                final int[][] result = new int[values.length][];
                for (int i = 0; i < values.length; i++) {
                    final String[] row = values[i];
                    final int[] rowValues = new int[values.length];
                    for (int j = 0; j <  row.length; j++) {
                        final String str = row[j];
                        final int value;
                        try {
                            value = Integer.parseInt(str);
                        } catch (final NumberFormatException ex) {
                            throw badDecimationValueException(str, path, i, j, "not a valid double literal");
                        }
                        if (value < 0) {
                            throw badDecimationValueException(str, path, i, j, "negatives are not allowed");
                        } else if (Double.isNaN(value)) {
                            throw badDecimationValueException(str, path, i, j, "NaN are not allowed");
                        } else if (!Double.isFinite(value)) {
                            throw badDecimationValueException(str, path, i, j, "must be finite");
                        }
                        rowValues[j] = value;
                        totalValues++;
                    }
                    result[i] = rowValues;
                }
                if (totalValues == 0) {
                    throw new UserException.BadInput("the input decimation matrix does contain any values:" + path);
                }
                return result;
            }
        }

        private static RuntimeException badDecimationValueException(final String str, final String path, final int i, final int j,
                                                                    final String details) {
            throw new UserException.BadInput(String.format("bad decimation value found in %s for period and repeats (%d, %d) with string (%s)%s",
                    path, i, j, str, details == null || details.isEmpty()? "": ": " + details));
        }

        public static long[][] calculateDecimationMask(final int[][] decimationMatrix) {
            Utils.nonNull(decimationMatrix);
            final long[][] result = new long[decimationMatrix.length][];
            for (int i = 0; i < result.length; i++) {
                final int[] row = decimationMatrix[i];
                result[i] = new long[row.length];
                for (int j = 0; j < row.length; j++) {
                    result[i][j] = (1 << row[j]) - 1;
                }
            }
            return result;
        }

        public long mask(final int period, final int repeats) {
            final int p = period >= decimationMask.length ? decimationMask.length - 1 : period;
            final long[] masks = decimationMask[p];
            if (masks.length == 0) {
                return 0;
            } else if (repeats >= masks.length) {
                return masks[masks.length - 1];
            } else {
                return masks[repeats];
            }
        }

        public boolean decimate(final int seqNumber, final int bestPeriod, final long bestPeriodRepeats) {
            if (counts.length <= bestPeriod) {
                return true;
            } else {
                final long[] periodCounts = counts[bestPeriod];
                if (periodCounts.length == 0) {
                    return true;
                } else {
                    final int effectiveRepeatCount
                            = (int) (bestPeriodRepeats < periodCounts.length ? bestPeriodRepeats : periodCounts.length - 1);
                    final long count = periodCounts[effectiveRepeatCount]++;
                    final long left = count + seqNumber;
                    final long right = decimationMask[bestPeriod][effectiveRepeatCount];
                    return ((int) left & (int) right) == 0 && ((left >> 32) & (right >> 32)) == 0;
                }
            }
        }
    }

    @Argument(fullName="decimation", doc="decimation per perior and repeat. It can be \"DEFAULT\" to use the default values (default), " +
            " \"NONE\" to deactivate decimation (potentially resulting in a very large output file) or indicate the path to a file" +
            " that contains the decimation matrix.", optional = true)
    private DecimationTable decimationTable = DecimationTable.DEFAULT;

    @Argument(fullName="output", shortName = "O")
    private String outputPath;

    @Argument(fullName="max-period", doc="maximum STR period sampled", optional = true, minValue = 1, maxValue = 10)
    private int maxPeriod = 8;

    @Argument(fullName="max-repeats", doc="maximum STR repeat sampled", optional = true, minValue = 1, maxValue = 20)
    private int maxRepeat = 20;

    @Override
    public boolean requiresReference() {
        return true;
    }


    @Override
    public void traverse() {
        final SAMSequenceDictionary dictionary = getReferenceDictionary();
        final ReferenceDataSource referenceDataSource = directlyAccessEngineReferenceDataSource();
        try (final PrintWriter outputWriter = new PrintWriter(BucketUtils.createFile(outputPath))) {
            outputWriter.println("#chrom\tstart\tend\tperiod\trepeat\tunit");
            for (final SAMSequenceRecord sequence : dictionary.getSequences()) {
                traverse(sequence.getSequenceIndex() + 1, sequence, outputWriter, decimationTable);
            }

        } catch (final RuntimeIOException | UncheckedIOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputPath, ex);
        }
    }

    private void traverse(final int seqNumber, final SAMSequenceRecord sequence, final PrintWriter outputWriter, final DecimationTable decimationTable) {
        final String id = sequence.getSequenceName();
        final NucleotideSequence fullSequence = LazyLoadingReferenceNucleotideSequence.of(directlyAccessEngineReferenceDataSource(), sequence.getSequenceName(),10000);
        final int maxPeriod = this.maxPeriod;
        final int length = sequence.getSequenceLength();
        long pos = 1;
        while (pos <= length) {
           final byte base = fullSequence.byteAt(pos);
           long beg, end, cmp;
           for (beg = pos - 1; beg >= 1 && fullSequence.byteAt(beg) == base; beg--);
           beg++;
           for (end = pos + 1; end <= length && fullSequence.byteAt(end) == base; end++);
           end--;
           int bestPeriod = 1;
           long bestPeriodRepeats = end - beg + 1;
           long bestEnd = end;
           long bestBeg = beg;
           for (int period = 2; period <= maxPeriod; period++) {
               for (beg = length - pos > period ? pos - 1 : length - period,
                    cmp = beg + period; beg >= 1 && fullSequence.byteAt(cmp) == fullSequence.byteAt(beg); beg--, cmp--);
               beg++;
               for (end = pos >= period ? pos + 1 : period + 1, cmp = end - period; end <= length && fullSequence.byteAt(end) == fullSequence.byteAt(cmp); end++, cmp++);
               end--;
               final long strLength = end - beg + 1;
               final long repeats =  strLength / period;
               if (repeats > bestPeriodRepeats) {
                   bestPeriod = period;
                   bestPeriodRepeats = repeats;
                   bestEnd = end;
                   bestBeg = beg;
               }
           }
           emitOrDecimateSTR(seqNumber, id, bestBeg, bestBeg + bestPeriod * bestPeriodRepeats - 1, bestPeriod, bestPeriodRepeats, outputWriter, decimationTable, fullSequence);
           pos = bestEnd + 1;
        }
    }

    private void emitOrDecimateSTR(final int seqNumber, final String seqId, long bestBeg, long bestEnd, int bestPeriod, long bestPeriodRepeats, final PrintWriter outputWriter, final DecimationTable decimationTable, final NucleotideSequence sequence) {
        if (decimationTable.decimate(seqNumber, bestPeriod, bestPeriodRepeats)) {
            final String unit = sequence.subString(bestEnd - bestPeriod + 1, bestEnd);
            outputWriter.println(Utils.join("\t", seqId, bestBeg, bestEnd, bestPeriod, bestPeriodRepeats, unit));
        }
    }
}
