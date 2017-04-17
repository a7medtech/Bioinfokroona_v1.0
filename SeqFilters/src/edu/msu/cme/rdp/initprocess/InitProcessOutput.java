/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.initprocess;

import edu.msu.cme.rdp.readseq.QSequence;
import edu.msu.cme.rdp.readseq.writers.FastaWriter;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.core.FastqCore;
import edu.msu.cme.rdp.readseq.writers.FastqWriter;
import edu.msu.cme.rdp.readseq.writers.SequenceWriter;
import edu.msu.cme.rdp.seqfilter.SeqFilter;
import edu.msu.cme.rdp.seqfilter.SeqFilterResult;
import edu.msu.cme.rdp.seqfilter.SeqFilteringResult;
import edu.msu.cme.rdp.seqfilter.output.SeqFilterOutput;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.geom.RectangularShape;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import org.apache.commons.lang.StringUtils;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.xy.XYBarPainter;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleEdge;

/**
 *
 * @author fishjord
 */
public class InitProcessOutput implements SeqFilterOutput {

    private FastqWriter fastqWriter;
    private FastaWriter fastaWriter;
    private int writtenSequences = 0;
    private Map<Integer, Integer> lengthMap = new HashMap();
    private Set<String> seqsMissingQualSeq = new HashSet();
    private Map<Integer, Integer> qualityMap = new LinkedHashMap();
    private Map<Integer, Integer> positionCounts = new HashMap();
    private static final DecimalFormat format = new DecimalFormat("####.##");

    public InitProcessOutput(File fastqOut, File fastaOut) throws IOException {
        this.fastaWriter = new FastaWriter(fastaOut);
        this.fastqWriter = new FastqWriter(fastqOut, FastqCore.Phred33QualFunction);
    }

    public void appendSequence(Sequence seq) {

        if (seq instanceof QSequence) {
            byte[] qual = ((QSequence) seq).getQuality();
            for (int index = 0; index < qual.length; index++) {
                byte q = qual[index];

                int positionCount = 0;
                int totalQual = 0;
                if (qualityMap.containsKey(index)) {
                    totalQual = qualityMap.get(index);
                    positionCount = positionCounts.get(index);
                }
                qualityMap.put(index, q + totalQual);
                positionCounts.put(index, positionCount + 1);
            }
            try {
                fastqWriter.writeSeq(seq);
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        } else {
            fastaWriter.writeSeq(seq);
        }

        int seqLength = seq.getSeqString().length();

        if (!lengthMap.containsKey(seqLength)) {
            lengthMap.put(seqLength, 0);
        }

        lengthMap.put(seqLength, lengthMap.get(seqLength) + 1);
        writtenSequences++;
    }

    public void close() throws IOException {
        fastaWriter.close();
        fastqWriter.close();
    }

    private void writeLengthStats(File lengthStatsFile, File lengthChartFile, int avgLength, int filteredSeqs) throws IOException {

        XYSeries histoData = new XYSeries("Sequence Length");
        PrintStream out = new PrintStream(lengthStatsFile);

        out.println("Sequence Length\tCount of Sequences with Length");
        for (int seqLength : lengthMap.keySet()) {
            out.println(seqLength + "\t" + lengthMap.get(seqLength));
            histoData.add(seqLength, lengthMap.get(seqLength));
        }

        out.close();

        JFreeChart chart = ChartFactory.createXYBarChart(
                "Initial Process Length Histogram", // chart title
                "Length in bases", // domain axis label
                false,
                "Number of sequences", // range axis label
                new XYSeriesCollection(histoData), // data
                PlotOrientation.VERTICAL, // orientation
                false, // include legend
                false, // tooltips?
                false // URLs?
                );

        chart.setSubtitles(Arrays.asList(new TextTitle("Avg Length: " + avgLength), new TextTitle("Total Seqs: " + (writtenSequences + filteredSeqs) + " After Trimming: " + writtenSequences)));

        XYBarRenderer renderer = (XYBarRenderer) chart.getXYPlot().getRenderer();
        renderer.setShadowVisible(false);
        renderer.setBarPainter(new XYBarPainter() {
            public void paintBar(Graphics2D arg0, XYBarRenderer arg1, int arg2, int arg3, RectangularShape arg4, RectangleEdge arg5) {
                Rectangle r = arg4.getBounds();
                arg0.setPaint(arg1.getItemPaint(arg2, arg3));
                arg0.fillRect((int) r.getX(), (int) r.getY(), (int) r.getWidth(), (int) r.getHeight());
            }

            public void paintBarShadow(Graphics2D arg0, XYBarRenderer arg1, int arg2, int arg3, RectangularShape arg4, RectangleEdge arg5, boolean arg6) {
                throw new UnsupportedOperationException("Not supported yet.");
            }
        });

        ChartUtilities.saveChartAsPNG(lengthChartFile, chart, 640, 640);
    }

    private void writeQualityStats(File qualStatsFile, File qualStatsChart, int avgLength, int filteredSeqs) throws IOException {

        XYSeries qualityScatterData = new XYSeries("Sequence Length");
        PrintStream out = new PrintStream(qualStatsFile);

        out.println("Position\tAverage Quality");
        for (int position : qualityMap.keySet()) {
            double avgQual = (double) qualityMap.get(position) / positionCounts.get(position);
            out.println((position + 1) + "\t" + avgQual);
            qualityScatterData.add(position + 1, avgQual);
        }

        out.close();

        JFreeChart chart = ChartFactory.createXYLineChart(
                "Average quality vs Position", // chart title
                "Position", // domain axis label
                "Average Quality", // range axis label
                new XYSeriesCollection(qualityScatterData), // data
                PlotOrientation.VERTICAL, // orientation
                false, // include legend
                false, // tooltips?
                false // URLs?
                );

        chart.setSubtitles(Arrays.asList(new TextTitle("Avg Length: " + avgLength), new TextTitle("Total Seqs: " + (writtenSequences + filteredSeqs) + " After Trimming: " + writtenSequences)));

        ChartUtilities.saveChartAsPNG(qualStatsChart, chart, 640, 640);
    }

    public void writeStats(File filteredSeqsFile, File summaryFile, File lengthStatsFile, File lengthChartFile, File qualStatsFile, File qualStatsChart, String tag, SeqFilteringResult seqFilterResult) throws IOException {

        int filteredSeqs = 0;
        PrintStream out = new PrintStream(filteredSeqsFile);
        Map<SeqFilter, Map<String, SeqFilterResult>> filteredSeqsMap = seqFilterResult.getFilteredSeqsMap();

        out.println("Sequence Id\tFilter Name\tReason for Filtering");
        for (SeqFilter filter : filteredSeqsMap.keySet()) {
            for (String seqid : filteredSeqsMap.get(filter).keySet()) {
                out.println(seqid + "\t" + filter.getName() + "\t" + filteredSeqsMap.get(filter).get(seqid).getErrorMessage());
                filteredSeqs++;
            }

            for (String seqid : seqsMissingQualSeq) {
                out.println(seqid + "\tNo Quality Sequence\tNot filtered, warning only");
            }
        }

        out.close();

        if (writtenSequences == 0) {
            return;
        }

        int avgLength = 0;
        double stdDev = 0;


        for (int length : lengthMap.keySet()) {
            avgLength += length * lengthMap.get(length);
        }
        avgLength /= writtenSequences;

        for (int length : lengthMap.keySet()) {
            stdDev += Math.pow(length - avgLength, 2) * lengthMap.get(length);
        }
        stdDev = (writtenSequences == 1) ? Double.NaN : Math.sqrt(stdDev / (writtenSequences - 1));

        out = new PrintStream(summaryFile);
        out.println("Tag\t" + tag);
        out.println("Total sequences\t" + (writtenSequences + filteredSeqs));
        out.println("Sequences after trimming\t" + writtenSequences);
        out.println("Average sequence length after trimming\t" + avgLength);
        out.println("Standard Deviation of sequence length\t" + format.format(stdDev));
        out.println();
        out.println("Filtered sequence summary");
        for (SeqFilter filter : filteredSeqsMap.keySet()) {
            out.println(filter.getName() + "\t" + filteredSeqsMap.get(filter).size());
        }
        out.close();

        writeLengthStats(lengthStatsFile, lengthChartFile, avgLength, filteredSeqs);

        if (!qualityMap.isEmpty()) {
            writeQualityStats(qualStatsFile, qualStatsChart, avgLength, filteredSeqs);
        }
    }
}
