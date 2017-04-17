/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.seqfilter.utils;

import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.writers.FastaWriter;
import edu.msu.cme.rdp.seqfilter.SeqFilterResult;
import edu.msu.cme.rdp.seqfilter.filters.PrimerFilter;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

/**
 *
 * @author fishjord
 */
public class YankPrimers {

    private static final Options opts = new Options();

    static {
        opts.addOption("f", "forward", true, "Forward primers (comma seperated list)");
        opts.addOption("r", "reverse", true, "Reverse primers (comma seperated list)");
        opts.addOption("F", "forward-dist", true, "Maximum forward primer edit distance (default=2)");
        opts.addOption("R", "reverse-dist", true, "Maximum reverse primer edit distance (default=0)");
        opts.addOption("i", "in", true, "Input file");
        opts.addOption("o", "out", true, "Output file");
    }

    public static void yankPrimers(File inFile, File outFile, String[] fowardPrimers, int maxForwardDist, String[] reversePrimers, int maxReverseDist, boolean reversible) throws IOException {
        SequenceReader seqReader = new SequenceReader(inFile);
        PrimerFilter filter = new PrimerFilter(Arrays.asList(fowardPrimers), Arrays.asList(reversePrimers), maxForwardDist, maxReverseDist, reversible, false);
        Sequence seq;

        FastaWriter out = new FastaWriter(outFile);

        try {
            while ((seq = seqReader.readNextSequence()) != null) {
                SeqFilterResult result = filter.filterSequence(seq);
                if (!result.failed()) {
                    out.writeSeq(result.getResultSeq());
                } else {
                    System.err.println(result.getErrorMessage());
                }
            }
        } finally {
            out.close();
        }

    }

    public static void main(String[] args) throws IOException {
        String[] forwardPrimers = null;
        String[] reversePrimers = null;
        File inFile = null;
        File outFile = null;
        int maxForward = 2;
        int maxReverse = 0;

        try {
            CommandLine line = new PosixParser().parse(opts, args);

            if(line.hasOption("forward")) {
                forwardPrimers = line.getOptionValue("forward").split(",");
            } else {
                throw new Exception("Must supply at least one forward primer");
            }

            if(line.hasOption("in")) {
                inFile = new File(line.getOptionValue("in"));
            } else {
                throw new Exception("Input file required");
            }

            if(line.hasOption("out")) {
                outFile = new File(line.getOptionValue("out"));
            } else {
                throw new Exception("Output file required");
            }

            if(line.hasOption("reverse")) {
                reversePrimers = line.getOptionValue("reverse").split(",");
            }

            if(line.hasOption("reverse-dist")) {
                maxForward = new Integer(line.getOptionValue("reverse-dist"));
            }

            if(line.hasOption("forward-dist")) {
                maxReverse = new Integer(line.getOptionValue("forward-dist"));
            }

        } catch (Exception e) {
            System.err.println("ERROR: " + e.getMessage());
            return;
        }

        YankPrimers.yankPrimers(inFile, outFile, forwardPrimers, maxForward, reversePrimers, maxReverse, false);
    }
}
