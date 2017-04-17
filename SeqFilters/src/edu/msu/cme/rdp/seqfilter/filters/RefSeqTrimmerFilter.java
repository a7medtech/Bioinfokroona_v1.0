/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.seqfilter.filters;

import edu.msu.cme.rdp.readseq.QSequence;
import edu.msu.cme.rdp.readseq.readers.QSeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.utils.IUBUtilities;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import edu.msu.cme.rdp.readseq.writers.FastaWriter;
import edu.msu.cme.rdp.seqfilter.SeqFilter;
import edu.msu.cme.rdp.seqfilter.SeqFilterResult;
import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author fishjord
 */
public class RefSeqTrimmerFilter implements SeqFilter {

    private int expectedSeqLength;
    private int trimStart;
    private int trimEnd;
    private Map<Integer, Integer> indexToRefPos = new HashMap();

    /**
     * Trim points are inclusive, starting with base 1 (NOT 0)
     *
     * @param refSeq
     * @param refTrimStart
     * @param refTrimEnd
     */
    public RefSeqTrimmerFilter(Sequence refSeq, int refTrimStart, int refTrimEnd) {
        this(refSeq, null, refTrimStart, refTrimEnd);
    }

    public RefSeqTrimmerFilter(Sequence refSeq, Sequence maskSeq, int refTrimStart, int refTrimEnd) {
        char[] refSeqStr = refSeq.getSeqString().toCharArray();
        char[] maskSeqStr = null;

        expectedSeqLength = refSeqStr.length;

        if (maskSeq != null) {
            maskSeqStr = maskSeq.getSeqString().toCharArray();
            if (refSeqStr.length != maskSeqStr.length) {
                throw new IllegalArgumentException("Ref sequence length [" + refSeqStr.length + "] differs from mask seq length [" + maskSeqStr.length + "]");
            }
        }

        int refIndex = 0;
        for (int index = 0; index < refSeqStr.length; index++) {
            if (!IUBUtilities.isGap(refSeqStr[index])) {
                refIndex++;
                indexToRefPos.put(index, refIndex);

                if (refIndex == refTrimStart) {
                    if (maskSeqStr != null && (maskSeqStr[index] == '0' || maskSeqStr[index] == '.')) {
                        System.err.println("Warning: ref sequence position " + refTrimStart + " is a non model position");
                    }

                    trimStart = index;
                }

                if (refIndex == refTrimEnd) {
                    if (maskSeqStr != null && (maskSeqStr[index] == '0' || maskSeqStr[index] == '.')) {
                        System.err.println("Warning: ref sequence position " + refTrimEnd + " is a non model position");
                    }

                    //We want to be inclusive...teehee
                    trimEnd = index + 1;
                }
            }
        }
    }

    public SeqFilterResult filterSequence(Sequence seq) {
        String newSeqString = seq.getSeqString();
        byte[] newQual = null;

        if (seq instanceof QSequence) {
            newQual = ((QSequence) seq).getQuality();
        }

        if (seq.getSeqString().length() != expectedSeqLength) {
            throw new IllegalArgumentException(seq.getSeqName() + "'s length [" + seq.getSeqString().length() + "] doesn't match expected [" + expectedSeqLength + "]");
        }

        newSeqString = newSeqString.substring(trimStart, trimEnd);
        String unalignedSeqStr = SeqUtils.getUnalignedSeqString(newSeqString);

        int seqStart = -1, seqEnd = 0;
        char[] bases = seq.getSeqString().toCharArray();
        for(int index = 0;index < bases.length;index++) {
            if(!IUBUtilities.isGap(bases[index])) {
                if(seqStart == -1) {
                    seqStart = index;
                }

                seqEnd = index;
            }
        }

        if (unalignedSeqStr.isEmpty() || seqStart > trimStart || seqEnd < trimEnd) {
            Integer refStart = indexToRefPos.get(seqStart);
            Integer refEnd = indexToRefPos.get(seqEnd);
            return new SeqFilterResult("Sequence didn't cover region: " + refStart + ", " + refEnd);
        }

        Sequence ret = null;
        if (newQual != null) {
            int qualStart = 0;
            for (int index = trimStart; index >= 0; index--) {
                if (!IUBUtilities.isGap(bases[index])) {
                    qualStart++;
                }
            }

            newQual = Arrays.copyOfRange(newQual, qualStart, qualStart + unalignedSeqStr.length());
            ret = new QSequence(seq.getSeqName(), seq.getDesc(), newSeqString, newQual);
        } else {
            ret = new Sequence(seq.getSeqName(), seq.getDesc(), newSeqString);
        }

        return new SeqFilterResult(seq);
    }

    public String getName() {
        return "Reference sequence trimmer [" + indexToRefPos.get(trimStart) + "-" + indexToRefPos.get(trimEnd) + "]";
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 4 && args.length != 5) {
            System.err.println("USAGE: RefSeqTrimmer <refseqid> <refseq_start> <refseq_end> <seqfile> [qual_file]");
            return;
        }

        String refseqid = args[0];
        int refStart = Integer.valueOf(args[1]);
        int refEnd = Integer.valueOf(args[2]);

        File seqFile = new File(args[3]);
        File qualFile = null;

        String maskSeqId = "#=GC_RF";
        FastaWriter seqOut = new FastaWriter("trimmed_" + seqFile.getName());
        FastaWriter qualOut = null;

        if (args.length == 5) {
            qualFile = new File(args[4]);
            qualOut = new FastaWriter("trimmed_" + qualFile.getName());
        }

        int inCount = 0;
        int outCount = 0;
        if (qualFile != null) {
            QSeqReader reader = new QSeqReader(seqFile, qualFile, false);
            Sequence maskSeq = maskSeq = reader.readSeq(maskSeqId);
            List<String> seqids = reader.getSeqIds();

            Sequence refSeq = reader.readSeq(refseqid);

            RefSeqTrimmerFilter seqTrimmer = new RefSeqTrimmerFilter(refSeq, maskSeq, refStart, refEnd);

            for (String seqid : seqids) {
                if (seqid.startsWith("#")) {
                    continue;
                }
                inCount++;

                Sequence seq = reader.readSeq(seqid);
                QSequence trimmedSeq = (QSequence) seqTrimmer.filterSequence(seq).getResultSeq();

                if (trimmedSeq == null) {
                    continue;
                }
                outCount++;

                seqOut.writeSeq(trimmedSeq);

                StringBuilder qualSeq = new StringBuilder();
                for (byte b : trimmedSeq.getQuality()) {
                    qualSeq.append(b).append("  ");
                }

                qualOut.writeSeq(trimmedSeq.getSeqName(), trimmedSeq.getDesc(), qualSeq.toString());
            }

        } else {
            throw new IllegalStateException("I'm lazy and haven't implemented this yet");
        }

        System.out.println("Read in " + inCount + " sequences");
        System.out.println("Wrote out " + outCount + " trimmed sequences");

        seqOut.close();
        if (qualOut == null) {
            qualOut.close();
        }
    }
}
