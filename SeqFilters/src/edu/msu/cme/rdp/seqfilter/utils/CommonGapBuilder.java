package edu.msu.cme.rdp.seqfilter.utils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import edu.msu.cme.rdp.readseq.MaskSequenceNotFoundException;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.IndexedSeqReader;

public class CommonGapBuilder {

    public static final String FASTA_REFERENCE_KEY = "#=GC_RF";

    public static Set<Integer> getGapSet(File file, String maskName) throws IOException, MaskSequenceNotFoundException {
        ArrayList<File> files = new ArrayList<File>();
        files.add(file);
        return getGapSet(files, maskName);
    }

    public static Set<Integer> getGapSet(File file) throws IOException, MaskSequenceNotFoundException {
        return getGapSet(file, FASTA_REFERENCE_KEY);
    }

    public static Set<Integer> getGapSet(List<File> files) throws IOException, MaskSequenceNotFoundException {
        return getGapSet(files, FASTA_REFERENCE_KEY);
    }

    public static Set<Integer> getGapSet(List<File> files, String maskName) throws IOException, MaskSequenceNotFoundException {
        Set<Integer> gapset = null;
        for (File f : files) {
            IndexedSeqReader reader = new IndexedSeqReader(f, maskName);
            Set<Integer> gaps = buildGapSet(reader);

            if (gapset != null) {
                gapset.retainAll(gaps);
            } else {
                gapset = gaps;
            }
        }
        return gapset;
    }

    /**
     *
     * @param seqList (a List of Sequence Object)
     * @return the set of common gap positions, including the model positions
     */
    public static Set<Integer> buildGapSet(List<Sequence> seqList) {
        Set<Integer> gapset = null;
        Iterator<Sequence> it = seqList.iterator();

        while (it.hasNext()) {
            Sequence seq = it.next();
            Set<Integer> gaps = idGapCols(seq.getSeqString());

            if (gapset != null) {
                gapset.retainAll(gaps);
            } else {
                gapset = gaps;
            }
        }
        return gapset;

    }

    /**
     *
     * @param IndexedSeqReader<T> reader
     * @return the set of common gap positions, including the model positions
     * @throws IOException
     */
    public static Set<Integer> buildGapSet(IndexedSeqReader reader) throws IOException {
        Set<Integer> gapset = null;
        for (String seqid : reader.getSeqIds()) {
            String seqstring = reader.readSeq(seqid).getSeqString();
            Set<Integer> gaps = idGapCols(seqstring);

            if (gapset != null) {
                gapset.retainAll(gaps);
            } else {
                gapset = gaps;
            }
        }
        return gapset;

    }

    private static Set<Integer> idGapCols(String seqstring) {
        Set<Integer> gaps = new HashSet<Integer>();
        char[] bases = seqstring.toCharArray();

        for (int index = 0; index < bases.length; index++) {
            char c = bases[index];
            if (c == '-' || c == '.') {
                gaps.add(index);
            }
        }

        return gaps;
    }

    public static void main(String[] args) throws IOException, MaskSequenceNotFoundException {
        ArrayList<File> files = new ArrayList<File>();
        for (int i = 0; i < args.length; i++) {
            files.add(new File(args[i]));
        }
        Set<Integer> gapset = CommonGapBuilder.getGapSet(files);
        CommonGapRemover remover = new CommonGapRemover(gapset);
        for (File f : files) {
            IndexedSeqReader reader = new IndexedSeqReader(f, CommonGapBuilder.FASTA_REFERENCE_KEY);
            for (String seqid : reader.getSeqIds()) {
                String noGapseq = remover.removeCommonGaps(reader.readSeq(seqid));
                System.err.println(">" + seqid + "\n" + noGapseq);
            }
        }

    }
}
