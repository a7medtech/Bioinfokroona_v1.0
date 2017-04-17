package edu.msu.cme.pyro.derep;

/*
 * Ribosomal Database Project II
 * Copyright 2009 Michigan State University Board of Trustees
 */
import edu.msu.cme.rdp.readseq.writers.FastaWriter;
import edu.msu.cme.rdp.readseq.MaskSequenceNotFoundException;
import edu.msu.cme.rdp.readseq.QSequence;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.SequenceParsingException;
import edu.msu.cme.rdp.readseq.readers.IndexedSeqReader;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.codec.digest.DigestUtils;

/**
 *
 * @author farrisry
 */
public class Dereplicator {

    Map<String, List<DerepSeq>> hashToSeq = new HashMap<String, List<DerepSeq>>();
    List<DerepSeq> seqs = new ArrayList<DerepSeq>();

    static class DerepSeq {

        List<String> ids = new ArrayList<String>();
        Sequence seq;

        public DerepSeq(Sequence seq) {
            this.seq = seq;
        }
    }

    public void addSeq(Sequence seq) {
        String bases = seq.getSeqString().toUpperCase();

        String shaHash = DigestUtils.shaHex(bases);

        List<DerepSeq> hashSeqs = hashToSeq.get(shaHash);
        if (hashSeqs == null) {
            hashSeqs = new ArrayList<DerepSeq>();
            hashToSeq.put(shaHash, hashSeqs);
        }

        DerepSeq exemplar = null;
        for (DerepSeq s : hashSeqs) {
            if (s.seq.getSeqString().equals(seq.getSeqString())) {
                exemplar = s;
                break;
            }
        }

        if (exemplar == null) {
            exemplar = new DerepSeq(seq);
            hashSeqs.add(exemplar);
            seqs.add(exemplar);
        }

        exemplar.ids.add(seq.getSeqName());
    }

    public Set<DerepSeq> getUniqueSeqs() {
        Set<DerepSeq> seqToIds = new HashSet<DerepSeq>();
        for (DerepSeq seq : seqs) {
            seqToIds.add(seq);
        }
        return seqToIds;
    }

    public Set<DerepSeq> getSortedUniqueSeqs() {
        Set<DerepSeq> seqToIds = new TreeSet<DerepSeq>(new SeqComparator());
        for (DerepSeq seq : seqs) {
            seqToIds.add(seq);
        }
        return seqToIds;
    }
    
    public static class SeqComparator implements Comparator {
        public int compare(Object a, Object b){
            if ( ((DerepSeq)a).ids.size() < ((DerepSeq)b).ids.size()){
                return 1;
            } else {
                return -1;
            }
        }
    }
    
    private static Set<Integer> idGapCols(Dereplicator derep) {
        Set<Integer> ret = null;

        for (DerepSeq seq : derep.getUniqueSeqs()) {
            Set<Integer> gaps = new HashSet();
            char[] bases = seq.seq.getSeqString().toCharArray();

            for (int index = 0; index < bases.length; index++) {
                char c = bases[index];
                if (c == '-' || c == '.') {
                    gaps.add(index);
                }
            }

            if (ret == null) {
                ret = gaps;
            } else {
                ret.retainAll(gaps);
            }
        }

        return ret;
    }

    private static String removeCommonGaps(Set<Integer> gaps, String seq) {
        StringBuilder ret = new StringBuilder();

        char[] bases = seq.toCharArray();
        for (int index = 0; index < bases.length; index++) {
            if (!gaps.contains(index)) {
                ret.append(bases[index]);
            }
        }

        return ret.toString();
    }

    private static enum DerepMode { unaligned, aligned, model_only, formatted_model };

    public static void main(String... args) throws IOException, SequenceParsingException, MaskSequenceNotFoundException {
        Options options = new Options();
        options.addOption("u", "unaligned", false, "Dereplicate unaligned sequences");
        options.addOption("a", "aligned", false, "Dereplicate aligned sequences");
        options.addOption("m", "model-only", true, "Dereplicate aligned sequences using mask sequence");
        options.addOption("f", "formatted", false, "Dereplicate formated (uppercase/- = comparable, lowercase/. = non-comparable) aligned sequences");
        options.addOption("g", "keep-common-gaps", false, "Don't remove common gaps in output sequences");

        options.addOption("q", "qual-out", true, "Write quality sequences to this file");
        options.addOption("o", "out", true, "Write sequences to this file");
        
        options.addOption("s", "sorted", false, "Sort sequence by number of members represented");

        FastaWriter qualOut = null;
        FastaWriter seqWriter = new FastaWriter(System.out);
        DerepMode mode = null;
        String maskId = null;
        boolean keepCommonGaps;
        boolean sorted = false;

        Dereplicator derep = new Dereplicator();
        SampleMapping<String> sampleMapping = new SampleMapping<String>();

        try {
            CommandLine line = new PosixParser().parse(options, args);

            if (line.hasOption("unaligned")) {
                mode = DerepMode.unaligned;
            } 
            if (line.hasOption("aligned")) {
                if ( mode != null){
                    throw new Exception("Already specified mode \"" + mode.toString() + "\", can not take option " + "--aligned");
                }
                mode = DerepMode.aligned;
            } 
            if (line.hasOption("model-only")) {
                if ( mode != null){
                    throw new Exception("Already specified mode \""  + mode.toString() + "\", can not take option " + "--model-only");
                }
                mode = DerepMode.model_only;
                maskId = line.getOptionValue("model-only");
                System.err.println("Using " + maskId + " as mask sequence");
            } 
            if(line.hasOption("formatted")) {
                if ( mode != null){
                    throw new Exception("Already specified mode \""  + mode.toString() + "\", can not take option " + "--formatted");
                }
                mode = DerepMode.formatted_model;
            } 
            if ( mode == null) {
                System.err.println("Warning, no derep mode specified, using unaligned");
                mode = DerepMode.unaligned;
            }

            if (line.hasOption("qual-out")) {
                qualOut = new FastaWriter(line.getOptionValue("qual-out"));
            }

            if (line.hasOption("out")) {
                seqWriter = new FastaWriter(line.getOptionValue("out"));
            }

            keepCommonGaps = line.hasOption("keep-common-gaps");
            sorted = line.hasOption("sorted");

            if (line.getArgs().length < 3) {
                throw new Exception("Too few arguments");
            } else if (mode == DerepMode.aligned && line.getArgs().length > 3) {
                throw new Exception("Cannot dereplicate multiple aligned files, try using a mask sequence instead");
            }

            List<Sequence> echoSeqs = new ArrayList();
            long startTime = System.currentTimeMillis();
            long totalSeqs = 0;
            long expectedSeqLength = -1;

            for (int index = 2; index < line.getArgs().length; index++) {
                String arg = line.getArgs()[index];
                String sampleName = new File(arg).getName();
                if (sampleName.contains(".")) {
                    sampleName = sampleName.substring(0, sampleName.lastIndexOf("."));
                }

                System.err.println("Processing " + arg);
                if (mode == DerepMode.model_only) {
                    IndexedSeqReader reader = new IndexedSeqReader(new File(arg), maskId);

                    for (String seqid : reader.getSeqIdSet()) {
                        Sequence seq = reader.readSeq(seqid);

                        if (seq.getSeqName().startsWith("#")) {
                            echoSeqs.add(seq);
                            continue;
                        }
                        totalSeqs++;

                        if(expectedSeqLength == -1) {
                            expectedSeqLength = seq.getSeqString().length();
                        } else if(expectedSeqLength != seq.getSeqString().length()) {
                            throw new IOException("Sequence " + seq.getSeqName() + "'s length (" + seq.getSeqString().length() + ") doesn't match expected length " + expectedSeqLength);
                        }

                        derep.addSeq(seq);
                        sampleMapping.addSeq(sampleName, seq.getSeqName());
                    }
                    reader.close();
                } else {
                    SequenceReader reader = new SequenceReader(new File(arg));
                    Sequence seq;
                    while ((seq = reader.readNextSequence()) != null) {

                        if (seq.getSeqName().startsWith("#")) {
                            echoSeqs.add(seq);
                            continue;
                        }

                        totalSeqs++;

                        if (mode == DerepMode.unaligned) {
                            seq = SeqUtils.getUnalignedSeq(seq);
                        } else if(mode == DerepMode.formatted_model) {
                            seq = SeqUtils.getMaskedBySeqString(seq);
                        }

                        derep.addSeq(seq);

                        sampleMapping.addSeq(sampleName, seq.getSeqName());
                    }
                    reader.close();
                }
            }

            File idMappingFile = new File(line.getArgs()[0]);
            File sampleMappingFile = new File(line.getArgs()[1]);

            sampleMapping.toStream(new PrintStream(new FileOutputStream(sampleMappingFile)));
            
            Set<DerepSeq> uniqueSeqs = null;
            if(!sorted){
                uniqueSeqs = derep.getUniqueSeqs();
            } else {
                uniqueSeqs = derep.getSortedUniqueSeqs();
            }
            int count = 0;
            IdMapping<Integer> idMapping = new IdMapping<Integer>();

            if (mode == DerepMode.unaligned) {
                for (DerepSeq seq : uniqueSeqs) {
                    List<String> ids = seq.ids;
                    idMapping.addIds(count++, ids);
                    seqWriter.writeSeq(seq.seq.getSeqName(), seq.seq.getDesc() + ";size=" + ids.size() + ";", seq.seq.getSeqString());

                    if (qualOut != null && seq.seq instanceof QSequence) {
                        String qseqStr = "";
                        QSequence qseq = (QSequence) seq.seq;
                        for (int index = 0; index < qseq.getQuality().length; index++) {
                            qseqStr += qseq.getQuality()[index] + "  ";
                        }

                        qualOut.writeSeq(seq.seq.getSeqName(), qseqStr);
                    }
                }
            } else {
                Set<Integer> allGapCols = idGapCols(derep);

                for (DerepSeq seq : uniqueSeqs) {
                    List<String> ids = seq.ids;
                    idMapping.addIds(count++, ids);

                    String seqString = seq.seq.getSeqString();

                    if (!keepCommonGaps) {
                        seqString = removeCommonGaps(allGapCols, seqString);
                    }
                    seqWriter.writeSeq(seq.seq.getSeqName(), seq.seq.getDesc() + ";size=" + ids.size() + ";", seqString);
                }

                if (mode == DerepMode.aligned) {
                    for (Sequence seq : echoSeqs) {
                        String seqString = seq.getSeqString();

                        if (!keepCommonGaps) {
                            seqString = removeCommonGaps(allGapCols, seqString);
                        }
                        seqWriter.writeSeq(seq.getSeqName(), seqString);
                    }
                }
            }

            idMapping.toStream(new PrintStream(new FileOutputStream(idMappingFile)));
            System.err.println("Total sequences: " + totalSeqs);
            System.err.println("Unique sequences: " + uniqueSeqs.size());
            System.err.println("Dereplication complete: " + (System.currentTimeMillis() - startTime));

        } catch (Exception e) {
            new HelpFormatter().printHelp("Dereplicator [options] <id-mapping-out> <sample-mapping-out> <seq-file>[,<qual-file>] ...", options);
            e.printStackTrace();
            System.err.println("Error: " + e.getMessage());
            System.exit(1);
        } finally {
            try {
                qualOut.close();
            } catch (Exception ignore) {
            }
            try {
                seqWriter.close();
            } catch (Exception ignore) {
            }
        }
    }
}
