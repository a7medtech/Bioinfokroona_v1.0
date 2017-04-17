/*
 * Copyright (C) 2012 Jordan Fish <fishjord at msu.edu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package edu.msu.cme.rdp.kmer;

import edu.msu.cme.rdp.kmer.io.KmerStart;
import edu.msu.cme.rdp.kmer.io.KmerStartsWriter;
import edu.msu.cme.rdp.kmer.trie.KmerGenerator;
import edu.msu.cme.rdp.kmer.trie.KmerTrie;
import edu.msu.cme.rdp.kmer.trie.KmerTrie.RefPos;
import edu.msu.cme.rdp.kmer.trie.KmerTrie.TrieLeaf;
import edu.msu.cme.rdp.readseq.SequenceType;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.utils.IUBUtilities;
import edu.msu.cme.rdp.readseq.utils.ProteinUtils;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import java.io.File;
import java.util.Date;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.ReentrantLock;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.cli.Options;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import java.io.IOException;
import java.util.ArrayList;
import org.apache.commons.cli.CommandLine;
import java.util.concurrent.TimeUnit;

/**
 *
 * @author fishjord
 */
public class KmerFilter {

    private static final Options options = new Options();

    static {
        options.addOption("o", "out", true, "Redirect output to file");
        options.addOption("a", "aligned", false, "Build trie from aligned sequences");
        options.addOption("T", "transl-table", true, "Translation table to use when translating nucleotide to protein sequences");
        options.addOption("t", "threads", true, "#Threads to use");
    }
    private static final ReentrantLock outputLock = new ReentrantLock();

    private static void processSeq(Sequence querySeq, List<String> refLabels, KmerTrie kmerTrie, KmerStartsWriter out, int wordSize, boolean translQuery, int translTable, boolean reverse) throws IOException {

        String seqString = querySeq.getSeqString();

        if (reverse) {
            seqString = IUBUtilities.reverseComplement(seqString);
        }

        List<char[]> tmpKmerList = KmerGenerator.getKmers(seqString, wordSize);
        char[][] kmers = tmpKmerList.toArray(new char[tmpKmerList.size()][]);
        char[][][] protKmers = null;

        if (translQuery) {
            protKmers = new char[3][][];
            for (int i = 0; i < 3; i++) {
                String frameSeq = seqString.substring(i);
                //System.out.println(frameSeq);
                frameSeq = ProteinUtils.getInstance().translateToProtein(frameSeq, true, translTable);
                //System.out.println(frameSeq);

                tmpKmerList = KmerGenerator.getKmers(frameSeq, wordSize / 3);
                protKmers[i] = tmpKmerList.toArray(new char[tmpKmerList.size()][]);
            }
        }

        int frame = 0;
        char[] kmer = null;
        char[] protKmer = null;
        TrieLeaf leaf = null;

        for (int kmerIndex = 0; kmerIndex < kmers.length; kmerIndex++) {
            kmer = kmers[kmerIndex];

            if (translQuery) {
                protKmer = null;
                if (protKmers[frame].length > kmerIndex / 3) {
                    protKmer = protKmers[frame][kmerIndex / 3];
                }
                if (protKmer == null) {
                    System.err.println("Skipping null prot kmer " + kmerIndex + " " + querySeq.getSeqName() + ((kmer != null) ? " " + new String(kmer) : ""));
                    continue;
                }
                leaf = kmerTrie.contains(protKmer);
            } else {
                leaf = kmerTrie.contains(kmer);
            }

            if (leaf != null) {
                outputLock.lock();
                try {
                    for (Integer refId : leaf.getRefSets()) {
                        for (RefPos refPos : leaf.getModelStarts(refId)) {

                            out.write(new KmerStart(refLabels.get(refId),
                                    querySeq.getSeqName(),
                                    refPos.seqid,
                                    new String(kmer),
                                    (reverse ? -(frame + 1) : (frame + 1)),
                                    refPos.modelPos,
                                    translQuery,
                                    (translQuery ? new String(protKmer) : null)));
                        }
                    }
                } finally {
                    outputLock.unlock();
                }
            }

            frame = ((frame + 1) >= 3) ? 0 : frame + 1;
        }
    }

    public static void main(String[] args) throws Exception {
        final KmerTrie kmerTrie;
        final SeqReader queryReader;
        final SequenceType querySeqType;
        final File queryFile;
        final KmerStartsWriter out;
        final boolean translQuery;
        final int wordSize;
        final int translTable;
        final boolean alignedSeqs;
        final List<String> refLabels = new ArrayList();
        final int maxThreads;

        try {
            CommandLine cmdLine = new PosixParser().parse(options, args);
            args = cmdLine.getArgs();

            if (args.length < 3) {
                throw new Exception("Unexpected number of arguments");
            }

            if (cmdLine.hasOption("out")) {
                out = new KmerStartsWriter(cmdLine.getOptionValue("out"));
            } else {
                out = new KmerStartsWriter(System.out);
            }

            if (cmdLine.hasOption("aligned")) {
                alignedSeqs = true;
            } else {
                alignedSeqs = false;
            }

            if (cmdLine.hasOption("transl-table")) {
                translTable = Integer.valueOf(cmdLine.getOptionValue("transl-table"));
            } else {
                translTable = 11;
            }

            if (cmdLine.hasOption("threads")) {
                maxThreads = Integer.valueOf(cmdLine.getOptionValue("threads"));
            } else {
                maxThreads = Runtime.getRuntime().availableProcessors();
            }


            queryFile = new File(args[1]);
            wordSize = Integer.valueOf(args[0]);
            SequenceType refSeqType = null;

            querySeqType = SeqUtils.guessSequenceType(queryFile);
            queryReader = new SequenceReader(queryFile);

            if (querySeqType == SequenceType.Protein) {
                throw new Exception("Expected nucl query sequences");
            }


            refSeqType = SeqUtils.guessSequenceType(new File(args[2].contains("=") ? args[2].split("=")[1] : args[2]));

            translQuery = refSeqType == SequenceType.Protein;

            if (translQuery && wordSize % 3 != 0) {
                throw new Exception("Word size must be a multiple of 3 for nucl ref seqs");
            }

            int trieWordSize;
            if (translQuery) {
                trieWordSize = wordSize / 3;
            } else {
                trieWordSize = wordSize;
            }
            kmerTrie = new KmerTrie(trieWordSize, translQuery);

            for (int index = 2; index < args.length; index++) {
                String refName;
                String refFileName = args[index];
                if (refFileName.contains("=")) {
                    String[] lexemes = refFileName.split("=");
                    refName = lexemes[0];
                    refFileName = lexemes[1];
                } else {
                    String tmpName = new File(refFileName).getName();
                    if (tmpName.contains(".")) {
                        refName = tmpName.substring(0, tmpName.lastIndexOf("."));
                    } else {
                        refName = tmpName;
                    }
                }

                File refFile = new File(refFileName);

                if (refSeqType != SeqUtils.guessSequenceType(refFile)) {
                    throw new Exception("Reference file " + refFile + " contains " + SeqUtils.guessFileFormat(refFile) + " sequences but expected " + refSeqType + " sequences");
                }

                SequenceReader seqReader = new SequenceReader(refFile);
                Sequence seq;

                while ((seq = seqReader.readNextSequence()) != null) {
                    if (seq.getSeqName().startsWith("#")) {
                        continue;
                    }
                    if (alignedSeqs) {
                        kmerTrie.addModelSequence(seq, refLabels.size());
                    } else {
                        kmerTrie.addSequence(seq, refLabels.size());
                    }
                }
                seqReader.close();

                refLabels.add(refName);
            }

        } catch (Exception e) {
            new HelpFormatter().printHelp("KmerSearch <word_size> <query_file> [name=]<ref_file> ...", options);
            System.err.println(e.getMessage());
            e.printStackTrace();
            System.exit(1);
            throw new RuntimeException("Stupid jvm");  //While this will never get thrown it is required to make sure javac doesn't get confused about uninitialized variables
        }

        long startTime = System.currentTimeMillis();
        long seqCount = 0;
        final int maxTasks = 25000;

        /*
         * if (args.length == 4) { maxThreads = Integer.valueOf(args[3]); } else {
         */

        //}

        System.err.println("Starting kmer mapping at " + new Date());
        System.err.println("*  Number of threads:       " + maxThreads);
        System.err.println("*  References:              " + refLabels);
        System.err.println("*  Reads file:              " + queryFile);
        System.err.println("*  Kmer length:             " + kmerTrie.getWordSize());

        final AtomicInteger processed = new AtomicInteger();
        final AtomicInteger outstandingTasks = new AtomicInteger();

        ExecutorService service = Executors.newFixedThreadPool(maxThreads);

        Sequence querySeq;

        while ((querySeq = queryReader.readNextSequence()) != null) {
            seqCount++;

            String seqString = querySeq.getSeqString();

            if (seqString.length() < 3) {
                System.err.println("Sequence " + querySeq.getSeqName() + "'s length is less than 3");
                continue;
            }

            final Sequence threadSeq = querySeq;

            Runnable r = new Runnable() {

                public void run() {
                    try {
                        processSeq(threadSeq, refLabels, kmerTrie, out, wordSize, translQuery, translTable, false);
                        processSeq(threadSeq, refLabels, kmerTrie, out, wordSize, translQuery, translTable, true);
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }

                    processed.incrementAndGet();
                    outstandingTasks.decrementAndGet();
                }
            };

            outstandingTasks.incrementAndGet();
            service.submit(r);

            while (outstandingTasks.get() >= maxTasks);

            if ((processed.get() + 1) % 1000000 == 0) {
                System.err.println("Processed " + processed + " sequences in " + (System.currentTimeMillis() - startTime) + " ms");
            }
        }

        service.shutdown();
        service.awaitTermination(1, TimeUnit.DAYS);

        System.err.println("Finished Processed " + processed + " sequences in " + (System.currentTimeMillis() - startTime) + " ms");

        out.close();
    }
}
