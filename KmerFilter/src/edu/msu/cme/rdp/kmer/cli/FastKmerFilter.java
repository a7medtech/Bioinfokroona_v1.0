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
package edu.msu.cme.rdp.kmer.cli;

import edu.msu.cme.rdp.kmer.Kmer;
import edu.msu.cme.rdp.kmer.io.KmerStart;
import edu.msu.cme.rdp.kmer.io.KmerStartsWriter;
import edu.msu.cme.rdp.kmer.set.KmerGenerator;
import edu.msu.cme.rdp.kmer.set.KmerSet;
import edu.msu.cme.rdp.kmer.set.NuclKmerGenerator;
import edu.msu.cme.rdp.kmer.set.ProtKmerGenerator;
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
import java.util.*;
import org.apache.commons.cli.CommandLine;
import java.util.concurrent.TimeUnit;

/**
 *
 * @author fishjord
 */
public class FastKmerFilter {

    private static class RefKmer {

        int modelPos;
        int refFileIndex;
        String refSeqid;

        @Override
        public boolean equals(Object obj) {
            if (obj == null) {
                return false;
            }
            if (getClass() != obj.getClass()) {
                return false;
            }
            final RefKmer other = (RefKmer) obj;
            if (this.modelPos != other.modelPos) {
                return false;
            }
            if (this.refFileIndex != other.refFileIndex) {
                return false;
            }
            return true;
        }

        @Override
        public int hashCode() {
            int hash = 3;
            hash = 37 * hash + this.modelPos;
            hash = 37 * hash + this.refFileIndex;
            return hash;
        }
    }
    private static final Options options = new Options();

    static {
        options.addOption("o", "out", true, "Redirect output to file");
        options.addOption("a", "aligned", false, "Build trie from aligned sequences");
        options.addOption("T", "transl-table", true, "Translation table to use when translating nucleotide to protein sequences");
        options.addOption("t", "threads", true, "#Threads to use");
    }
    private static final ReentrantLock outputLock = new ReentrantLock();

    private static void processSeq(Sequence querySeq, List<String> refLabels, KmerSet<Set<RefKmer>> kmerSet, KmerStartsWriter out, int wordSize, boolean translQuery, int translTable, boolean reverse) throws IOException {

        String seqString = querySeq.getSeqString();

        if (reverse) {
            seqString = IUBUtilities.reverseComplement(seqString);
        }

        KmerGenerator[] kmerGens;

        if (translQuery) {
            kmerGens = new KmerGenerator[3];
            for (int i = 0; i < 3; i++) {
                String frameSeq = seqString.substring(i);
                frameSeq = ProteinUtils.getInstance().translateToProtein(frameSeq, true, translTable);

                kmerGens[i] = new ProtKmerGenerator(frameSeq, wordSize / 3);
            }
        } else {
            kmerGens = new KmerGenerator[1];
            kmerGens[0] = new NuclKmerGenerator(seqString, wordSize);
        }

        int frame = 0;
        Kmer kmer;
        Set<RefKmer> leaves = null;
        List<char[]> nuclKmers = null;
        String nuclKmer;
        String protKmer = null;
        int nuclPos;

        for (int gen = 0; gen < kmerGens.length; gen++) {
            while (kmerGens[gen].hasNext()) {
                kmer = kmerGens[gen].next();
                leaves = kmerSet.get(kmer.getLongKmers());

                if (leaves != null) {
                    outputLock.lock();
                    try {
                        if (nuclKmers == null) {
                            nuclKmers = edu.msu.cme.rdp.kmer.trie.KmerGenerator.getKmers(seqString, wordSize);
                        }

                        if (translQuery) {

                            nuclPos = (kmerGens[gen].getPosition() - 1) * 3 + gen;
                            nuclKmer = new String(nuclKmers.get(nuclPos));
                            protKmer = kmer.toString();
                        } else {
                            nuclPos = kmerGens[gen].getPosition() - 1;
                            nuclKmer = new String(nuclKmers.get(nuclPos));

                        }
                        for (RefKmer refKmer : leaves) {

                            out.write(new KmerStart(refLabels.get(refKmer.refFileIndex),
                                    querySeq.getSeqName(),
                                    refKmer.refSeqid,
                                    nuclKmer,
                                    (reverse ? -(frame + 1) : (frame + 1)),
                                    refKmer.modelPos,
                                    translQuery,
                                    (translQuery? protKmer : null)));
                        }
                    } finally {
                        outputLock.unlock();
                    }
                }
            }
            frame++;
        }
    }

    public static void main(String[] args) throws Exception {
        final KmerSet<Set<RefKmer>> kmerSet;
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
        final int trieWordSize;

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

            if (translQuery) {
                trieWordSize = wordSize / 3;
            } else {
                trieWordSize = wordSize;
            }
            kmerSet = new KmerSet<Set<RefKmer>>();//new KmerTrie(trieWordSize, translQuery);

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

                    KmerGenerator kmers;
		try{
                    if (translQuery) { //protein ref
                        kmers = new ProtKmerGenerator(seq.getSeqString(), trieWordSize, alignedSeqs);
                    } else {
                        kmers = new NuclKmerGenerator(seq.getSeqString(), trieWordSize, alignedSeqs);
                    }
                    while (kmers.hasNext()) {
                        Kmer temp = kmers.next();
                        long[] next = temp.getLongKmers();
                        Set<RefKmer> refKmers = kmerSet.get(next);
                        if (refKmers == null) {
                            refKmers = new HashSet();
                            kmerSet.add(next, refKmers);
                        }

                        RefKmer kmerRef = new RefKmer();
                        kmerRef.modelPos = kmers.getPosition();
                        kmerRef.refFileIndex = refLabels.size();
                        kmerRef.refSeqid = seq.getSeqName();
                        refKmers.add(kmerRef);
                    }
		}catch(IllegalArgumentException ex){
                        //System.err.println(seq.getSeqName()+ " " + ex.getMessage());
                }
		}
                seqReader.close();

                refLabels.add(refName);
            }

        } catch (Exception e) {
            new HelpFormatter().printHelp("KmerSearch <kmerSize> <query_file> [name=]<ref_file> ...\nkmerSize should be multiple of 3, (recommend 45, minimum 30, maximum 63) ", options);
            e.printStackTrace();
            System.exit(1);
            throw new RuntimeException("Stupid jvm");  //While this will never get thrown it is required to make sure javac doesn't get confused about uninitialized variables
        }

        long startTime = System.currentTimeMillis();
        long seqCount = 0;
        final int maxTasks = 25000;

        System.err.println("Starting kmer mapping at " + new Date());
        System.err.println("*  Number of threads:       " + maxThreads);
        System.err.println("*  References:              " + refLabels);
        System.err.println("*  Reads file:              " + queryFile);
        System.err.println("*  Kmer length:             " + trieWordSize);
        System.err.println("*  Kmer Refset Size:        " + kmerSet.size());

        final AtomicInteger processed = new AtomicInteger();
        final AtomicInteger outstandingTasks = new AtomicInteger();

        ExecutorService service = Executors.newFixedThreadPool(maxThreads);

        Sequence querySeq;

        while ((querySeq = queryReader.readNextSequence()) != null) {
            seqCount++;

            String seqString = querySeq.getSeqString();

            if ((!translQuery && seqString.length() < wordSize) || (translQuery && seqString.length() < wordSize + 2)) {
                //System.err.println(querySeq.getSeqName() + "\t" + seqString.length());
                continue;
            }

            final Sequence threadSeq = querySeq;

            Runnable r = new Runnable() {

                public void run() {
                    try {
                        processSeq(threadSeq, refLabels, kmerSet, out, wordSize, translQuery, translTable, false);
                        processSeq(threadSeq, refLabels, kmerSet, out, wordSize, translQuery, translTable, true);

                        processed.incrementAndGet();
                        outstandingTasks.decrementAndGet();
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
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
