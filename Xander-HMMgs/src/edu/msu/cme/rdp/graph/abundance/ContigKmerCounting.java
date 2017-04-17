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
package edu.msu.cme.rdp.graph.abundance;

import edu.msu.cme.rdp.kmer.trie.KmerTrie;
import edu.msu.cme.rdp.kmer.trie.KmerTrie.TrieLeaf;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.Date;
import edu.msu.cme.rdp.readseq.utils.IUBUtilities;

/**
 *
 * @author fishjord
 */
public class ContigKmerCounting {

    private KmerTrie trie;

    public ContigKmerCounting(File contigFile, int k) throws IOException {
        SequenceReader reader = new SequenceReader(contigFile);
        trie = KmerTrie.buildTrie(reader, k);
        reader.close();
    }

    public void processSeq(Sequence seq) {
        char[] bases = seq.getSeqString().toCharArray();
        char[] rcbases = IUBUtilities.reverseComplement(seq.getSeqString()).toCharArray();
        for (int index = 0; index <= bases.length - trie.getWordSize(); index++) {
            TrieLeaf leaf = trie.contains(bases, index);

            if (leaf != null) {
                leaf.incQueryCount();
            }

            leaf = trie.contains(rcbases, index);
            if (leaf != null) {
                leaf.incQueryCount();
            }
        }
    }

    public void printResults(File contigFile, PrintStream out) throws IOException {
        SequenceReader reader = new SequenceReader(contigFile);
        Sequence seq;

        while ((seq = reader.readNextSequence()) != null) {
            out.print(seq.getSeqName());

            char[] bases = seq.getSeqString().toCharArray();
            for (int index = 0; index <= bases.length - trie.getWordSize(); index++) {
                TrieLeaf leaf = trie.contains(bases, index);

                out.print("\t");
                if (leaf == null) {
                    System.err.println("ERROR: " + seq.getSeqName() + " contains a kmer not in the trie..." + new String(Arrays.copyOfRange(bases, index, index + trie.getWordSize())));
                    out.print(0);
                } else {
                    out.print(leaf.getQueryCount());
                }
            }
            out.println();
        }

        reader.close();
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 3 && args.length != 4) {
            System.err.println("USAGE: ContigKmerCounting <nucl_contig_file> <reads_file> <k> [#threads]");
            System.exit(1);
        }

        final File nuclContigs = new File(args[0]);
        final File readsFile = new File(args[1]);
        final int k = Integer.valueOf(args[2]);
        final int maxThreads;
        final int maxTasks = 1000;

        if(args.length == 4){
            maxThreads = Integer.valueOf(args[3]);
        } else {
            maxThreads = Runtime.getRuntime().availableProcessors();
        }
        
        System.err.println("Starting kmer mapping at " + new Date());
        System.err.println("*  Number of threads:       " + maxThreads);
        System.err.println("*  Max outstanding tasks:   " + maxTasks);
        System.err.println("*  Nucleotide contigs file: " + nuclContigs);
        System.err.println("*  Reads file:              " + readsFile);
        System.err.println("*  Kmer length:             " + k);

        long startTime = System.currentTimeMillis();
        final ContigKmerCounting kmerCounter = new ContigKmerCounting(nuclContigs, k);
        System.err.println("Kmer trie built in " + (System.currentTimeMillis() - startTime) + " ms");
        
        final AtomicInteger processed = new AtomicInteger();
        final AtomicInteger outstandingTasks = new AtomicInteger();

        
        SequenceReader reader = new SequenceReader(readsFile);
        Sequence seq;
        
        ExecutorService service = Executors.newFixedThreadPool(maxThreads);

        startTime = System.currentTimeMillis();
        while ((seq = reader.readNextSequence()) != null) {
            final Sequence threadSeq = seq;

            Runnable r = new Runnable() {
                public void run() {
                    //System.err.println("Processing sequence " + threadSeq.getSeqName() + " in thread " + Thread.currentThread().getName());
                    kmerCounter.processSeq(threadSeq);
                    //System.err.println("Processed count " + processed);
                    //System.err.println("Outstanding count count " + outstandingTasks);
                    processed.incrementAndGet();
                    outstandingTasks.decrementAndGet();
                }
            };
            
            outstandingTasks.incrementAndGet();
            service.submit(r);
            
            //System.err.println("Submitting " + threadSeq.getSeqName() + ", outstanding tasks= " + outstandingTasks);
            
            while(outstandingTasks.get() >= maxTasks);
            
            if (processed.get() % 1000000 == 0) {
                System.err.println("Processed " + processed + " sequences in " + (System.currentTimeMillis() - startTime) + " ms");
            }
        }

        reader.close();
        
        service.shutdown();

        System.err.println("Processed " + processed + " sequences in " + (System.currentTimeMillis() - startTime) + " ms");

        kmerCounter.printResults(nuclContigs, System.out);
        //System.err.println("Unique kmers in contigs: " + kmerCounter.trie.uniqueWords());
        System.err.println("Processing complete");
    }
}
