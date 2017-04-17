package edu.msu.cme.rdp.kmer.cli;

import java.io.*;
import java.util.*;
import edu.msu.cme.rdp.readseq.readers.*;
import edu.msu.cme.rdp.kmer.set.*;

public class ExactCounting {
    public static void main(String[] args) throws Exception {
	if(args.length != 2) {
	    System.err.println("USAGE: ExactCounting <seq_file> <k>");
	    System.exit(1);
	}

	File inFile = new File(args[0]);
	File outFile = new File(args[0] + "_kmer_count.txt");
	int k = Integer.valueOf(args[1]);

	KmerSet<Integer> kmers = new KmerSet();
	NuclKmerGenerator kmerGen = null;
	
	SeqReader reader = new SequenceReader(inFile);
	Sequence seq;
	long[] kmer;
	Integer cnt;

	long startTime = System.currentTimeMillis();
	while((seq = reader.readNextSequence()) != null) {
	    kmerGen = new NuclKmerGenerator(seq.getSeqString(), k);
	    while(kmerGen.hasNext()) {
		kmer = kmerGen.next().getLongKmers();

		cnt = kmers.get(kmer);
		if(cnt == null) {
		    cnt = 0;
		}

		kmers.add(kmer, cnt + 1);
	    }
	}
	System.err.println("Kmers loaded in " + (System.currentTimeMillis() - startTime) / 1000.0 + "s");
	kmers.printStats();

	if(kmerGen == null) {
	    System.err.println("No sequences in file");
	    System.exit(1);
	}

	Map<Integer, Integer> occurHist = new HashMap();
	PrintStream out = new PrintStream(outFile);
	for(long[] thisKmer : kmers.getKeys()) {
	    cnt = kmers.get(thisKmer);
	    out.println(thisKmer.toString() + "\t" + cnt);

	    if(!occurHist.containsKey(cnt)) {
		occurHist.put(cnt, 0);
	    }

	    occurHist.put(cnt, occurHist.get(cnt) + 1);
	}

	List<Integer> counts = new ArrayList(occurHist.keySet());
	Collections.sort(counts);

	for(int c : counts) {
	    System.out.println(c + "\t" + occurHist.get(c));
	}
    }
}
