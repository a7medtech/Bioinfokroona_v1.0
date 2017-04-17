/*
 * Copyright (C) 2014 wangqion
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
import edu.msu.cme.rdp.kmer.set.NuclKmerGenerator;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.readers.core.SeqReaderCore;
import edu.msu.cme.rdp.readseq.stat.StdevCal;
import edu.msu.cme.rdp.readseq.utils.IUBUtilities;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;
import java.util.Locale;

/**
 *
 * @author wangqion
 */
public class KmerCoverage {
       private static final Options options = new Options();

    static {
	// set Locale to English to use dot as number separator
	Locale.setDefault(Locale.ENGLISH);
        options.addOption("m", "match_reads_out", true, "output the reads containing matching kmers");
        options.addOption("t", "threads", true, "#Threads to use. (default 1)");
    }
 
    private static final String dformat = "%1$.3f";
    private final int kmerSize;
         
    private ConcurrentHashMap<Integer, Contig> contigMap = new ConcurrentHashMap<Integer, Contig>();
    private ConcurrentHashMap<Kmer, KmerAbund>[] kmerMaps = new ConcurrentHashMap[2];   // the number of times kmer occurred in the contigs
    private AtomicInteger totalReads = new AtomicInteger();
    private boolean adjustCount = false;
    
    public class Contig{
        String name;
        // the sum of the weighted count of each kmer (note each match may have different weights since the kmers can be shared by mutliple contigs)      
        double coverage[];        
        
        public Contig(String n, int l){
            name = n;
            coverage = new double[l];
        }
    }
    
    public class ContigCoverage{
        int contigIdx;
        int startPos; // the starting position on the contig, offset is 0
        
        public ContigCoverage(int nameIdx, int p){
            contigIdx = nameIdx;
            startPos = p;
        }
    }
    
    
    public class KmerAbund{
        AtomicInteger count = new AtomicInteger(0);
        ArrayList<ContigCoverage> contigList = new ArrayList<ContigCoverage>();
        
    }
            
            
    /**
     * 
     * @param kmerSize
     * @param contigReader one contig file
     * @throws IOException 
     */
    public KmerCoverage(int kmerSize, SequenceReader contigReader) throws IOException{
        kmerMaps[0] = new ConcurrentHashMap<Kmer, KmerAbund>(); // kmer map for the forward direction
        kmerMaps[1] = new ConcurrentHashMap<Kmer, KmerAbund>(); // kmer map for the reverse direction
        
        this.kmerSize = kmerSize;
        processContigFile(contigReader);
       
    }
    
    /**
     * This is for JUNIT test
     * @param kmerSize
     * @param contigReader
     * @param readsReader
     * @param match_reads_out
     * @throws IOException 
     */
    public KmerCoverage(int kmerSize, SequenceReader contigReader, SeqReaderCore readsReader, PrintStream outStream) throws IOException{
        kmerMaps[0] = new ConcurrentHashMap<Kmer, KmerAbund>(); // kmer map for the forward direction
        kmerMaps[1] = new ConcurrentHashMap<Kmer, KmerAbund>(); // kmer map for the reverse direction
        
        this.kmerSize = kmerSize;

        processContigFile(contigReader);
        Sequence seq;
        while ( (seq = readsReader.readNextSequence()) !=null){
            if ( seq.getSeqString().length() < kmerSize){
                continue;
            }
            processReads(seq, outStream);            
        }
        readsReader.close();
        if ( outStream != null){
            outStream.close();
        }
    }
    
    public ConcurrentHashMap<Integer, Contig> getContigMap(){
        adjustCount();
        return contigMap;
    }
    
    public int getTotalContigs(){
        return contigMap.size();
    }
    /**
     * find the kmers in the contigs
     * @param reader
     * @throws IOException 
     */
    private void processContigFile(SequenceReader reader) throws IOException{
        Sequence seq;
        NuclKmerGenerator kmerGenerator;
        Kmer kmer;
        int contigIdx = 0;
        while ( (seq = reader.readNextSequence()) !=null){
            if ( seq.getSeqString().length() < kmerSize){
                continue;
            }
            // use int to represent seqname in case contig names are too long
            contigMap.put(contigIdx, new Contig(seq.getSeqName(), seq.getSeqString().length() - kmerSize +1));
            //forward direction
            kmerGenerator = new NuclKmerGenerator(seq.getSeqString(), kmerSize);
            while ( kmerGenerator.hasNext()){
                kmer = kmerGenerator.next();
                KmerAbund kmerAbund = kmerMaps[0].get(kmer);

                if ( kmerAbund == null) {
                    kmerAbund = new KmerAbund();
                    kmerMaps[0].put(kmer, kmerAbund);
                }
                kmerAbund.contigList.add(new ContigCoverage(contigIdx, kmerGenerator.getPosition() -1 ));
            }  
            
            // reverse direction
            kmerGenerator = new NuclKmerGenerator(IUBUtilities.reverseComplement(seq.getSeqString()), kmerSize);
            while ( kmerGenerator.hasNext()){
                kmer = kmerGenerator.next();
                KmerAbund kmerAbund = kmerMaps[1].get(kmer);

                if ( kmerAbund == null) {
                    kmerAbund = new KmerAbund();
                    kmerMaps[1].put(kmer, kmerAbund);
                }
                kmerAbund.contigList.add(new ContigCoverage(contigIdx, seq.getSeqString().length() - kmerGenerator.getPosition() -kmerSize +1));
            }  
            contigIdx ++;
        }
        reader.close();
    }
    
    /**
     * This need to be thread safe
     * @param seq
     * @param outStream
     * @throws IOException 
     */
    private void processReads(Sequence seq, PrintStream outStream) throws IOException{
        NuclKmerGenerator kmerGenerator;
        Kmer kmer;
        
        boolean found = false;            
        kmerGenerator = new NuclKmerGenerator(seq.getSeqString(), kmerSize);
        while ( kmerGenerator.hasNext()){
            
            kmer = kmerGenerator.next();  
            for ( int i = 0; i < kmerMaps.length; i ++){  // for forward and reverse direction
                KmerAbund kmerAbund = kmerMaps[i].get(kmer);
                if ( kmerAbund != null) {                   
                    // increment the count
                    kmerAbund.count.addAndGet(1); 
                    found = true;
                }     
                
            }              
        }   
        if ( found ){
            totalReads.incrementAndGet();
        }
        if ( outStream != null && found){
            writeSeq(seq, outStream);
        }        
    }
    
    private synchronized void writeSeq(Sequence seq, PrintStream outStream){
         outStream.println(">" + seq.getSeqName() + "\n" + seq.getSeqString());
    }
    
    /**
     * This should be only called once.
     */
    private synchronized void adjustCount(){
       if ( adjustCount ) return;
         // need to adjust the count
       for ( int i = 0; i < kmerMaps.length; i ++){  // for forward and reverse direction
            for ( Kmer kmer: kmerMaps[i].keySet()){
                KmerAbund kmerAbund = kmerMaps[i].get(kmer);
                // we assign an eqaul value to all the contigs containing the kmer.
                double weightedCount = (double)kmerAbund.count.get()/ (double)kmerAbund.contigList.size();
                for ( ContigCoverage contigCov: kmerAbund.contigList){
                     Contig contig = contigMap.get(contigCov.contigIdx);
                     contig.coverage[contigCov.startPos] += weightedCount ;
                } 
                
            }              
       }  
       adjustCount = true;
    }
  
    
   public void printCovereage(OutputStream coverage_out, OutputStream abundance_out) throws IOException{
       adjustCount();
        // print out the weighted kmer coverage
       // we found mean coverage matched the previous biological observation
        PrintStream coverage_outStream = new PrintStream(coverage_out);
        coverage_outStream.println("#total reads: " + totalReads.intValue());
        coverage_outStream.println("#use mean_cov to adjust the contig abundance, not median_cov ");
        coverage_outStream.println("#seqid\tmean_cov\tmedian_cov\ttotal_pos\tcovered_pos\tcovered_ratio");
        
        for ( Contig contig: contigMap.values()){
            ArrayList<Double> counts = new ArrayList<Double>();
            int coveredPos = 0;
            for ( int pos = 0; pos < contig.coverage.length; pos++){
                if ( contig.coverage[pos] > 0){
                    coveredPos++;
                }
                counts.add(contig.coverage[pos]);
            }
            if ( coveredPos > 0){
                coverage_outStream.println(contig.name + "\t" + String.format(dformat, StdevCal.calMean(counts))
                        +"\t" + String.format(dformat, (StdevCal.calMedian(counts)))      
                + "\t" + counts.size() + "\t" + coveredPos
                + "\t" + String.format(dformat,(double)coveredPos / (double)contig.coverage.length));
            }else { // no coverage
                coverage_outStream.println(contig.name +"\t" + 0 + "\t" + 0 + "\t" + contig.coverage.length + "\t" + 0 + "\t" + 0);
            }
        }
        coverage_outStream.close();
        
       // print kmer abundance
        HashMap<Integer, Integer> abundanceCountMap = new HashMap<Integer, Integer>(); // the frequeny of the kmer abundance         
        PrintStream abundance_outStream = new PrintStream(abundance_out);        
        // need to merge the counts from forward and reverse together.
        HashSet<Kmer> kmerSet = new HashSet<Kmer>();
        kmerSet.addAll(kmerMaps[0].keySet());
        for ( Kmer kmer: kmerSet){
            AtomicInteger abundance = kmerMaps[0].get(kmer).count;
                  
            String reverseKmerStr = IUBUtilities.reverseComplement(kmer.decodeLong(kmer.getLongKmers()));
            Kmer reverseKmer = (new NuclKmerGenerator(reverseKmerStr, this.kmerSize)).next();            
            KmerAbund kmerAbund = kmerMaps[1].get(reverseKmer);

            if ( kmerAbund != null){
                abundance.addAndGet(kmerAbund.count.get());
            }
             
            Integer count = abundanceCountMap.get(abundance.get());
            if ( count == null){
                abundanceCountMap.put(abundance.get(), 1);
            }else {
                abundanceCountMap.put(abundance.get(), count + 1);
            }
        }        
        
        abundance_outStream.println("kmer_abundance\tfrequency");
        for ( Integer abundance: abundanceCountMap.keySet()){
            abundance_outStream.println(abundance + "\t" + abundanceCountMap.get(abundance));
        }
        abundance_outStream.close();
    }
    
       
    /**
     * This program maps the kmers from reads to kmers on each contig,
     * writes the mean, median coverage of each contig to a file
     * writes the kmer abundance to a file
     * @param args
     * @throws IOException 
     */
    public static void main(String[] args) throws IOException, InterruptedException {
        int kmerSize = 45;
        final int maxThreads;
        final int maxTasks = 1000;
        final PrintStream match_reads_out ;
        try {
            CommandLine cmdLine = new PosixParser().parse(options, args);
            args = cmdLine.getArgs();
            if (args.length < 5) {
                throw new Exception("Unexpected number of arguments");
            }
            kmerSize = Integer.parseInt(args[0]);
            if ( kmerSize > Kmer.max_nucl_kmer_size ){
                throw new Exception("kmerSize should be less than " + Kmer.max_nucl_kmer_size);
            }        
            if (cmdLine.hasOption("match_reads_out")) {
                 match_reads_out = new PrintStream(cmdLine.getOptionValue("match_reads_out"));
            } else {
                match_reads_out = null;
            }
             if (cmdLine.hasOption("threads")) {
                maxThreads = Integer.valueOf(cmdLine.getOptionValue("threads"));
                if ( maxThreads >= Runtime.getRuntime().availableProcessors()) {
                   System.err.println(" Runtime.getRuntime().availableProcessors() " + Runtime.getRuntime().availableProcessors()); 
                }
                
            } else {
                maxThreads = 1;
            }
                   
            final KmerCoverage kmerCoverage = new KmerCoverage( kmerSize, new SequenceReader(new File(args[1])));
            if ( kmerCoverage.getTotalContigs() == 0){
                System.out.println("Found 0 contig with length >= kmer size " + kmerSize + " in input file " + args[1] + ". Exit program.");
                return;
            }
            final AtomicInteger outstandingTasks = new AtomicInteger();        
            ExecutorService service = Executors.newFixedThreadPool(maxThreads);

            Sequence seq;
            
            // parse one file at a time
            for (int index = 4; index < args.length; index++) {

                SequenceReader reader = new SequenceReader(new File(args[index]));
                while ( (seq = reader.readNextSequence()) !=null){
                    if ( seq.getSeqString().length() < kmerSize){
                        continue;
                    }
                    final Sequence threadSeq = seq;
                   
                    Runnable r = new Runnable() {

                        public void run() {
                            try {
                                kmerCoverage.processReads(threadSeq, match_reads_out);
                                outstandingTasks.decrementAndGet();
                            } catch (Exception e) {
                                e.printStackTrace();
                            }
                        }
                    };

                    outstandingTasks.incrementAndGet();
                    service.submit(r);
                    
                    while (outstandingTasks.get() >= maxTasks);                

                }
                reader.close();
            }
            service.shutdown();
            service.awaitTermination(1, TimeUnit.DAYS);        

            
            kmerCoverage.printCovereage(new FileOutputStream(new File(args[2])), new FileOutputStream(new File(args[3])));
            if ( match_reads_out != null){
                match_reads_out.close();
            }
        }catch (Exception e) {            
            new HelpFormatter().printHelp("KmerCoverage <kmerSize> <query_file> <coverage_out> <abundance_out> <reads_file> <reads_file>...\nmaximum kmerSize " + Kmer.max_nucl_kmer_size , options);
            e.printStackTrace();
            System.exit(1);
        }
    }
}
