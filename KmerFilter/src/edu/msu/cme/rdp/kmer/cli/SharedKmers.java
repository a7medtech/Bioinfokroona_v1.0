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
import java.io.File;
import java.io.IOException;
import java.util.HashMap;


/**
 *
 * @author wangqion
 */
public class SharedKmers {
    
    private int kmerSize;
    private HashMap<Kmer, Integer> kmerMaps = new HashMap<Kmer, Integer>();  // contains kmer maps for forward and reverse direction
   
    public SharedKmers(int kmer, SequenceReader reader1, SequenceReader reader2) throws IOException{
        this.kmerSize = kmer;
        processRead1File(reader1);
        processRead2File(reader2);
    }
    
    private void processRead1File(SequenceReader reader) throws IOException{
        Sequence seq;
        NuclKmerGenerator kmerGenerator;
        Kmer kmer;
        while ( (seq = reader.readNextSequence()) !=null){
            if ( seq.getSeqString().length() < kmerSize){
                continue;
            }
           
            //only check forward direction
            kmerGenerator = new NuclKmerGenerator(seq.getSeqString(), kmerSize);
            while ( kmerGenerator.hasNext()){
                kmer = kmerGenerator.next();
                if ( !kmerMaps.containsKey(kmer)){
                    kmerMaps.put(kmer, 0);
                }
            }           
        }
        reader.close();
    }
    
    private void processRead2File(SequenceReader reader) throws IOException{
        Sequence seq;
        NuclKmerGenerator kmerGenerator;
        Kmer kmer;
        int notFound = 0;
        while ( (seq = reader.readNextSequence()) !=null){
            if ( seq.getSeqString().length() < kmerSize){
                continue;
            }
           
            //only check forward direction
            kmerGenerator = new NuclKmerGenerator(seq.getSeqString(), kmerSize);
            while ( kmerGenerator.hasNext()){
                kmer = kmerGenerator.next();
                Integer count = kmerMaps.get(kmer);
                if ( count == null){
                    notFound ++;
                } else {
                    kmerMaps.put(kmer, count.intValue() +1);
                }
            }           
        }
        reader.close();
        
        int totalKmers = notFound;
        // get the shared kmer counts
        for ( Integer count: kmerMaps.values()){
            if (count.intValue() == 0){
                notFound ++;
            }
            totalKmers++;
        }
        
        System.out.println("total\t" + totalKmers + "\tshared\t" + (totalKmers -notFound) + "\tpct\t" +  (double)(totalKmers -notFound)/(double)totalKmers);
    }
    
    /**
     * This program takes two files (assuming reads are in correct direction) and compare how many kmers are shared between two files
     * @param args
     * @throws Exception 
     */
    public static void main(String[] args) throws Exception {
        int kmerSize = Integer.parseInt(args[0]);
        if ( kmerSize > Kmer.max_nucl_kmer_size ){
            throw new Exception("kmerSize should be less than " + Kmer.max_nucl_kmer_size);
        }        

        System.out.print(args[1] + "\t" + args[2] + "\t");
        SharedKmers theObj = new SharedKmers(kmerSize, new SequenceReader(new File(args[1])), new SequenceReader(new File(args[2])));
       
    }
}
