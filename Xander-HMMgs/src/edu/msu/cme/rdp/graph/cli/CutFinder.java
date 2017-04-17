/*
 * Copyright (C) 2012 Michigan State University <rdpstaff at msu.edu>
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

package edu.msu.cme.rdp.graph.cli;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.io.PrintStream;

import edu.msu.cme.rdp.graph.filter.BloomFilter;
import edu.msu.cme.rdp.readseq.SequenceFormat;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import edu.msu.cme.rdp.kmer.Kmer;
import edu.msu.cme.rdp.kmer.NuclKmer;

/**
 * Finds missing kmers in the bloom filter, given reference file(s)
 * containing the known sequence(s).
 * 
 * @author gilmanma
 */
public class CutFinder {
    
    /**
     * A kmer from a reference sequence that will be searched for in the bloom
     * filter
     */
    private static class RefKmer {
        /**
         * The nucleotide string representing the kmer
         */
        private final String str;
        
        /**
         * The location of the kmer relative to the start of the reference
         * sequence
         */
        private final int modelPos;
        
        /**
         * Class constructor
         * 
         * @param str       The nucleotide string representing the kmer
         * @param modelPos  The location of the kmer relative to the start of the reference
         *                  sequence
         */
        public RefKmer(String str, int modelPos) {
            this.str = str;
            this.modelPos = modelPos;
        }
        
        public String getString() {
            return str;
        }
        
        public int getModelPos() {
            return modelPos;
        }
        
        /**
         * Prints the model position and nucleotide string representation of
         * the kmer to the provided output stream
         * 
         * @param out   PrintStream to which to output the information
         *              (usually System.out or System.err)
         */
        public void print(PrintStream out) {
            out.println(String.valueOf(modelPos) + '\t' + str);
        }
    }
    
    private static BloomFilter bloom;
    private static BloomFilter.GraphState bloomState;
    
    /**
     * Loads a file containing the reference sequence(s)
     * 
     * @param fileName      String of the file name
     * @return              a file object of the reference file
     * @throws IOException 
     */
    private static File loadRefFile(String fileName) throws IOException {
        File f = new File(fileName);
        if(!f.exists()) {
            System.err.println("Could not find file: " + fileName);
            return null;
        }

        SequenceFormat format = SeqUtils.guessFileFormat(f);

        if (format == SequenceFormat.UNKNOWN || format == SequenceFormat.EMPTY) {
            return null;
        }
        
        return f;
    }
    
    /**
     * Scan over the sequence string and walk along the bloom filter, printing
     * out all kmers which are not found in the filter.
     * 
     * @param refSeq    reference Sequence to get nucleotides from
     * @param kmerSize  the size of the kmers to extract and check for
     */
    private static void searchKmers(Sequence refSeq, int kmerSize) {
        String seqString = refSeq.getSeqString();
        // create a Kmer; makes it easier to scan the sequence string
        Kmer kmer = new NuclKmer(seqString.substring(0, kmerSize).toCharArray());
        // initialize the bloom filter to the starting kmer
        bloomState.setState(kmer.toString().toCharArray());
        // beginning at the end of the starting kmer, advance along the
        // sequence string
        for(int i = kmerSize; i < seqString.length(); ++i) {
            // create a new reference kmer based on the current kmer
            // makes it easier to print (and store, should that functionality
            // be implemented at a later date)
            RefKmer refKmer = new RefKmer(kmer.toString(), i);
            if (!bloomState.hasCurrent()) {
                refKmer.print(System.out);
            }
            // append the next character to the end of the kmer and advance
            // along the bloom filter
            // **Warning** kmer.shiftLeft and bloomState.shiftLeft go in
            // OPPOSITE directions, which is why bloomState.shiftRight is used
            kmer = kmer.shiftLeft(seqString.charAt(i));
            bloomState.shiftRight(seqString.charAt(i));
        }
    }
    
    /**
     * Main function.  Calling this is how the class is intended to be used.
     * 
     * @param args  Command line arguments
     * @throws IOException 
     */
    public static void main(String args[]) throws IOException {
        if(args.length < 3) {
            System.err.println("USAGE: CutFinder <kmer size> <bloom filter> <query files>");
            System.exit(1);
        }
        
        bloom = BloomFilter.fromFile(new File(args[1]));
        bloomState = bloom.new GraphState();
        
        System.out.println("--< Missing kmers >--");
        
        // extract the reference files from the command line arguments
        String[] seqFileNames = Arrays.copyOfRange(args, 2, args.length);
        for(String fileName : seqFileNames) {
            File f = loadRefFile(fileName);
            if (f == null) {
                continue;
            }
            
            SequenceReader refReader = new SequenceReader(f);
            Sequence seq;
            while ((seq = refReader.readNextSequence()) != null) {
                System.out.println('>' + seq.getSeqName());
                searchKmers(seq, Integer.valueOf(args[0]));
            }
            refReader.close();
        }
    }
}