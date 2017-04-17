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
package edu.msu.cme.rdp.kmer.set;

import edu.msu.cme.rdp.kmer.Kmer;
import edu.msu.cme.rdp.kmer.ProtKmer;
import edu.msu.cme.rdp.readseq.utils.ProtBinMapping;
import static edu.msu.cme.rdp.readseq.utils.SeqUtils.proteinAlphabet;

/**
 *
 * @author fishjord
 */
public class ProtKmerGenerator implements KmerGenerator {

    private final char[] bases;
    private final int k;
    private Kmer next;

    private int index;     // index in the seqstring
    private int position;  // model position of the kmer found, may not be the kmer returned
    private int curModelPosition; // the model position of the current returning kmer
    private boolean modelOnly = false;
   

    public ProtKmerGenerator(String seq, int k) {
        this(seq, k, false);
    }

    public ProtKmerGenerator(String seq, int k, boolean modelOnly) {
        if (k > 24) {
            throw new IllegalArgumentException("K-mer size cannot be larger than 24");
        }

        if (seq.length() < k) {
            throw new IllegalArgumentException("Sequence length is less than the kmer length");
        }

        this.bases = seq.toCharArray();
        this.k = k;
        this.modelOnly = modelOnly;
        index = 0;
        position = 1;
        next = getFirstKmer(0);
    }


    public boolean hasNext() {
        
        return next != null;
    }

    public Kmer next() {
        Kmer ret = next;
        curModelPosition = position;
        findNextKmer(k -1);
        return ret;
    }

    private Kmer getFirstKmer(int klength){         
        char[] kmerStr = new char[k];
        while (index < bases.length) {
            char base = bases[index++];

            if (modelOnly && (Character.isLowerCase(base) || base == '-' || base == 'X' || base == 'x')) {
                if(base == '-' || base == 'X') {
                    position++;
                }
                klength = 0;
            } else {
                if (!modelOnly || (modelOnly && (base != '.' && proteinAlphabet.contains(base) && base != '*'))) {
                    if (ProtBinMapping.asciiMap[base] == -1) {
                        throw new IllegalArgumentException("Unknown prot base " + base);
                    }
                    kmerStr[klength] = base;
                    position++;
                    klength++;
                }
                
                if (klength == k) {
                    curModelPosition = position;
                    return new ProtKmer(kmerStr);
                }
            }

        }
        
        return null;
    }
    
    private void findNextKmer(int klength){  
        if (next == null){
            return ;
        }

        while (index < bases.length) {
            char base = bases[index++];

            if (modelOnly && (Character.isLowerCase(base) || base == '-' || base == 'X' || base == 'x')) {
                if(base == '-' || base == 'X') {
                    position++;
                }
                klength = 0;
            } else {
                // make sure it's a valid amino acid
                if (!modelOnly ||  (modelOnly && (base != '.' && proteinAlphabet.contains(base) && base != '*')) ) {
                    if (ProtBinMapping.asciiMap[base] == -1) {
                        throw new IllegalArgumentException("Unknown prot base " + base);
                    }
                    next = next.shiftLeft(base);
                    position++;
                    klength++;
                }
                if ( klength == k){
                    return ;
                }
            }

        }
        
        if (klength != k) {  // not a valid kmer
            next = null; 
        }
        
    }
    
    
    /**
     * Get the position of the FIRST CHARACTER in this kmer
     * in the sequence
     * @return
     */
    public int getPosition() {
        return curModelPosition - k ;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Not supported yet.");
    }


    public static void main(String[] args) {
        String seq = "tmamrqcalygkggigkstttqnlvaalaemgkkvmiv";
        ProtKmerGenerator kmer = new ProtKmerGenerator(seq, 20);
        int i = 0;

        while (kmer.hasNext()) {
            Kmer v = kmer.next();
           
            for (int index = 0; index < i; index++) {
                System.out.print(" ");
            }
            System.out.println(v.toString() + " " + kmer.getPosition());
            i++;
        }
        
        seq = "TMAMR-CalGK.GGIGKSTTTQNLVAALAEMGKKVM";
        kmer = new ProtKmerGenerator(seq, 20, true);
        i = 0;
        while (kmer.hasNext()) {
            Kmer v = kmer.next();
           
            for (int index = 0; index < i; index++) {
                System.out.print(" ");
            }
            System.out.println(v.toString() + " " + kmer.getPosition());
            i++;
        }
    }
}
