
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

import edu.msu.cme.rdp.readseq.utils.ProtBinMapping;
import java.util.Arrays;

/**
 *
 * @author wangqion
 */
public class ProtKmer extends Kmer {
       
    /**
     *
     * @param charKmer
     */
    public ProtKmer(char[] charKmer) {
        this(charKmer.length);
        initialize(charKmer);
    }


    private ProtKmer(int k) {
        if (k > max_prot_kmer_size) {  // 2 * itemsPerBucket
            throw new IllegalArgumentException("k-mer lenght must be <= " + max_prot_kmer_size);
        }

        setUp(); 
        this.k = k;
        l = (int) Math.ceil(k / ((double)itemsPerBucket));
        lastFill = itemsPerBucket - (itemsPerBucket * l - k);
        long tmp = 0;
        for (int index = 0; index < lastFill; index++) {
            tmp = (tmp << bitsToshift) | charMask;
        }
        lastMask = tmp;
    }

    public ProtKmer(ProtKmer k) {        
        this.k = k.k;
        this.l = k.l;
        this.kmers = Arrays.copyOf( k.kmers, k.kmers.length);

        this.lastFill = k.lastFill;
        this.lastMask = k.lastMask;
        setUp();
    }
    
    private void setUp(){
        bitsToshift = 0x5;   // 5 bits for each amino acid
        charMask = 0x1F ;  // 11111
        itemsPerBucket = 12;  // maximum number of amino acids holded by a long of 64 bits 
    }
        
    public Kmer shiftRight(char c) {
        ProtKmer ret = new ProtKmer(this);
        return shiftRight(charToByte(c), ret);
    }

    public Kmer shiftLeft(char c) {
        ProtKmer ret = new ProtKmer(this);
        return shiftLeft(charToByte(c), ret);
    }
    
    public Kmer shiftRight(byte b){
        ProtKmer ret = new ProtKmer(this);
        return shiftRight(b, ret);
    }
    
    public Kmer shiftLeft(byte b){
       ProtKmer ret = new ProtKmer(this);
       return shiftLeft(b, ret); 
    }
    
    public byte charToByte(char c) {
        return ProtBinMapping.asciiMap[c];
    }

    public char intToChar(int i){
        return ProtBinMapping.intToChar[i];
    }
}