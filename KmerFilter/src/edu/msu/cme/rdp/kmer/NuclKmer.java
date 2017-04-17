
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

import edu.msu.cme.rdp.readseq.utils.NuclBinMapping;
import java.util.Arrays;
import java.util.Random;

/**
 *
 * @author wangqion
 */
public class NuclKmer extends Kmer{
   
    /**
     *
     * @param charKmer
     */
    public NuclKmer(char[] charKmer) {
        this(charKmer.length);
        initialize(charKmer);
    }

    private NuclKmer(int k) {
        if (k > max_nucl_kmer_size) {
            throw new IllegalArgumentException("k-mer lenght must be <= " + max_nucl_kmer_size);
        }

        this.k = k;
        l = (int) Math.ceil(k / ((double)itemsPerBucket));
        lastFill = itemsPerBucket - (itemsPerBucket * l - k);
        long tmp = 0;
        for (int index = 0; index < lastFill; index++) {
            tmp = (tmp << bitsToshift) | charMask;
        }
        lastMask = tmp;

    }
    
    public NuclKmer(NuclKmer k) {        
        this.k = k.k;
        this.l = k.l;
        this.kmers = Arrays.copyOf( k.kmers, k.kmers.length);
        
        this.lastFill = k.lastFill;
        this.lastMask = k.lastMask;
    }
    
    public NuclKmer(long kmer1, int k) {
        this(kmer1, 0, k);
        if(k > itemsPerBucket) {
            throw new IllegalArgumentException("this constructor can only accept k sizes <= " + itemsPerBucket);
        }
    }
    
    public NuclKmer(long kmer1, long kmer2, int k) {
        this(k);
        this.kmers[0] = kmer1;
        this.kmers[1] = kmer2;
    }
        
    
    public Kmer shiftRight(char c) {
        NuclKmer ret = new NuclKmer(this);
        return shiftRight(charToByte(c), ret);
    }

    public Kmer shiftLeft(char c) {
        NuclKmer ret = new NuclKmer(this);
        return shiftLeft(charToByte(c), ret);
    }
    
    public Kmer shiftRight(byte b){
        NuclKmer ret = new NuclKmer(this);
        return shiftRight(b, ret);
    }
    
    public Kmer shiftLeft(byte b){
       NuclKmer ret = new NuclKmer(this);
       return shiftLeft(b, ret); 
    }
    
    public byte charToByte(char c) {
        return NuclBinMapping.validateLookup[c];
    }

    public char intToChar(int i){
        return NuclBinMapping.intToChar[i];
    }

    public static Kmer randomKmer(int l) {
        Random rand = new Random();
        char[] ckmer = new char[l];
        for(int index = 0;index < ckmer.length;index++) {
            ckmer[index] = NuclBinMapping.intToChar[rand.nextInt(4)];
        }

        return new NuclKmer(ckmer);
    }
     
   
}
