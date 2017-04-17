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

package edu.msu.cme.rdp.kmer;


/**
 *
 * @author fishjord
 */
public abstract class Kmer {
    public static final int max_nucl_kmer_size = 64;
    public static final int max_prot_kmer_size = 24;

    protected long[] kmers = new long[2];
    protected int k;
    protected int l;
    protected int lastFill;
    protected long lastMask;

	/* The following three default are for nucleotides, since it's used most frequently by xander package */
    protected int bitsToshift = 0x2;   // 2 bits for each nucleotides 
    protected int charMask = 0x3 ;  // 11
    protected int itemsPerBucket = 32; // maximum number of nucleotides hold by a long of 64 bits

    public abstract byte charToByte(char c);
    public abstract char intToChar(int i); 
    public abstract Kmer shiftRight(char c);
    public abstract Kmer shiftLeft(char c);
    public abstract Kmer shiftRight(byte b);
    public abstract Kmer shiftLeft(byte b);
     
    public int length() {
        return k;
    } 

    public int packedLength() {
        return l;
    }

    public long getPart(int index) {
       return kmers[index];
    }

    public long[] getLongKmers(){
        return kmers;
    }
    
    public void initialize(char[] charKmer) {
        byte b;

        /*
         * So this is mildly complicated We want to store the kmer in
         * little-endian byte order so that the shifts make logical sense (right
         * = >> and left = <<), but that means storing the kmer kind of
         * backwards in the array...
         */
        for (int index = 0; index < charKmer.length; index++) {
            b = charToByte(charKmer[index]);
            if (b == -1) {
                throw new IllegalArgumentException("Kmer contains one or more invalid bases (" + new String(charKmer) + ")");
            }

            if (index < itemsPerBucket) {
                kmers[0] = (kmers[0] << bitsToshift) | (b & charMask);
            } else {
                kmers[1] = (kmers[1] << bitsToshift) | (b & charMask);
            }
            
        }
    }
   
    protected Kmer shiftRight(byte b, Kmer ret) {

        if (l == 2) {
            byte overflow = (byte) (kmers[0] & charMask);
            ret.kmers[0] = (kmers[0] >>> bitsToshift) | ((long) b << ( (itemsPerBucket-1) *  bitsToshift));
            
            ret.kmers[1] = (kmers[1] >>> bitsToshift) | ((long) overflow << ((ret.lastFill - 1) * bitsToshift));
            ret.kmers[1] &= lastMask;
        } else {
            ret.kmers[0] = (kmers[0] >>> bitsToshift) | ((long) b << ((ret.lastFill - 1) * bitsToshift ));
            ret.kmers[0] &= lastMask;
        }

        return ret;
    }

    protected Kmer shiftLeft(byte b, Kmer ret) {       
        if (l == 2) {
            byte overflow = (byte) ((kmers[1] >>> ((lastFill - 1) * bitsToshift)) & charMask);
            ret.kmers[1] = (kmers[1] << bitsToshift ) | ((long) b & charMask);
            ret.kmers[1] &= lastMask;
            ret.kmers[0] = (kmers[0] << bitsToshift ) | (overflow);
        } else {
            ret.kmers[0] = (kmers[0] << bitsToshift) | ((long) b & charMask);
            ret.kmers[0] &= lastMask;
        }
        return ret;
    }

    @Override
    public String toString() {
        StringBuilder buf = new StringBuilder();
        for (int seg = l; seg > 0; seg--) {
            int index = itemsPerBucket;
            if (seg == l) {
                index = lastFill;
            }
            long thisk = (seg == 2)? kmers[1] : kmers[0];
            for (; index > 0; index--) {
                buf.append( intToChar((int) (thisk & charMask)) );
                thisk = thisk >> bitsToshift;
            }
        }

        return buf.reverse().toString();
    }
    
    
    public String decodeLong(long[] inkmer) {
        StringBuilder buf = new StringBuilder();
        for (int seg = l; seg > 0; seg--) {
            int index = itemsPerBucket;
            if (seg == l) {
                index = lastFill;
            }
            long thisk = (seg == 2)? inkmer[1] : inkmer[0];
            for (; index > 0; index--) {
                buf.append( intToChar((int) (thisk & charMask)) );
                thisk = thisk >> bitsToshift;
            }
        }

        return buf.reverse().toString();
    }
    

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Kmer other = (Kmer) obj;
        if (this.kmers[0] != other.kmers[0] || this.kmers[1] != other.kmers[1]) {
            return false;
        }
        return true;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 59 * hash + (int) (this.kmers[0] ^ (this.kmers[1] >>> itemsPerBucket));
        hash = 59 * hash + (int) (this.kmers[0] ^ (this.kmers[1] >>> itemsPerBucket));
        return hash;
    }
    
}
