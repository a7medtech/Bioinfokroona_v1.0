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

import java.util.Arrays;
import java.util.Random;

/**
 *
 * @author fishjord
 */
public class InfinateKmer {

    public static final int a = 0, c = 1, g = 2, t = 3;
    public static final char[] intToChar = new char[]{'a', 'c', 'g', 't'};
    public static final byte[] validateLookup = new byte[127];
    public static final byte[] complementLookup = new byte[127];

    static {
        validateLookup['A'] = a;
        validateLookup['C'] = c;
        validateLookup['G'] = g;
        validateLookup['T'] = t;
        validateLookup['a'] = a;
        validateLookup['c'] = c;
        validateLookup['g'] = g;
        validateLookup['t'] = t;
        validateLookup['U'] = t;
        validateLookup['u'] = t;

        Arrays.fill(complementLookup, (byte) -1);

        complementLookup['A'] = t;  // A
        complementLookup['T'] = a;  // T
        complementLookup['U'] = a;  // U
        complementLookup['G'] = c;  // G
        complementLookup['C'] = g;  // C

        complementLookup[a] = t;  // a
        complementLookup[t] = a;  // t
        complementLookup[g] = c;  // g
        complementLookup[c] = g;  // c
    }
    private final long[] kmer;
    private final int k;
    private final int lastFill;
    private final long lastMask;

    /**
     *
     * acgt
     *    a
     *   ac
     *  acg
     * acgt
     *
     * @param charKmer
     */
    public InfinateKmer(char[] charKmer) {
        this(charKmer.length);

        int i = 0;
        int l = 0;
        byte b;

        /*
         * So this is mildly complicated We want to store the kmer in
         * little-endian byte order so that the shifts make logical sense (right
         * = >> and left = <<), but that means storing the kmer kind of
         * backwards in the array...
         */
        for (int index = 0; index < charKmer.length; index++) {
            b = validateLookup[charKmer[index]];
            if (b == -1) {
                throw new IllegalArgumentException("Kmer contains one or more invalid bases (" + new String(charKmer) + ")");
            }

            kmer[l] = (kmer[l] << 2) | (b & 0x3);

            i++;
            if (i >= 32) {
                i = 0;
                l++;
            }
        }
    }

    public InfinateKmer(long kmer, int k) {
        this(k);
        this.kmer[0] = kmer;
    }

    private InfinateKmer(int k) {
        this.k = k;
        int length = (int) Math.ceil(k / 32.0);
        lastFill = 32 - (32 * length - k);
        long tmp = 0;
        for(int index = 0;index < lastFill;index++) {
            tmp = (tmp << 2) | 3;
        }
        lastMask = tmp;
        kmer = new long[length];
    }

    private InfinateKmer(InfinateKmer k) {
        this.k = k.k;
        this.kmer = new long[k.kmer.length];//Arrays.copyOf(k.kmer, k.kmer.length);
        this.lastFill = k.lastFill;
        this.lastMask = k.lastMask;
    }

    public int length() {
        return k;
    }

    public int packedLength() {
        return kmer.length;
    }

    public long getPart(int index) {
        return kmer[index];
    }

    public InfinateKmer shiftRight(char c) {
        return shiftRight(validateLookup[c]);
    }

    public InfinateKmer shiftLeft(char c) {
        return shiftLeft(validateLookup[c]);
    }

    public InfinateKmer shiftRight(byte b) {
        InfinateKmer ret = new InfinateKmer(this);
        byte overflow = b;
        byte tmp;
        for(int index = 0;index < kmer.length - 1;index++) {
            tmp = (byte)(kmer[index] & 0x3);
            ret.kmer[index] = (kmer[index] >>> 2) | ((long)overflow << 62);
            overflow = tmp;
        }
        ret.kmer[ret.kmer.length - 1] = (kmer[kmer.length - 1] >>> 2) | ((long)overflow << ((ret.lastFill - 1) * 2));
        ret.kmer[ret.kmer.length - 1] &= lastMask;

        return ret;
    }

    public InfinateKmer shiftLeft(byte b) {
        InfinateKmer ret = new InfinateKmer(this);
        byte overflow = b;
        byte tmp;

        tmp = (byte)((kmer[kmer.length - 1] >>> ((lastFill - 1) * 2)) & 0x3);
        ret.kmer[kmer.length - 1] = (kmer[kmer.length - 1] << 2) | ((long)overflow & 0x3);
        ret.kmer[kmer.length - 1] &= lastMask;
        overflow = tmp;

        for(int index = kmer.length - 2;index >= 0;index--) {
            tmp = (byte)((kmer[index] >>> 62) & 0x3);
            ret.kmer[index] = (kmer[index] << 2) | (overflow);
            overflow = tmp;
        }

        return ret;
    }

    @Override
    public String toString() {
        StringBuilder buf = new StringBuilder();
        for (int seg = kmer.length - 1; seg >=0;seg--) {
            int index = 32;
            if (seg == kmer.length - 1) {
                index = lastFill;
            }
            long l = kmer[seg];
            for (; index > 0; index--) {
                buf.append(intToChar[(int) (l & 0x3)]);
                l = l >> 2;
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
        final InfinateKmer other = (InfinateKmer) obj;
        if (!Arrays.equals(this.kmer, other.kmer)) {
            return false;
        }
        return true;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 19 * hash + Arrays.hashCode(this.kmer);
        return hash;
    }

    public static void main(String[] args) {
        char[] t = new char[32];
        Arrays.fill(t, 'g');
        new InfinateKmer(t);
    }
}
