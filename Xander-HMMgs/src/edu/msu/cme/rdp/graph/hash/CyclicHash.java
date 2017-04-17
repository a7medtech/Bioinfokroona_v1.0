package edu.msu.cme.rdp.graph.hash;

import java.io.Serializable;

public class CyclicHash implements Hash, Serializable {

    public final static int wordsize = 64;
    private int n;
    private NucleotideHash hasher = NucleotideHash.getInstance();
    static final long serialVersionUID = -8788171152437524877L;

    // myn is the length in characters of the blocks you want to hash
    public CyclicHash(int myn) {
        n = myn;
        if (n > wordsize) {
            throw new IllegalArgumentException();
        }

    }

    private long fastleftshiftn(long x) {
        return (x << n) | (x >>> (wordsize - n));
    }

    private static long fastleftshift1(long x) {
        return (x << 1) | (x >>> (wordsize - 1));
    }

    private static long fastrightshift1(long x) {
        return (x >>> 1 | (x << (wordsize - 1)));
    }

    private long fastrightshiftn(long x) {
        //??
        return (x >>> n | (x << (wordsize - n)));
    }

    // add new character (useful to initiate the hasher)
    // to get a strongly universal hash value, you have to ignore the last or first (n-1) bits.
    public long eatRight(long hashvalue, int c) {
        hashvalue = fastleftshift1(hashvalue);
        hashvalue ^= hasher.hashvalues[c];
        return hashvalue;
    }

    public long eatLeft(long hashvalue, int c) {
        hashvalue = fastrightshift1(hashvalue ^ fastleftshiftn(hasher.hashvalues[c]));
        return hashvalue;
    }

    public long getInitialHashvalue(String s) {
        long hashvalue = 0;
        for (int i = 0; i < n; i++) {
            // add new character (useful to initiate the hasher)
            // to get a strongly universal hash value, you have to ignore the last or first (n-1) bits.
            hashvalue = fastleftshift1(hashvalue);
            int idx = 0;
            switch (s.charAt(i)) {
                case 'A':
                case 'a':
                    idx = 0;
                    break;
                case 'C':
                case 'c':
                    idx = 1;
                    break;
                case 'G':
                case 'g':
                    idx = 2;
                    break;
                case 'T':
                case 't':
                case 'U':
                case 'u':
                    idx = 3;
                    break;
                default:
                    throw new IllegalArgumentException("Can't map character " + s.charAt(i));
            }
            hashvalue ^= hasher.getHashvalue(idx);
        }
        return hashvalue;

    }

    public long replaceRight(long hashvalue, int oldC, int newC) {
        hashvalue = hashvalue ^ hasher.hashvalues[oldC] ^ hasher.hashvalues[newC];
        return hashvalue;
    }

    public long replaceLeft(long hashvalue, int oldC, int newC) {
        hashvalue = hashvalue ^ fastrightshift1(fastleftshiftn(hasher.hashvalues[oldC] ^ hasher.hashvalues[newC]));
        return hashvalue;
    }

    /**
     * removing the right most char, add a char to the left most
     *
     * @param hashvalue
     * @param outchar
     * @param inchar
     * @return
     */
    public long updateLeft(long hashvalue, int outchar, int inchar) {
        long z = hasher.getHashvalue(outchar);
        hashvalue = fastrightshift1(hashvalue ^ z ^ fastleftshiftn(hasher.hashvalues[inchar]));
        return hashvalue;
    }

    /**
     * removing the left most char, add a char to the right most
     *
     * @param hashvalue
     * @param outchar
     * @param inchar
     * @return
     */
    public long updateRight(long hashvalue, int outchar, int inchar) {
        long z = fastleftshiftn(hasher.hashvalues[outchar]);
        hashvalue = fastleftshift1(hashvalue) ^ z ^ hasher.hashvalues[inchar];
        return hashvalue;
    }
}