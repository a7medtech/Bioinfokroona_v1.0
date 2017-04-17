package edu.msu.cme.rdp.graph.hash;

import java.io.Serializable;
import java.util.Random;

// not to be used directly
public class NucleotideHash implements Serializable {

    public long hashvalues[] = new long[4];
    static final long serialVersionUID = -8788171152437524878L;

    public NucleotideHash() {
        Random r = new Random();
        for (int k = 0; k < hashvalues.length; ++k) {
            hashvalues[k] = r.nextLong();
        }
    }

    public static NucleotideHash getInstance() {
        return charhash;
    }

    public long getHashvalue(int c) {
        if(c > 3 || c < 0) {
            throw new IllegalArgumentException("Expected an integer from 0-3");
        }
        return hashvalues[c];
    }
    static NucleotideHash charhash = new NucleotideHash();
}
