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
package edu.msu.cme.rdp.kmer.trie;

import edu.msu.cme.rdp.readseq.SequenceType;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

/**
 *
 * @author fishjord
 */
public class ModelPositionKmerGenerator implements Iterator<char[]>, Iterable<char[]> {

    private int index = 0;
    private int modelPosition = 0;
    private int k;
    private int maxIndex;
    private char[] bases;
    private char[] nextKmer = null;
    private final Set<Character> ambigChars;

    public ModelPositionKmerGenerator(String seq, int k, SequenceType seqType) {
        this.k = k;
        this.bases = seq.toCharArray();
        this.maxIndex = bases.length - k + 1;
        
        if(seqType == SequenceType.Protein) {
            ambigChars = SeqUtils.proteinAmbiguity;
        } else if(seqType == SequenceType.Nucleotide) {
            ambigChars = SeqUtils.rrnaAmbiguity;            
        } else {
            throw new IllegalArgumentException("Unknown sequence type " + seqType);
        }
    }

    public static List<char[]> getKmers(String seq, int k, SequenceType seqType) {
        List<char[]> ret = new ArrayList();

        for (char[] kmer : new ModelPositionKmerGenerator(seq, k, seqType)) {
            ret.add(kmer);
        }

        return ret;
    }

    public Iterator<char[]> iterator() {
        return this;
    }

    public boolean hasNext() {
        nextKmer = (nextKmer == null) ? nextKmer() : nextKmer;

        return nextKmer != null;
    }

    private char[] nextKmer() {

        if (index >= maxIndex) {
            return null;
        }

        while (index < bases.length && (bases[index] == '.')) {
            index++;
        }

        char[] kmer = new char[k];
        int currIndex = index;
        int kmerIndex = 0;
        int modelOffset = modelPosition;

        while (currIndex < bases.length && kmerIndex < k) {
            char base = bases[currIndex++];

            if (Character.isUpperCase(base) || base == '-') {
                modelOffset++;
            }
            
            boolean ambig = false;

            if (base == '-' || Character.isLowerCase(base) || ambigChars.contains(base)) {
                kmerIndex = 0;
                index = currIndex;
                modelPosition = modelOffset;
            } else if (base != '.') {
                kmer[kmerIndex++] = base;
            }
        }

        while (index < bases.length && (bases[index] == '.')) {
            index++;
        }

        index++;
        modelPosition++;

        if (kmerIndex == k) {
            return kmer;
        } else {
            return null;
        }
    }

    public int getModelPosition() {
        return modelPosition;
    }

    public char[] next() {
        hasNext();
        char[] ret = nextKmer;
        nextKmer = null;
        return ret;
    }

    public void remove() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
