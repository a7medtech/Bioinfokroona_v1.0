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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

/**
 *
 * @author fishjord
 */
public class KmerGenerator implements Iterator<char[]>, Iterable<char[]> {
    
    private int index = 0;
    private int k;
    private int maxIndex;
    private char[] bases;
    
    public KmerGenerator(String seq, int k) {
        this.k = k;
        this.bases = seq.toLowerCase().toCharArray();
        this.maxIndex = bases.length - k + 1;
    }

    public static List<char[]> getKmers(String seq, int k) {
        List<char[]> ret = new ArrayList();

        for(char[] kmer : new KmerGenerator(seq, k)) {
            ret.add(kmer);
        }

        return ret;
    }

    public Iterator<char[]> iterator() {
        return this;
    }

    public boolean hasNext() {
        return index < maxIndex;
    }

    public char[] next() {
        if(!hasNext()) {
            return null;
        }
        
        int at = index;
        index++;
        
        return Arrays.copyOfRange(bases, at, at + k);
    }

    public void remove() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
