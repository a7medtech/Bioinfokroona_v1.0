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
import edu.msu.cme.rdp.kmer.NuclKmer;
import java.util.Iterator;

/**
 *
 * @author Jordan Fish <fishjord at msu.edu>
 */
public class KmerIterator implements Iterator<Kmer> {
    private final char[] bases;
    private Kmer next;
    private int pos;

    public KmerIterator(String seq, int k) {
        this.bases = seq.toCharArray();
        String firstKmer = seq.substring(0, k);
        next = new NuclKmer(firstKmer.toCharArray());
        pos = k;
    }

    public boolean hasNext() {
        return next != null;
    }

    public Kmer next() {
        Kmer ret = next;

        if (pos < bases.length) {
            next = next.shiftLeft(bases[pos]);
            pos++;
        } else {
	    next = null;
	}
        return ret;
    }

    public void remove() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
