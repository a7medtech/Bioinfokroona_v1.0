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
import edu.msu.cme.rdp.readseq.utils.NuclBinMapping;
import edu.msu.cme.rdp.readseq.utils.ProtBinMapping;

/**
 *
 * @author fishjord
 */
public class NuclKmerGenerator implements KmerGenerator {

    private final char[] bases;
    private final int k;
    private Kmer next;
    private int index;     // index in the seqstring
    private int position;  // model position of the kmer found, may not be the kmer returned
    private int curModelPosition; // the model position of the current returning kmer
    private boolean modelOnly = false;

    public NuclKmerGenerator(String seq, int k) {
        this(seq, k, false);
    }


    public NuclKmerGenerator(String seq, int k, boolean modelOnly) {
        if (k > Kmer.max_nucl_kmer_size) {
            throw new IllegalArgumentException("K-mer size cannot be larger than " + Kmer.max_nucl_kmer_size);
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
        findNextKmer(k - 1);
        return ret;
    }

    private Kmer getFirstKmer(int klength) {
        char[] kmerStr = new char[k];
        while (index < bases.length) {
            char base = bases[index++];


            if (modelOnly && (Character.isLowerCase(base) || base == '-' )) {
                if (base == '-') {
                    position++;
                }
                klength = 0;
            } else {
                if (base == 'N' || base == 'n') {
                    klength = 0;
                    position++;
                } else if (base != '.') {
                    if (NuclBinMapping.validateLookup[base] == -1) {
                        throw new IllegalArgumentException("Only nucleotide bases excepted, not '" + base + "'");
                    }

                    kmerStr[klength] = base;
                    position++;
                    klength++;
                }

                if (klength == k) {
                    curModelPosition = position;
                    return new NuclKmer(kmerStr);
                }
            }
        }

        return null;
    }

    private void findNextKmer(int klength) {
        if (next == null) {
            return;
        }

        while (index < bases.length) {
            char base = bases[index++];

            if (modelOnly && (Character.isLowerCase(base) || base == '-')) {
                if (base == '-') {
                    position++;
                }
                klength = 0;
            } else {
                if (base == 'N' || base == 'n') {
                    klength = 0;
                    position++;
                } else if (base != '.') {
                    if (NuclBinMapping.validateLookup[base] == -1) {
                        throw new IllegalArgumentException("Only nucleotide bases excepted, not '" + base + "'");
                    }

                    next = next.shiftLeft(base);
                    position++;
                    klength++;

                    if (klength == k) {
                        return;
                    }
                }
            }
        }
        if (klength != k) {  // not a valid kmer
            next = null;
        }

    }


    /**
     * Get the position of the FIRST CHARACTER in this kmer in the sequence
     *
     * @return
     */
    public int getPosition() {
        return curModelPosition - k;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    public static void main(String[] args) {
        String seq = "GTTCACGCGCTCGAACTGGTTCCCGGCAAGGGCGCTCAGCTCGCCCGCAGCGCCGGCACCCGCGTACNCAGGTTCAGGCTCGAACTGGTTCCCGGCAAGGGCGCTCAGCTCGCCCGCAGCGCCGGCACCCGCG";
        NuclKmerGenerator kmer = new NuclKmerGenerator(seq, 60);
        int i = 0;

        while (kmer.hasNext()) {
            Kmer v = kmer.next();
           
            for (int index = 0; index < i; index++) {
                System.out.print(" ");
            }
            System.out.println(v.toString() + " " + kmer.getPosition());
            i++;
        }
        
        seq = "GTTCAC-CG.CTCGAACTGGTTCCCGGCAAGGGCGCTCAGCTCGCCCGCAGCGCCGGCACCCGCGTACNCAGGTTCAG";
        kmer = new NuclKmerGenerator(seq, 60, true);
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
