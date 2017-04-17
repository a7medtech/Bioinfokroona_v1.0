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
package edu.msu.cme.rdp.graph.filter;

import edu.msu.cme.rdp.kmer.NuclKmer;
import edu.msu.cme.rdp.readseq.utils.ProteinUtils;
import edu.msu.cme.rdp.readseq.utils.ProteinUtils.AminoAcid;
import java.io.Serializable;

/**
 *
 * @author fishjord
 */
public class NextCodon implements Serializable {

    static final AminoAcid[][][] bacteriaCodonMapping = ProteinUtils.getInstance().getTranslationTable(11);

    private int codon;
    private char aminoAcid;

    public NextCodon(boolean forward, byte... codon) {

        if(codon.length != 3) {
            throw new IllegalArgumentException("Codon length should be exactly 3");
        }

        this.codon = ((codon[0] << 4) | (codon[1] << 2) | (codon[2])) & 127;
        try {
            if(forward) {
                this.aminoAcid = bacteriaCodonMapping[codon[0]][codon[1]][codon[2]].getAminoAcid();
            } else {
                this.aminoAcid = bacteriaCodonMapping[codon[2]][codon[1]][codon[0]].getAminoAcid();
            }
        } catch(NullPointerException e) {
            System.out.println(codon[0] + "" + codon[1] + "" + codon[2] + " = " + bacteriaCodonMapping[codon[0]][codon[1]][codon[2]]);
            System.exit(1);
            throw e;
        }
    }

    public int getCodon() {
        return codon;
    }

    public char getAminoAcid() {
        return aminoAcid;
    }

    @Override
    public String toString() {
        return aminoAcid + " " + new NuclKmer(codon, 3) + " " +  Long.toBinaryString(codon);
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final NextCodon other = (NextCodon) obj;
        if (this.codon != other.codon) {
            return false;
        }
        if (this.aminoAcid != other.aminoAcid) {
            return false;
        }
        return true;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 37 * hash + this.codon;
        hash = 37 * hash + this.aminoAcid;
        return hash;
    }
}
