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

import edu.msu.cme.rdp.kmer.Kmer;
import edu.msu.cme.rdp.kmer.NuclKmer;
import org.junit.*;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class NextCodonTest {
    /**
     * Test of getCodon method, of class NextCodon.
     */
    @Test
    public void testGetCodon() {

        char[] in = {'a', 'c', 'g', 't'};

        Kmer kmer = new NuclKmer(in);
        assertEquals(27, kmer.getPart(0));
        assertEquals("acgt", kmer.toString());

        NextCodon codon = new NextCodon(true, new byte[]{0, 1, 2});
        assertEquals('t', codon.getAminoAcid());
        assertEquals(6, codon.getCodon());
        assertEquals("acg", new NuclKmer(codon.getCodon(), 3).toString());
    }

    @Test
    public void testSanity() {
        int a = 0, c = 1, g = 2, t = 3;
        assertEquals((Character)'s', NextCodon.bacteriaCodonMapping[a][g][t].getAminoAcid());
        assertEquals((Character)'f', NextCodon.bacteriaCodonMapping[t][t][t].getAminoAcid());
        assertEquals((Character)'l', NextCodon.bacteriaCodonMapping[t][t][a].getAminoAcid());
        assertEquals((Character)'g', NextCodon.bacteriaCodonMapping[g][g][t].getAminoAcid());
    }
}
