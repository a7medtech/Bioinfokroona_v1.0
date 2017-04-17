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
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class SearchTest {

    int numBits = 1;
    
    @Test
    public void testWalking() {
        int hashSizeLog2 = 20;
        int hashCount = 3;
        int kmerSize = 10;
        int bitsetSizeLog2 = 16;

        BloomFilter filter = new BloomFilter(hashSizeLog2, hashCount, kmerSize, bitsetSizeLog2, numBits);
        BloomFilter.GraphBuilder graphBuilder = filter.new GraphBuilder();
        // first kmer: aaattgaaga
        String seq = "aaattgaagagtttgatcatggct";
        String testMer;

        graphBuilder.addString(seq.toCharArray());
        testMer = "aaattgaaga";
        BloomFilter.RightCodonFacade codonFacade = filter.new RightCodonFacade(testMer);
        long startFwd = codonFacade.getFwdHash();
        long startRc = codonFacade.getRcHash();
        Kmer startKmer = new NuclKmer(testMer.toCharArray());

        Byte[] expected = new Byte[]{codonFacade.getNextNucl(), codonFacade.getNextNucl(), codonFacade.getNextNucl()};

        codonFacade.jumpTo(startKmer, startFwd, startRc);
        assertEquals(expected[0], codonFacade.getNextNucl());
        assertEquals(expected[1], codonFacade.getNextNucl());
        assertEquals(expected[2], codonFacade.getNextNucl());
    }

    @Test
    public void testWalkingNuclShift() {
        int hashSizeLog2 = 20;
        int hashCount = 3;
        int kmerSize = 10;
        int bitsetSizeLog2 = 16;

        BloomFilter filter = new BloomFilter(hashSizeLog2, hashCount, kmerSize, bitsetSizeLog2, numBits);
        BloomFilter.GraphBuilder graphBuilder = filter.new GraphBuilder();
        // first kmer: aaattgaaga
        String seq = "aaattgaagagtttgatcatggct";
        String testMer;

        graphBuilder.addString(seq.toCharArray());
        testMer = "aaattgaaga";
        BloomFilter.RightCodonFacade codonFacade = filter.new RightCodonFacade(testMer);
        Kmer startKmer = new NuclKmer(testMer.toCharArray());

        codonFacade.getNextNucl();
        long startFwd = codonFacade.getFwdHash();
        long startRc = codonFacade.getRcHash();
        Kmer kmer = startKmer.shiftLeft((byte) 2);

        Byte expected = codonFacade.getNextNucl();

        codonFacade.jumpTo(kmer, startFwd, startRc);
        assertEquals(expected, codonFacade.getNextNucl());
    }

    @Test
    public void testWalkingProtShift() {
        int hashSizeLog2 = 20;
        int hashCount = 3;
        int kmerSize = 10;
        int bitsetSizeLog2 = 16;

        BloomFilter filter = new BloomFilter(hashSizeLog2, hashCount, kmerSize, bitsetSizeLog2, numBits);
        BloomFilter.GraphBuilder graphBuilder = filter.new GraphBuilder();
        // first kmer: aaattgaaga
        String seq = "aaattgaagagtttgatcatggct";
        String testMer;

        graphBuilder.addString(seq.toCharArray());
        testMer = "aaattgaaga";
        BloomFilter.RightCodonFacade codonFacade = filter.new RightCodonFacade(testMer);
        Kmer startKmer = new NuclKmer(testMer.toCharArray());

        assertEquals(testMer, startKmer.toString());

        NextCodon nc = codonFacade.getNextCodon();
        long startFwd = codonFacade.getFwdHash();
        long startRc = codonFacade.getRcHash();

        int codon = nc.getCodon();

        byte b1 = (byte) (codon & 0x3);
        byte b2 = (byte) (codon >> 2 & 0x3);
        byte b3 = (byte) (codon >> 4 & 0x3);

        Kmer kmer = startKmer/*.shiftLeft(b3)*/.shiftLeft(b2).shiftLeft(b1);

        //System.err.println(Long.toBinaryString(nc.getCodon()));
        //System.err.println(codonFacade.getPathString());
        assertEquals("agt", new NuclKmer(nc.getCodon(), 3).toString());
        assertEquals("attgaagagt", kmer.toString());

        NextCodon expected = codonFacade.getNextCodon();
        //System.err.println(codonFacade.getPathString());

        codonFacade.jumpTo(kmer, startFwd, startRc);
        //System.err.println(codonFacade.getPathString());
        nc = codonFacade.getNextCodon();
        //System.err.println(codonFacade.getPathString());
        assertEquals(expected, nc);
    }
}
