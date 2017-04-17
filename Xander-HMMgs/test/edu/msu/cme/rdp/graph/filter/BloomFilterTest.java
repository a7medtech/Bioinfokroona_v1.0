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

import edu.msu.cme.rdp.graph.filter.BloomFilter.GraphState;
import edu.msu.cme.rdp.kmer.Kmer;
import edu.msu.cme.rdp.kmer.NuclKmer;
import edu.msu.cme.rdp.readseq.utils.NuclBinMapping;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author wangqion
 */
public class BloomFilterTest {

    int numBits = 1;
    
    public BloomFilterTest() {
    }

    /**
     * Test of addNode method, of class BloomFilter.
     */
    @Test
    public void testAddString() {
        System.err.println("test GraphBuilder AddString()");
        int hashSizeLog2 = 20;
        int hashCount = 3;
        int kmerSize = 10;
        int bitsetSizeLog2 = 16;        

        BloomFilter filter = new BloomFilter(hashSizeLog2, hashCount, kmerSize, bitsetSizeLog2, numBits);
        BloomFilter.GraphBuilder graphBuilder = filter.new GraphBuilder();
        String seq = "aaattgaagagtttgatcatggctcagattgaacgctggcggca";

        graphBuilder.addString(seq.toCharArray());

        //check the first kmer
        graphBuilder.setState(seq.substring(0, kmerSize).toCharArray());
        boolean wasSet = graphBuilder.hasCurrent();
        assertTrue(wasSet);

        // check the last kmer
        graphBuilder.setState(seq.substring(seq.length() - kmerSize, seq.length()).toCharArray());
        wasSet = graphBuilder.hasCurrent();
        assertTrue(wasSet);

        for (int i = seq.length() - kmerSize - 1; i >= 0; i--) {
            graphBuilder.shiftLeft(seq.charAt(i));
            assertTrue(graphBuilder.hasCurrent());
        }
    }

    /**
     * Test of getNextCodon and getSibCodon method, of class BloomFilter.
     */
    @Test
    public void testGetNextRightCodon() {
        System.err.println("test GraphBuilder getNextCodon(), getSibCodon()");
        int hashSizeLog2 = 20;
        int hashCount = 3;
        int kmerSize = 10;
        int bitsetSizeLog2 = 16;

        BloomFilter filter = new BloomFilter(hashSizeLog2, hashCount, kmerSize, bitsetSizeLog2, numBits);
        BloomFilter.GraphBuilder graphBuilder = filter.new GraphBuilder();
        // first kmer: aaattgaaga
        String seq = "aaattgaagaa";
        String seq1 = "aaattgaagagtttgatcatggct";
        String seq2 = "aaattgaagaaatgcat";
        String seq3 = "aaattgaagagttagtat";
        String testMer;

        graphBuilder.addString(seq.toCharArray());
        // test frame 0
        testMer = "aaattgaaga";
        BloomFilter.RightCodonFacade codonFacade = filter.new RightCodonFacade(testMer);
        //System.out.println(codonFacade.getPathString());
        char nextCodon;
        //char nextCodon = codonFacade.getNextCodon().getAminoAcid();  // not aa found
        NextCodon tmp = codonFacade.getNextCodon();
        //System.out.println(codonFacade.getPathString());
        //assertEquals(tmp, null);

        //codonFacade = filter.new RightCodonFacade(testMer, 2);
        //nextCodon = codonFacade.getNextCodon().getAminoAcid();  // found gaa
        //assertEquals(nextCodon, 'e');
        //assertEquals(codonFacade.getPathString(), "a".toUpperCase());

        // test a kmer not in bloomfilter
        try {
            codonFacade = filter.new RightCodonFacade("aaattgaagg");
            fail("should throw IllegalArgumentException because kmer is not in bloomfilter");
        } catch (IllegalArgumentException e) {
        }

        graphBuilder.addString(seq1.toCharArray());
        graphBuilder.addString(seq2.toCharArray());
        graphBuilder.addString(seq3.toCharArray());

        // test frame 0
        codonFacade = filter.new RightCodonFacade(testMer);
        nextCodon = codonFacade.getNextCodon().getAminoAcid();  // first seq2 aaa
        assertEquals(nextCodon, 'k');
        nextCodon = codonFacade.getNextCodon().getAminoAcid();    // tgc
        //System.err.println(codonFacade.getPathString());
        assertEquals(nextCodon, 'c');

        //nextCodon = codonFacade.getNextCodon().getAminoAcid();     // no more, stay on tgc
        assertEquals(codonFacade.getNextCodon(), null);

        assertEquals("aatgc", codonFacade.getPathString());
        //nextCodon = codonFacade.getSibCodon().getAminoAcid();      // back to aaa
        assertEquals(codonFacade.getSibCodon(), null);

        nextCodon = codonFacade.getSibCodon().getAminoAcid();      // find agt
        assertEquals(nextCodon, 's');
        nextCodon = codonFacade.getNextCodon().getAminoAcid();     // find seq3 tag
        assertEquals(nextCodon, '*');
        nextCodon = codonFacade.getNextCodon().getAminoAcid();     // find tat
        assertEquals(nextCodon, 'y');

        //nextCodon = codonFacade.getSibCodon().getAminoAcid();      // back to tag
        assertEquals(codonFacade.getSibCodon(), null);

        nextCodon = codonFacade.getSibCodon().getAminoAcid();      // find seq1 ttg
        assertEquals(nextCodon, 'l');
        nextCodon = codonFacade.getNextCodon().getAminoAcid();     // atc
        assertEquals(nextCodon, 'i');
        assertEquals("gtttgatc", codonFacade.getPathString());

        // test from frame 1
        codonFacade = filter.new RightCodonFacade(testMer);
        nextCodon = codonFacade.getNextCodon().getAminoAcid();  // first seq2 aat
        assertEquals(nextCodon, 'k');
        nextCodon = codonFacade.getNextCodon().getAminoAcid();    // gca
        assertEquals(nextCodon, 'c');

        //nextCodon = codonFacade.getNextCodon().getAminoAcid();     // no more, stay on gca
        assertEquals(codonFacade.getNextCodon(), null);


        try {
            codonFacade = filter.new RightCodonFacade(testMer.substring(1));
            fail("should throw IllegalArgumentException because input length not equal to k-mer length");
        } catch (IllegalArgumentException e) {
        }

        try {
            codonFacade = filter.new RightCodonFacade("naattgaagg");
            fail("should throw InvalidDNABaseException because of invalide bases");
        } catch (InvalidDNABaseException e) {
        }
    }

    /**
     * Test of getNextCodon and getSibCodon method, of class BloomFilter.
     */
    @Test
    public void testLongRight() {
        int hashSizeLog2 = 20;
        int hashCount = 3;
        int kmerSize = 63;
        int bitsetSizeLog2 = 16;

        BloomFilter filter = new BloomFilter(hashSizeLog2, hashCount, kmerSize, bitsetSizeLog2, numBits);
        BloomFilter.GraphBuilder graphBuilder = filter.new GraphBuilder();
        // first kmer: aaattgaaga
        String seq =     "ATGTCTTTGCGCCAGATTGCGTTCTACGGTAAGGGCGGTATCGGAAAGTCCACCACCTCCCAGAACACCCTGGCCGCGCTGGTCGAGCTGGATCAGAAGATCCTGATCGTCGGCTGCGATCCGAAGGCCGACTCGACCCGCCTGATCCTGCACGCCAAGGCGCAGGACACCGTGCTGCACCTCGCCGCCGAAGCCGGCTCGGTCGAGGATCTGGAACTCGAGGACGTTCTCAAGATCGGCTACAAGGGCATCAAGTGCGTCGAGTCCGGCGGTCCGGAGCCGGGGGTCGGCTGCGCCGGCCGCGGCGTGATCACCTCGATCAACTTCCTCGAAGAGAACGGCGCCTACGACGACGTGGACTACGTCTCCTACGACGTGCTGGGCGACGTGGTGTGCGGCGGTTTCGCCATGCCCATCCGCGAGAACAAGGCCCAGGAAATCTACATCGTCATGTCCGGTGAGATGATGGCGCTCTACGCCGCCAACAACATCGCCAAGGGCATTCTGAAGTACGCGCACAGCGGCGGCGTGCGCCTCGGCGGCCTGATCTGCAACGAGCGCCAGACCGACAAGGAAATCGACCTCGCCTCGGCCCTGGCCGCCCGCCTCGGCACCCAGCTCATCCACTTCGTGCCGCGCGACAACATCGTGCAGCACGCCGAGCTGCGCCGCATGACCGTGATCGAGTACGCGCCGGACAGCCAGCAGGCCCAGGAATACCGCCAGCTCGCCAACAAGGTCCACGCGAACAAGGGCAAGGGCACCATCCCGACCCCGATCACGATGGAAGAGCTGGAGGAGATGCTGATGGACTTCGGCATCATGAAGTCGGAGGAGCAGCAGCTCGCCGAGCTCCAGGCCAAGGAAGCCGCCAAGGCCTGA";
        String testMer = "ATGTCTTTGCGCCAGATTGCGTTCTACGGTAAGGGCGGTATCGGAAAGTCCACCACCTCCCAG";

        graphBuilder.addString(seq.toCharArray());
        // test frame 0
        BloomFilter.RightCodonFacade codonFacade = filter.new RightCodonFacade(testMer);
        assertEquals((byte)NuclBinMapping.a, (byte)codonFacade.getNextNucl());
        assertEquals((byte)NuclBinMapping.a, (byte)codonFacade.getNextNucl());
        assertEquals((byte)NuclBinMapping.c, (byte)codonFacade.getNextNucl());
        assertEquals((byte)NuclBinMapping.a, (byte)codonFacade.getNextNucl());
        assertEquals((byte)NuclBinMapping.c, (byte)codonFacade.getNextNucl());
        assertEquals((byte)NuclBinMapping.c, (byte)codonFacade.getNextNucl());
        assertEquals((byte)NuclBinMapping.c, (byte)codonFacade.getNextNucl());
        assertEquals((byte)NuclBinMapping.t, (byte)codonFacade.getNextNucl());
        assertEquals((byte)NuclBinMapping.g, (byte)codonFacade.getNextNucl());
        assertEquals((byte)NuclBinMapping.g, (byte)codonFacade.getNextNucl());
        assertNull(codonFacade.getSibNucl());
        assertNull(codonFacade.getSibNucl());
        assertNull(codonFacade.getSibNucl());
        assertNull(codonFacade.getSibNucl());
        assertNull(codonFacade.getSibNucl());
        assertNull(codonFacade.getSibNucl());
        assertNull(codonFacade.getSibNucl());
        assertNull(codonFacade.getSibNucl());
        assertNull(codonFacade.getSibNucl());
        assertNull(codonFacade.getSibNucl());

    }

    /**
     * Test of getNextCodon and getSibCodon method, of class BloomFilter.
     */
    @Test
    public void testLongRightProt() {
        int hashSizeLog2 = 20;
        int hashCount = 3;
        int kmerSize = 63;
        int bitsetSizeLog2 = 16;

        BloomFilter filter = new BloomFilter(hashSizeLog2, hashCount, kmerSize, bitsetSizeLog2, numBits);
        BloomFilter.GraphBuilder graphBuilder = filter.new GraphBuilder();
        // first kmer: aaattgaaga
        String seq =     "ATGTCTTTGCGCCAGATTGCGTTCTACGGTAAGGGCGGTATCGGAAAGTCCACCACCTCCCAGAACACCCTGGCCGCGCTGGTCGAGCTGGATCAGAAGATCCTGATCGTCGGCTGCGATCCGAAGGCCGACTCGACCCGCCTGATCCTGCACGCCAAGGCGCAGGACACCGTGCTGCACCTCGCCGCCGAAGCCGGCTCGGTCGAGGATCTGGAACTCGAGGACGTTCTCAAGATCGGCTACAAGGGCATCAAGTGCGTCGAGTCCGGCGGTCCGGAGCCGGGGGTCGGCTGCGCCGGCCGCGGCGTGATCACCTCGATCAACTTCCTCGAAGAGAACGGCGCCTACGACGACGTGGACTACGTCTCCTACGACGTGCTGGGCGACGTGGTGTGCGGCGGTTTCGCCATGCCCATCCGCGAGAACAAGGCCCAGGAAATCTACATCGTCATGTCCGGTGAGATGATGGCGCTCTACGCCGCCAACAACATCGCCAAGGGCATTCTGAAGTACGCGCACAGCGGCGGCGTGCGCCTCGGCGGCCTGATCTGCAACGAGCGCCAGACCGACAAGGAAATCGACCTCGCCTCGGCCCTGGCCGCCCGCCTCGGCACCCAGCTCATCCACTTCGTGCCGCGCGACAACATCGTGCAGCACGCCGAGCTGCGCCGCATGACCGTGATCGAGTACGCGCCGGACAGCCAGCAGGCCCAGGAATACCGCCAGCTCGCCAACAAGGTCCACGCGAACAAGGGCAAGGGCACCATCCCGACCCCGATCACGATGGAAGAGCTGGAGGAGATGCTGATGGACTTCGGCATCATGAAGTCGGAGGAGCAGCAGCTCGCCGAGCTCCAGGCCAAGGAAGCCGCCAAGGCCTGA";
        String testMer = "ATGTCTTTGCGCCAGATTGCGTTCTACGGTAAGGGCGGTATCGGAAAGTCCACCACCTCCCAG";

        NextCodon nc;

        graphBuilder.addString(seq.toCharArray());
        // test frame 0
        BloomFilter.RightCodonFacade codonFacade = filter.new RightCodonFacade(testMer);
        nc = codonFacade.getNextCodon();
        assertEquals("Expected n not " + nc.getAminoAcid(), 'n', nc.getAminoAcid());
        nc = codonFacade.getNextCodon();
        assertEquals("Expected t not " + nc.getAminoAcid(), 't', nc.getAminoAcid());
        nc = codonFacade.getNextCodon();
        assertEquals("Expected l not " + nc.getAminoAcid(), 'l', nc.getAminoAcid());
        assertNull(codonFacade.getSibCodon());
        assertNull(codonFacade.getSibCodon());
        assertNull(codonFacade.getSibCodon());

    }

    /**
     * Test of getNextCodon and getSibCodon method, of class BloomFilter.
     */
    @Test
    public void testGetNextLeftCodon() {
        System.err.println("test GraphBuilder getNextCodon(), getSibCodon()");
        int hashSizeLog2 = 20;
        int hashCount = 3;
        int kmerSize = 10;
        int bitsetSizeLog2 = 16;

        BloomFilter filter = new BloomFilter(hashSizeLog2, hashCount, kmerSize, bitsetSizeLog2, numBits);
        BloomFilter.GraphBuilder graphBuilder = filter.new GraphBuilder();
        // first kmer: aaattgaaga
        String seq = "aagaagttaaa";
        String seq1 = "tcggtactagtttgagaagttaaa";
        String seq2 = "tacgtaaagaagttaaa";
        String seq3 = "tatgattgagaagttaaa";
        String testMer;

        graphBuilder.addString(seq.toCharArray());
        // test frame 0
        testMer = "agaagttaaa";
        BloomFilter.LeftCodonFacade codonFacade = filter.new LeftCodonFacade(testMer);

        char nextCodon;
        //char nextCodon = codonFacade.getNextCodon().getAminoAcid();  // not aa found
        assertEquals(null, codonFacade.getNextCodon());

        /*codonFacade = filter.new LeftCodonFacade(testMer, 2);
        nextCodon = codonFacade.getNextCodon().getAminoAcid();  // found gaa
        assertEquals("Expected k, not " + nextCodon, 'k', nextCodon);
        assertEquals(codonFacade.getPathString(), "a".toUpperCase());*/

        // test a kmer not in bloomfilter
        try {
            codonFacade = filter.new LeftCodonFacade("aaattgaagg");
            fail("should throw IllegalArgumentException because kmer is not in bloomfilter");
        } catch (IllegalArgumentException e) {
        }

        graphBuilder.addString(seq1.toCharArray());
        graphBuilder.addString(seq2.toCharArray());
        graphBuilder.addString(seq3.toCharArray());

        // test frame 0
        codonFacade = filter.new LeftCodonFacade(testMer);

        nextCodon = codonFacade.getNextCodon().getAminoAcid();  // first seq2 aaa
        assertEquals("Expected k, not " + nextCodon, 'k', nextCodon);

        nextCodon = codonFacade.getNextCodon().getAminoAcid();    // tgc
        assertEquals("Expected r, not " + nextCodon, 'r', nextCodon);

        //nextCodon = codonFacade.getNextCodon().getAminoAcid();     // no more, stay on tgc
        assertEquals(codonFacade.getNextCodon(), null);

        assertEquals("cgtaa", codonFacade.getPathString());

        //nextCodon = codonFacade.getSibCodon().getAminoAcid();      // back to aaa
        assertEquals(codonFacade.getSibCodon(), null);

        nextCodon = codonFacade.getSibCodon().getAminoAcid();      // find agt
        assertEquals("Expected *, not " + nextCodon, '*', nextCodon);

        nextCodon = codonFacade.getNextCodon().getAminoAcid();     // find seq3 tag
        assertEquals("Expected d, not " + nextCodon, 'd', nextCodon);

        nextCodon = codonFacade.getNextCodon().getAminoAcid();     // find tat
        assertEquals("Expected y, not " + nextCodon, 'y', nextCodon);

        //nextCodon = codonFacade.getSibCodon().getAminoAcid();      // back to tag
        assertEquals(codonFacade.getSibCodon(), null);

        nextCodon = codonFacade.getSibCodon().getAminoAcid();      // find seq1 ttg
        assertEquals("Expected v, not " + nextCodon, 'v', nextCodon);

        nextCodon = codonFacade.getNextCodon().getAminoAcid();     // atc
        assertEquals("Expected l, not " + nextCodon, 'l', nextCodon);
        assertEquals("ctagtttg", codonFacade.getPathString());

        // test from frame 1
        codonFacade = filter.new LeftCodonFacade(testMer);
        nextCodon = codonFacade.getNextCodon().getAminoAcid();  // first seq2 aat
        assertEquals("Expected k, not " + nextCodon, 'k', nextCodon);
        nextCodon = codonFacade.getNextCodon().getAminoAcid();    // gca
        assertEquals("Expected r, not " + nextCodon, 'r', nextCodon);

        //nextCodon = codonFacade.getNextCodon().getAminoAcid();     // no more, stay on gca
        assertEquals(codonFacade.getNextCodon(), null);

        try {
            codonFacade = filter.new LeftCodonFacade(testMer.substring(1));
            fail("should throw IllegalArgumentException because input length not equal to k-mer length");
        } catch (IllegalArgumentException e) {
        }

        try {
            codonFacade = filter.new LeftCodonFacade("naattgaagg");
            fail("should throw InvalidDNABaseException because of invalide bases");
        } catch (InvalidDNABaseException e) {
        }
    }

    /**
     * Test of bloomfilter method, of class BloomFilter.
     */
    @Test
    public void testReload() throws IOException, ClassNotFoundException {
        int hashSizeLog2 = 20;
        int hashCount = 3;
        int kmerSize = 10;
        int bitsetSizeLog2 = 16;
        BloomFilter filter = new BloomFilter(hashSizeLog2, hashCount, kmerSize, bitsetSizeLog2, numBits);
        BloomFilter.GraphBuilder graphBuilder = filter.new GraphBuilder();
        // first kmer: aaattgaaga
        String seq1 = "aaattgaagagtttgatcatggct";
        String seq2 = "aaattgaagaaatgcat";
        String seq3 = "aaattgaagagttagtat";
        graphBuilder.addString(seq1.toCharArray());
        graphBuilder.addString(seq2.toCharArray());
        graphBuilder.addString(seq3.toCharArray());

        File testFile = new File("xx");
        try {
            ObjectOutputStream oos = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(testFile)));

            oos.writeObject(filter);
            oos.close();
            ObjectInputStream ois = new ObjectInputStream(new BufferedInputStream(new FileInputStream(testFile)));
            BloomFilter newFilter = (BloomFilter) ois.readObject();

            // test frame 0
            BloomFilter.RightCodonFacade codonFacade = newFilter.new RightCodonFacade(seq1.substring(0, kmerSize));
            char nextCodon = codonFacade.getNextCodon().getAminoAcid();  // first seq2 aaa
            assertEquals(nextCodon, 'k');
            nextCodon = codonFacade.getNextCodon().getAminoAcid();    // tgc
            assertEquals(nextCodon, 'c');

            //nextCodon = codonFacade.getNextCodon().getAminoAcid();     // no more, stay on tgc
            assertEquals(codonFacade.getNextCodon(), null);

            assertEquals("aatgc", codonFacade.getPathString());

            //nextCodon = codonFacade.getSibCodon().getAminoAcid();      // back to aaa
            assertEquals(codonFacade.getSibCodon(), null);

            nextCodon = codonFacade.getSibCodon().getAminoAcid();      // find agt
            assertEquals(nextCodon, 's');
            nextCodon = codonFacade.getNextCodon().getAminoAcid();     // find seq3 tag
            assertEquals(nextCodon, '*');
            nextCodon = codonFacade.getNextCodon().getAminoAcid();     // find tat
            assertEquals(nextCodon, 'y');

            //nextCodon = codonFacade.getSibCodon().getAminoAcid();      // back to tag
            assertEquals(codonFacade.getSibCodon(), null);

            nextCodon = codonFacade.getSibCodon().getAminoAcid();      // find seq1 ttg
            assertEquals(nextCodon, 'l');
            nextCodon = codonFacade.getNextCodon().getAminoAcid();     // atc
            assertEquals(nextCodon, 'i');
            assertEquals("gtttgatc", codonFacade.getPathString());

            // test from frame 1
            codonFacade = newFilter.new RightCodonFacade(seq1.substring(0, kmerSize));
            nextCodon = codonFacade.getNextCodon().getAminoAcid();  // first seq2 aat
            assertEquals(nextCodon, 'k');
            nextCodon = codonFacade.getNextCodon().getAminoAcid();    // gca
            assertEquals(nextCodon, 'c');

            //nextCodon = codonFacade.getNextCodon().getAminoAcid();     // no more, stay on gca
            assertEquals(codonFacade.getNextCodon(), null);
        } finally {
            if (testFile.exists()) {
                testFile.delete();
            }
        }

    }

    @Test
    public void testRightNucl() {
        int hashSizeLog2 = 20;
        int hashCount = 3;
        int kmerSize = 4;
        int bitsetSizeLog2 = 16;

        BloomFilter filter = new BloomFilter(hashSizeLog2, hashCount, kmerSize, bitsetSizeLog2, numBits);
        BloomFilter.GraphBuilder graphBuilder = filter.new GraphBuilder();
        // first kmer: aaattgaaga
        String seq =  "aaattgacgaa";
        String seq1 = "aaattgtagag";
        //String seq2 = "aaattgaagaaatgcat";
        //String seq3 = "aaattgaagagttagtat";
        String testMer;

        graphBuilder.addString(seq.toCharArray());
        // test frame 0
        testMer = "attg";
        BloomFilter.RightCodonFacade codonFacade = filter.new RightCodonFacade(testMer);
        char nextCodon;

        nextCodon = NuclBinMapping.intToChar[codonFacade.getNextNucl()];
        assertEquals("Expected A not " + nextCodon, 'a', (char)nextCodon);
        nextCodon = NuclBinMapping.intToChar[codonFacade.getNextNucl()];
        assertEquals("Expected C not " + nextCodon, 'c', (char)nextCodon);
        nextCodon = NuclBinMapping.intToChar[codonFacade.getNextNucl()];
        assertEquals("Expected G not " + nextCodon, 'g', (char)nextCodon);
        nextCodon = NuclBinMapping.intToChar[codonFacade.getNextNucl()];
        assertEquals("Expected A not " + nextCodon, 'a', (char)nextCodon);
        nextCodon = NuclBinMapping.intToChar[codonFacade.getNextNucl()];
        assertEquals("Expected A not " + nextCodon, 'a', (char)nextCodon);
        assertNull(codonFacade.getNextNucl());

        graphBuilder.addString(seq1.toCharArray());
        codonFacade = filter.new RightCodonFacade(testMer);

        nextCodon = NuclBinMapping.intToChar[codonFacade.getNextNucl()];
        assertEquals("Expected A not " + nextCodon, 'a', (char)nextCodon);
        nextCodon = NuclBinMapping.intToChar[codonFacade.getSibNucl()];
        assertEquals("Expected T not " + nextCodon, 't', (char)nextCodon);

        nextCodon = NuclBinMapping.intToChar[codonFacade.getNextNucl()];
        assertEquals("Expected A not " + nextCodon, 'a', (char)nextCodon);
        nextCodon = NuclBinMapping.intToChar[codonFacade.getNextNucl()];
        assertEquals("Expected G not " + nextCodon, 'g', (char)nextCodon);
        nextCodon = NuclBinMapping.intToChar[codonFacade.getNextNucl()];
        assertEquals("Expected A not " + nextCodon, 'a', (char)nextCodon);
        nextCodon = NuclBinMapping.intToChar[codonFacade.getNextNucl()];
        assertEquals("Expected G not " + nextCodon, 'g', (char)nextCodon);
        assertNull(codonFacade.getNextNucl());
    }

    @Test
    public void testLeftNucl() {
        int hashSizeLog2 = 20;
        int hashCount = 3;
        int kmerSize = 4;
        int bitsetSizeLog2 = 16;

        BloomFilter filter = new BloomFilter(hashSizeLog2, hashCount, kmerSize, bitsetSizeLog2, numBits);
        BloomFilter.GraphBuilder graphBuilder = filter.new GraphBuilder();
        // first kmer: aaattgaaga
        String seq =  "gatttgacgaa";
        String seq1 = "acattgaagag";
        //String seq2 = "aaattgaagaaatgcat";
        //String seq3 = "aaattgaagagttagtat";
        String testMer;

        graphBuilder.addString(seq.toCharArray());
        // test frame 0
        testMer = "ttga";
        BloomFilter.LeftCodonFacade codonFacade = filter.new LeftCodonFacade(testMer);
        Character nextCodon;

        nextCodon = NuclBinMapping.intToChar[codonFacade.getNextNucl()];
        assertEquals("Expected T not " + nextCodon, 't', (char)nextCodon);
        nextCodon = NuclBinMapping.intToChar[codonFacade.getNextNucl()];
        assertEquals("Expected A not " + nextCodon, 'a', (char)nextCodon);
        nextCodon = NuclBinMapping.intToChar[codonFacade.getNextNucl()];
        assertEquals("Expected G not " + nextCodon, 'g', (char)nextCodon);
        assertNull(codonFacade.getNextNucl());

        graphBuilder.addString(seq1.toCharArray());
        codonFacade = filter.new LeftCodonFacade(testMer);

        nextCodon = NuclBinMapping.intToChar[codonFacade.getNextNucl()];
        assertEquals("Expected A not " + nextCodon, 'a', (char)nextCodon);
        nextCodon = NuclBinMapping.intToChar[codonFacade.getSibNucl()];
        assertEquals("Expected T not " + nextCodon, 't', (char)nextCodon);

        nextCodon = NuclBinMapping.intToChar[codonFacade.getNextNucl()];
        assertEquals("Expected A not " + nextCodon, 'a', (char)nextCodon);
        nextCodon = NuclBinMapping.intToChar[codonFacade.getNextNucl()];
        assertEquals("Expected G not " + nextCodon, 'g', (char)nextCodon);
        assertNull(codonFacade.getNextNucl());
    }

    @Test
    public void testExaustive() {
        int hashSizeLog2 = 20;
        int hashCount = 3;
        int kmerSize = 4;
        int bitsetSizeLog2 = 16;

        BloomFilter filter = new BloomFilter(hashSizeLog2, hashCount, kmerSize, bitsetSizeLog2, numBits);
        BloomFilter.GraphBuilder graphBuilder = filter.new GraphBuilder();
        // first kmer: aaattgaaga
        String seq =  "acgaaa";
        String seq1 = "acgagc";
        String testMer;

        graphBuilder.addString(seq.toCharArray());
        graphBuilder.addString(seq1.toCharArray());
        // test frame 0
        testMer = "acga";
        BloomFilter.RightCodonFacade codonFacade = filter.new RightCodonFacade(testMer);

        char c;

        c = NuclBinMapping.intToChar[codonFacade.getNextNucl()];
        //System.out.println(codonFacade.getPathString());
        assertEquals('a', c);
        c = NuclBinMapping.intToChar[codonFacade.getNextNucl()];
        //System.out.println(codonFacade.getPathString());
        assertEquals('a', c);
        assertNull(codonFacade.getNextNucl());
        //System.out.println(codonFacade.getPathString());
        assertNull(codonFacade.getSibNucl());
        //System.out.println(codonFacade.getPathString());

        c = NuclBinMapping.intToChar[codonFacade.getSibNucl()];
        //System.out.println(codonFacade.getPathString());
        assertEquals('g', c);
        //System.out.println(codonFacade.getPathString());
        c = NuclBinMapping.intToChar[codonFacade.getNextNucl()];
        //System.out.println(codonFacade.getPathString());
        assertEquals("Expected c not " + c, 'c', c);
        assertNull(codonFacade.getNextNucl());
        //System.out.println(codonFacade.getPathString());
        assertNull(codonFacade.getSibNucl());
        //System.out.println(codonFacade.getPathString());
        assertNull(codonFacade.getSibNucl());
        //System.out.println(codonFacade.getPathString());

        assertFalse(codonFacade.hasMoreNucl());
    }

    @Test
    public void testExaustiveProt() {
        int hashSizeLog2 = 20;
        int hashCount = 3;
        int kmerSize = 6;
        int bitsetSizeLog2 = 16;

        BloomFilter filter = new BloomFilter(hashSizeLog2, hashCount, kmerSize, bitsetSizeLog2, numBits);
        BloomFilter.GraphBuilder graphBuilder = filter.new GraphBuilder();
        // first kmer: aaattgaaga
        String seq =  "acgatcagttta";
        String seq1 = "acgatcagtgta";
        String seq2 = "acgatcagtgta";
        String seq3 = "acgatcgctgac";
        String testMer;

        graphBuilder.addString(seq.toCharArray());
        graphBuilder.addString(seq1.toCharArray());
        // test frame 0
        testMer = "acgatc";
        BloomFilter.RightCodonFacade codonFacade = filter.new RightCodonFacade(testMer);


    }
    
    @Test
    public void testMercykmer(){
        int hashSizeLog2 = 30;
        int hashCount = 4;
        int kmerSize = 24;
        int bitsetSizeLog2 = 30;

        BloomFilter filter = new BloomFilter(hashSizeLog2, hashCount, kmerSize, bitsetSizeLog2, 2);
        
        String[] seqs = new String[3];
        seqs[0] = "ttgaaacagattgcattttacggaaaTggagggattggaaagtcaactacagtatg";
        seqs[1] = "ttgaaacagattgcattttacggaaaCggagggattggaaagtcaactacagtatC";
        seqs[2] = "ttgaaacagattgcattttacggaaaGggagggattggaaagtcaacta";
        
        String testMer1 = "aaacagattgcattttacggaaat";
        String testMer2 = "tggagggattggaaagtcaactac";
        String testMer3 = "aacagattgcattttacggaaagg";
        Kmer kmer1 = new NuclKmer(testMer1.toCharArray());
        Kmer kmer2 = new NuclKmer(testMer2.toCharArray());
        Kmer kmer3 = new NuclKmer(testMer3.toCharArray());
        
        BloomFilter.GraphBuilder graphBuilder = filter.new GraphBuilder();
        for ( int i = 0; i < seqs.length; i++){
            graphBuilder.addString(seqs[i].toCharArray());
        }
        // before mercy kmer
        BloomFilter.GraphState bloomState =  filter.new GraphState();
        bloomState.setState(kmer1.toString().toCharArray());
        assertEquals(kmer1.toString(), testMer1);
        assertEquals(bloomState.getMinCurrentCount(), 1);
        bloomState.setState(kmer2.toString().toCharArray());
        assertEquals(kmer2.toString(), testMer2);
        assertEquals(bloomState.getMinCurrentCount(), 1);
        
        // after mercy kmer
        BloomFilter.GraphMercyKmer mercyKmerChecker = filter.new GraphMercyKmer();
        for ( int i = 0; i < seqs.length; i++){
             mercyKmerChecker.checkMercyKmer(seqs[i].toCharArray());
        }
        
        bloomState =  filter.new GraphState();
        bloomState.setState(kmer1.toString().toCharArray());
        assertEquals(kmer1.toString(), testMer1);
        assertEquals(bloomState.getMinCurrentCount(), 2);
        bloomState.setState(kmer2.toString().toCharArray());
        assertEquals(kmer2.toString(), testMer2);
        assertEquals(bloomState.getMinCurrentCount(), 2);
        bloomState.setState(kmer3.toString().toCharArray());
        assertEquals(kmer3.toString(), testMer3);
        assertEquals(bloomState.getMinCurrentCount(), 1); 
        
        
        // 2nd test
        filter = new BloomFilter(hashSizeLog2, hashCount, kmerSize, bitsetSizeLog2, 2);

        graphBuilder = filter.new GraphBuilder();
        for ( int i = 0; i < seqs.length; i++){
            graphBuilder.addString(seqs[i].toCharArray());
        }
        // if we add seqs[2] twice, singleton kmers in seqs[0] and seq[1] should not be promoted to mery kmers.
        graphBuilder.addString(seqs[2].toCharArray());
        
        // after mercy kmer
        mercyKmerChecker = filter.new GraphMercyKmer();
        for ( int i = 0; i < seqs.length; i++){
             mercyKmerChecker.checkMercyKmer(seqs[i].toCharArray());
        }
        
        bloomState =  filter.new GraphState();
        bloomState.setState(kmer1.toString().toCharArray());
        assertEquals(kmer1.toString(), testMer1);
        assertEquals(bloomState.getMinCurrentCount(), 1); 
    }
}
