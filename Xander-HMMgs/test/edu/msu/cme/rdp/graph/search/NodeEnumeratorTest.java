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
package edu.msu.cme.rdp.graph.search;

import edu.msu.cme.rdp.alignment.hmm.ProfileHMM;
import edu.msu.cme.rdp.alignment.hmm.TSC;
import edu.msu.cme.rdp.alignment.hmm.XSC;
import edu.msu.cme.rdp.alignment.hmm.XSTATES;
import edu.msu.cme.rdp.graph.filter.BloomFilter;
import edu.msu.cme.rdp.graph.filter.BloomFilter.GraphBuilder;
import edu.msu.cme.rdp.graph.filter.BloomFilter.RightCodonFacade;
import edu.msu.cme.rdp.graph.filter.NextCodon;
import edu.msu.cme.rdp.graph.search.HMMGraphSearch.PartialResult;
import edu.msu.cme.rdp.graph.search.heuristic.weight.HeuristicWeight;
import edu.msu.cme.rdp.graph.search.heuristic.weight.StaticHeuristicWeight;
import edu.msu.cme.rdp.kmer.Kmer;
import edu.msu.cme.rdp.kmer.NuclKmer;
import edu.msu.cme.rdp.readseq.SequenceType;
import edu.msu.cme.rdp.readseq.utils.NuclBinMapping;
import java.util.*;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class NodeEnumeratorTest {

    int numBits = 1;
    
    private static class MockProfileHMM extends ProfileHMM {

        private final int k, m;
        private final SequenceType t;

        public MockProfileHMM(int k, int m, SequenceType t) {
            this.k = k;
            this.m = m;
            this.t = t;
        }

        @Override
        public int K() {
            return k;
        }

        @Override
        public int M() {
            return m;
        }

        @Override
        public SequenceType getAlphabet() {
            return t;
        }

        @Override
        public double getMaxMatchEmission(int i) {
            return -.25;
        }

        @Override
        public double isc(int k, char b) {
            return -1;
        }

        @Override
        public double isc(int k, int b) {
            return -1;
        }

        @Override
        public double msc(int k, char b) {
            return 0;
        }

        @Override
        public double msc(int k, int b) {
            return 0;
        }

        @Override
        public void reconfigureLength(int L) {
        }

        @Override
        public void rescaleMatchEmission(int k, char b, double scale) {
        }

        @Override
        public double tsc(int k, TSC trans) {
            return -.5;
        }

        @Override
        public double[] tsc(TSC trans) {
            return null;
        }

        @Override
        public double xsc(XSTATES xstate, XSC trans) {
            return -1;
        }

        @Override
        public void xsc(XSTATES xstate, XSC trans, double val) {
        }
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
        AStarNode start = new AStarNode(null, new NuclKmer(testMer.toCharArray()), codonFacade.getFwdHash(), codonFacade.getRcHash(), 0, 'm');
        AStarNode curr;
        Set<AStarNode> neighbors;
        // this needs test
        HeuristicWeight hweight = new StaticHeuristicWeight(0.0);
        NodeEnumerator ne = new NodeEnumerator(new MockProfileHMM(20, 10, SequenceType.Protein), hweight);

        curr = start;
        for ( int i = 0; i <= 7; i++){
            neighbors = ne.enumerateNodes(curr, codonFacade);
            assertEquals((curr.state == 'm')? 3 : 2, neighbors.size());
            curr = neighbors.iterator().next();
            nc = kmerToCodon(curr.kmer);
            //System.out.println("curr char " + (char) nc.getAminoAcid() + " " + curr.state + " " + curr.stateNo + " from " + curr.discoveredFrom.state + " " + curr.discoveredFrom.stateNo  + " " + curr.kmer.toString());

            assertEquals("Expected q not " + nc.getAminoAcid(), 'q', nc.getAminoAcid());
        }
        
        neighbors = ne.enumerateNodes(curr, codonFacade);
        assertEquals((curr.state == 'm')? 3 : 2, neighbors.size());
        curr = neighbors.iterator().next();
        nc = kmerToCodon(curr.kmer);
        assertEquals("Expected n not " + nc.getAminoAcid(), 'n', nc.getAminoAcid());

        neighbors = ne.enumerateNodes(curr, codonFacade);
        assertEquals((curr.state == 'm')? 3 : 2, neighbors.size());
        curr = neighbors.iterator().next();
        nc = kmerToCodon(curr.kmer);
        assertEquals("Expected t not " + nc.getAminoAcid(), 't', nc.getAminoAcid());

        neighbors = ne.enumerateNodes(curr, codonFacade);
        assertEquals((curr.state == 'm')? 3 : 2, neighbors.size());
        curr = neighbors.iterator().next();
        nc = kmerToCodon(curr.kmer);
        assertEquals("Expected l not " + nc.getAminoAcid(), 'l', nc.getAminoAcid());
               
    }

    private static NextCodon kmerToCodon(Kmer kmer) {
        String s = kmer.toString();
        char[] charmer = s.substring(s.length() - 3).toCharArray();

        return new NextCodon(true, NuclBinMapping.validateLookup[charmer[0]], NuclBinMapping.validateLookup[charmer[1]], NuclBinMapping.validateLookup[charmer[2]]);
    }

    @Test
    public void testNuclEnumerate() {
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

        AStarNode start = new AStarNode(null, new NuclKmer(testMer.toCharArray()), codonFacade.getFwdHash(), codonFacade.getRcHash(), 0, 'm');

        // this needs test
        HeuristicWeight hweight = new StaticHeuristicWeight(0.0);
        NodeEnumerator ne = new NodeEnumerator(new MockProfileHMM(4, 10, SequenceType.Nucleotide), hweight);
        List<PriorityQueue<AStarNode>> neighbors = new ArrayList();
        neighbors.add(new PriorityQueue());
        neighbors.get(0).add(start);

        for (int index = 0; index < 5; index++) {
            PriorityQueue<AStarNode> rank = neighbors.get(neighbors.size() - 1);
            PriorityQueue<AStarNode> next = new PriorityQueue();
            for (AStarNode node : rank) {
                next.addAll(ne.enumerateNodes(node, codonFacade));
            }

            neighbors.add(next);
        }

        PriorityQueue<AStarNode> rank = neighbors.get(neighbors.size() - 1);

        int i = 0;
        for (AStarNode node : rank) {
            PartialResult result = HMMGraphSearch.partialResultFromGoal(node, true, false, 10, 1);
            assertTrue(nuclExpected.contains(result.alignment));
            nuclExpected.remove(result.alignment);
        }

        assertTrue(nuclExpected.isEmpty());
    }

    @Test
    public void testProtEnumerate() {
        int hashSizeLog2 = 20;
        int hashCount = 3;
        int kmerSize = 6;
        int bitsetSizeLog2 = 16;

        BloomFilter filter = new BloomFilter(hashSizeLog2, hashCount, kmerSize, bitsetSizeLog2, numBits);
        BloomFilter.GraphBuilder graphBuilder = filter.new GraphBuilder();
        // first kmer: aaattgaaga
        String seq = "aaattgaagagtttgatcatggct";
        String testMer;

        graphBuilder.addString(seq.toCharArray());
        testMer = "aaattg";
        BloomFilter.RightCodonFacade codonFacade = filter.new RightCodonFacade(testMer);

        AStarNode start = new AStarNode(null, new NuclKmer(testMer.toCharArray()), codonFacade.getFwdHash(), codonFacade.getRcHash(), 0, 'm');

        // this needs test
        HeuristicWeight hweight = new StaticHeuristicWeight(0.0);
        NodeEnumerator ne = new NodeEnumerator(new MockProfileHMM(20, 10, SequenceType.Protein), hweight);
        List<PriorityQueue<AStarNode>> neighbors = new ArrayList();
        neighbors.add(new PriorityQueue());
        neighbors.get(0).add(start);

        for (int index = 0; index < 10; index++) {
            PriorityQueue<AStarNode> rank = neighbors.get(neighbors.size() - 1);
            PriorityQueue<AStarNode> next = new PriorityQueue();
            for (AStarNode node : rank) {
                next.addAll(ne.enumerateNodes(node, codonFacade));
            }

            neighbors.add(next);
        }

        PriorityQueue<AStarNode> rank = neighbors.get(neighbors.size() - 1);

        int i = 0;
        for (AStarNode node : rank) {
            PartialResult result = HMMGraphSearch.partialResultFromGoal(node, true, false, 10, 1);
            //System.err.println(result.alignment);
        }
    }

    private static final Set<String> nuclExpected = new HashSet(Arrays.asList(new String[]{"-----",
        "-G---",
        "G----",
        "---G-",
        "--G--",
        "G-T--",
        "-GTT-",
        "GTT--",
        "GTTT-",
        "GT---",
        "-GT--",
        "GT---",
        "-G-T-",
        "--GT-",
        "gT---",
        "-GtT-",
        "gT-T-",
        "GtT--",
        "--G-T",
        "gTT--",
        "G-TT-",
        "gTTT-",
        "gttT-",
        "G--T-",
        "GttT-",
        "-G-TT",
        "-G--T",
        "GT-T-",
        "GTT-T",
        "G-T-T",
        "GtT-T",
        "gT--T",
        "----G",
        "GTtT-",
        "gTT-T",
        "gTtT-",
        "---Gt",
        "G-TtT",
        "G-TTT",
        "GtTT-",
        "G-TTt",
        "G---T",
        "gtT--",
        "gtTT-",
        "gTTtG",
        "GT--T",
        "-G-Tt",
        "G--TT",
        "-GT-T",
        "GtttG",
        "GT-TT",
        "GT-Tt",
        "-GTTT",
        "--GTt",
        "-GTTt",
        "--GTT",
        "-GttT",
        "-Gttt",
        "--Gtt",
        "--GtT",
        "-GTtT",
        "-GTtt",
        "-GtTT",
        "gT-TT",
        "-GtTt",
        "GTtTG",
        "gT-Tt",
        "GTtTg",
        "GTttG",
        "GTttg",
        "gTtTG",
        "---GT",
        "gTtTg",
        "gTttG",
        "gTttg",
        "GtTTG",
        "G-Ttt",
        "GtTTg",
        "gtT-T",
        "GTTTG",
        "GTTTg",
        "GTTtG",
        "GTTtg",
        "gTTTG",
        "gtTTG",
        "gTTTg",
        "GtTtG",
        "GtTtg",
        "gtTTg",
        "gTTtg",
        "G--Tt",
        "gttTG",
        "gttTg",
        "gtttG",
        "gtttg",
        "GttTG",
        "gtTtG",
        "GttTg",
        "gtTtg",
        "Gtttg"}));
}
