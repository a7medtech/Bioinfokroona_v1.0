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
package edu.msu.cme.rdp.graph.readcorrection;

/**
 *
 * @author fishjord
 */
public class ReadCorrector {

    /*
    private static ScoringMatrix sm = ScoringMatrix.getDefaultNuclMatrix();

    private static class CandidateAlignment implements Comparable<CandidateAlignment> {

        private String alignment;
        private int score;

        public int compareTo(CandidateAlignment o) {
            return o.score - score;
        }
    }

    private static CandidateAlignment align(BloomFilter filter, char[] seq, char[] kmer, int index) {
        AStarSearchNode leftStart = new AStarSearchNode(null, new char[]{kmer[0]}, index, 'm', kmer);
        AStarSearchNode rightStart = new AStarSearchNode(null, new char[]{kmer[kmer.length - 1]}, index + kmer.length - 1, 'm', kmer);

        //System.out.println("Generating left");
        AStarSearchNode leftGoal = subalign(filter.new LeftCodonFacade(kmer, 0), leftStart, new NodeEnumerator(false, seq), seq.length);
        //System.out.println("Left generated");
        //System.out.println("Generating right");
        AStarSearchNode rightGoal = subalign(filter.new RightCodonFacade(kmer, 0), rightStart, new NodeEnumerator(true, seq), seq.length);
        //System.out.println("Right generated");

        if (leftGoal == null || rightGoal == null) {
            return null;
        }

        String align = extractString(leftGoal, false) + new String(kmer) + extractString(rightGoal, true);
        int score = 0;
        for (int i = 0; i < align.length(); i++) {
            score += sm.score(align.charAt(i), seq[i]);
        }

        CandidateAlignment aln = new CandidateAlignment();
        aln.alignment = align;
        aln.score = score;

        return aln;
    }

    private static AStarSearchNode subalign(CodonWalker walker, AStarSearchNode startingVertex, NodeEnumerator enumerator, int seqLength) {
        HashHeap<AStarSearchNode> open = new HashHeap();
        Set<AStarSearchNode> closed = new HashSet();
        open.add(startingVertex);
        int max = 0;

        AStarSearchNode curr;
        while (!open.isEmpty()) {
            curr = open.pop();
            closed.add(curr);
            /*if(curr.stateNo > max) {
                max = curr.stateNo;
                System.out.println(closed.size() + "\t" + open.size() + "\t" + curr.deletes + "\t" + curr);
            }*//*

            if (curr.stateNo == seqLength - 1 || curr.stateNo == 0) {
                return curr;
            }

            for (AStarSearchNode next : enumerator.enumerateNodes(curr, walker)) {
                if (!closed.contains(next)
                        && (!open.contains(next) || next.score > open.get(next).score)) {
                    open.add(next);
                }
            }
        }

        return null;
    }

    private static String extractString(AStarSearchNode goal, boolean forward) {
        StringBuilder ret = new StringBuilder();

        while (goal.discoveredFrom != null) {  //Ignore the first node, as it is part of the kmer
            char b = goal.emission[0];

            if (goal.state == 'm') {
                b = Character.toUpperCase(b);
            } else if (goal.state == 'i') {
                b = Character.toLowerCase(b);
            }

            if (forward) {
                ret.insert(0, b);
            } else {
                ret.append(b);
            }

            goal = goal.discoveredFrom;
        }

        return ret.toString();
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 2) {
            System.err.println("USAGE: ReadCorrector <bloom_filter> <reads>");
            System.exit(1);
        }

        File bloomFile = new File(args[0]);
        File readsFile = new File(args[1]);
        File correctedFile = new File("corrected_" + readsFile.getName());

        BloomFilter bloom;

        ObjectInputStream ois = new ObjectInputStream(new BufferedInputStream(new FileInputStream(bloomFile)));
        bloom = (BloomFilter) ois.readObject();
        ois.close();

        SeqReader reader = new SequenceReader(readsFile);
        Sequence seq;

        FastaWriter out = new FastaWriter(correctedFile);
        int k = bloom.getKmerSize();
        GraphState gs = bloom.new GraphState();

        HashHeap<CandidateAlignment> alignments = new HashHeap();

        while ((seq = reader.readNextSequence()) != null) {
            char[] bases = seq.getSeqString().toCharArray();

            for (int index = 0; index < (bases.length - k + 1); index++) {
                char[] kmer = Arrays.copyOfRange(bases, index, index + k);

                gs.setState(kmer);

                if (gs.hasCurrent()) {
                    CandidateAlignment aln = align(bloom, bases, kmer, index);
                    if (aln != null) {
                        System.err.println(seq.getSeqName() + "\t" + index + "\t" + aln.score);
                        alignments.add(aln);
                    } else {
                        System.err.println(seq.getSeqName() + "\t" + index + "\t-");
                    }
                }
            }

            CandidateAlignment best = alignments.top();
            if (best != null) {
                System.out.println(seq.getSeqName() + "\t" + best.score);
                out.writeSeq(seq.getSeqName(), best.alignment);
            } else {
                System.out.println(seq.getSeqName() + "\t(none)");
            }
        }

        out.close();
    }*/
}
