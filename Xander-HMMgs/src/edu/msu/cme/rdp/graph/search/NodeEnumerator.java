
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

import edu.msu.cme.rdp.alignment.hmm.MostProbableHCostHMM;
import edu.msu.cme.rdp.alignment.hmm.ProfileHMM;
import edu.msu.cme.rdp.graph.filter.NextCodon;
import static edu.msu.cme.rdp.alignment.hmm.TSC.*;
import edu.msu.cme.rdp.graph.filter.CodonWalker;
import edu.msu.cme.rdp.graph.search.heuristic.weight.HeuristicWeight;
import edu.msu.cme.rdp.kmer.Kmer;
import edu.msu.cme.rdp.readseq.SequenceType;
import edu.msu.cme.rdp.readseq.utils.NuclBinMapping;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author fishjord (edited by gilmanma)
 */
public class NodeEnumerator {

    private HeuristicWeight hweight;
    private AStarNode next;
    private NextCodon protEmission = null;
    private Byte nextNucl = null;
    private char emission;
    private double matchTrans;
    private double insTrans;
    private double delTrans;
    private Kmer nextKmer;
    private int nextState;
    private final ProfileHMM hmm;
    private final boolean protSearch;
    private final MostProbableHCostHMM hcost;
    private final long mask = (1L << 61) - 1;

    public NodeEnumerator(ProfileHMM hmm, HeuristicWeight hweight) {
        this.hmm = hmm;
        this.protSearch = hmm.getAlphabet() == SequenceType.Protein;
        this.hcost = hmm.getHCost();
        this.hweight = hweight;
    }

    public Set<AStarNode> enumerateNodes(AStarNode curr, CodonWalker walker) {
        return enumerateNodes(curr, walker, null);
    }
    
    /**
     * 
     * @param curr
     * @param walker
     * @param childNode if child node is not null, only return the node matching the child node
     * @return 
     */
    public Set<AStarNode> enumerateNodes(AStarNode curr, CodonWalker walker, AStarNode childNode) {
        Set<AStarNode> ret = new HashSet();

        nextState = curr.stateNo + 1;

        switch (curr.state) {
            case 'm':
                matchTrans = hmm.tsc(curr.stateNo, MM);
                insTrans = hmm.tsc(curr.stateNo, MI);
                delTrans = hmm.tsc(curr.stateNo, MD);
                break;
            case 'd':
                matchTrans = hmm.tsc(curr.stateNo, DM);
                insTrans = Double.NEGATIVE_INFINITY;
                delTrans = hmm.tsc(curr.stateNo, DD);
                break;
            case 'i':
                matchTrans = hmm.tsc(curr.stateNo, IM);
                insTrans = hmm.tsc(curr.stateNo, II);
                delTrans = Double.NEGATIVE_INFINITY;
                break;
            default:
                throw new RuntimeException("I hate you.");
        }

        walker.jumpTo(curr.kmer, curr.fwdHash, curr.rcHash);

        if (protSearch) {
            protEmission = walker.getNextCodon();
        } else {
            nextNucl = walker.getNextNucl();
        }

        double maxMatchEmission = hmm.getMaxMatchEmission(nextState);
        while ((protSearch && protEmission != null) || (!protSearch && nextNucl != null)) {
            if (protSearch) {
                if (protEmission.getAminoAcid() == '*') {
                    protEmission = walker.getSibCodon();
                    continue;
                }

                int codon = protEmission.getCodon() & 127;
                byte b1 = (byte) (codon & 0x3);
                byte b2 = (byte) (codon >> 2 & 0x3);
                byte b3 = (byte) (codon >> 4 & 0x3);
                nextKmer = curr.kmer.shiftLeft(b3).shiftLeft(b2).shiftLeft(b1);
                emission = protEmission.getAminoAcid();
            } else {
                nextKmer = curr.kmer.shiftLeft((byte) (nextNucl & 3));
                emission = NuclBinMapping.intToChar[nextNucl];
            }

            if ( childNode != null && !childNode.kmer.equals(nextKmer)){
                if (protSearch) {
                    protEmission = walker.getSibCodon();
                } else {
                    nextNucl = walker.getSibNucl();
                }
                continue;
            }
            
            final long fwdHash = walker.getFwdHash();
            final long rcHash = walker.getRcHash();

            /**
             * ************************************
             *
             * MATCH NODE
             *
             *************************************
             */
            next = new AStarNode(curr, nextKmer, fwdHash, rcHash, nextState, 'm');

            next.realScore = curr.realScore + matchTrans + hmm.msc(nextState, emission);
            if(next.realScore >= curr.maxScore) {
                next.maxScore = next.realScore;
                next.negativeCount = 0;
            } else {
                next.maxScore = curr.maxScore;
                next.negativeCount = curr.negativeCount + 1;
            }
            next.emission = emission;
            next.thisNodeScore = matchTrans + hmm.msc(nextState, emission) - maxMatchEmission;
            next.length = curr.length + 1;
            next.score = (curr.score + next.thisNodeScore);
            next.fval = (int) (HMMGraphSearch.INT_SCALE * (next.score + hweight.w(next) * hcost.computeHeuristicCost('m', nextState)));
            next.indels = curr.indels;

            if ( childNode!= null && childNode.equals(next)){
                ret.add(next);
                return ret;
            }else {
                ret.add(next);
            }

            /**
             * ************************************
             *
             * INSERT NODE
             *
             *************************************
             */
            if (curr.state != 'd') { //Transitions from delete to insert aren't allowed, don't waste time computing stuff
                next = new AStarNode(curr, nextKmer, fwdHash, rcHash, curr.stateNo /*
                         * Inserts don't advance the state
                         */, 'i');
                next.realScore = curr.realScore + insTrans + hmm.isc(nextState, emission);
                next.maxScore = curr.maxScore;
                next.negativeCount = curr.negativeCount + 1;
                next.emission = emission;
                next.thisNodeScore = insTrans + hmm.isc(nextState, emission);
                next.length = curr.length + 1;
                next.score = (curr.score + next.thisNodeScore);
                next.fval = (int) (HMMGraphSearch.INT_SCALE * (next.score + hweight.w(next) * hcost.computeHeuristicCost('i', curr.stateNo)));
                next.indels = curr.indels + 1;

                if ( childNode!= null && childNode.equals(next)){
                    ret.add(next);
                    return ret;
                }else {
                    ret.add(next);
                }
            }

            if (protSearch) {
                protEmission = walker.getSibCodon();
            } else {
                nextNucl = walker.getSibNucl();
            }
        }

        /**
         * ************************************
         *
         * DELETE NODE
         *
         *************************************
         */
        if (curr.state != 'i') {
            next = new AStarNode(curr, curr.kmer, curr.fwdHash, curr.rcHash, nextState, 'd');

            next.realScore = curr.realScore + delTrans;
            next.maxScore = curr.maxScore;
            next.negativeCount = curr.negativeCount + 1;
            next.emission = '-';
            next.thisNodeScore = delTrans - maxMatchEmission;
            next.length = curr.length;
            next.score = (curr.score + next.thisNodeScore);
            next.fval = (int) (HMMGraphSearch.INT_SCALE * (next.score + hweight.w(next) * hcost.computeHeuristicCost('d', nextState)));
            next.indels = curr.indels + 1;

            if ( childNode!= null && childNode.equals(next)){
                ret.add(next);
                return ret;
            }else {
                ret.add(next);
            }
        }

        return ret;
    }
}
