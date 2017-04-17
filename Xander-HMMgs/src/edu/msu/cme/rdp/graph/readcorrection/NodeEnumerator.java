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

import edu.msu.cme.rdp.graph.search.*;
import java.util.Arrays;
import edu.msu.cme.rdp.alignment.pairwise.ScoringMatrix;
import edu.msu.cme.rdp.graph.filter.CodonWalker;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author fishjord
 */
public class NodeEnumerator {

    private static final ScoringMatrix sm = ScoringMatrix.getDefaultNuclMatrix();

    /*static {
        try {
            sm = new ScoringMatrix(ScoringMatrix.class.getResourceAsStream("/data/NUC.4.4"), -100, -1);
        } catch(IOException e) {
            throw new RuntimeException(e);
        }
    }*/

    private int bestMatch = 5;
    private AStarNode next;
    private char[] nextKmer;
    private int index;
    private final boolean forward;
    private final char[] seq;
    private Character nextNucl;

    public NodeEnumerator(boolean forward, char[] seq) {
        this.forward = forward;
        this.seq = seq;
    }

    public Set<AStarNode> enumerateNodes(AStarNode curr, CodonWalker walker) {
        Set<AStarNode> ret = new HashSet();

        //walker.jumpTo(curr.kmer);
        //nextNucl = walker.getNextNucl();

        int remaining;
        if (forward) {
            index = curr.stateNo + 1;
            remaining = seq.length - index;
        } else {
            index = curr.stateNo - 1;
            remaining = index;
        }

        if (nextNucl != null) {

            while (nextNucl != null) {
                nextKmer = makeNextKmer(forward, nextNucl, curr);

                /**
                 * ************************************
                 *
                 * MATCH NODE
                 *
                 *************************************
                 */
                //next = new AStarSearchNode(curr, new char[]{nextNucl}, index, 'm', nextKmer);

                next.score = curr.score + sm.score(nextNucl, seq[index]);
                next.fval = (int) (HMMGraphSearch.INT_SCALE * (next.score + remaining * bestMatch));

                ret.add(next);

                /**
                 * ************************************
                 *
                 * Insert NODE
                 *
                 *************************************
                 */
                //next = new AStarSearchNode(curr, new char[]{nextNucl}, index, 'i', nextKmer);

                next.score = curr.score + sm.getGapOpen();
                next.fval = (int) (HMMGraphSearch.INT_SCALE * (next.score + remaining * bestMatch));

                ret.add(next);

                //nextNucl = walker.getSibNucl();
            }

            /**
             * ************************************
             *
             * Delete NODE
             *
             *************************************
             */
            //next = new AStarSearchNode(curr, new char[]{'-'}, curr.stateNo, 'd', curr.kmer);

            next.score = curr.score + sm.getGapOpen();
            next.fval = (int) (HMMGraphSearch.INT_SCALE * (next.score + remaining * bestMatch));

            ret.add(next);
        }

        return ret;
    }

    private static char[] makeNextKmer(boolean forward, char b, AStarNode curr) {
        return null;
        /*char[] kmer;

        if (forward) {
            kmer = Arrays.copyOfRange(curr.kmer, 1, curr.kmer.length + 1);
            kmer[kmer.length - 1] = b;
        } else {
            kmer = new char[curr.kmer.length];

            kmer[0] = b;

            for (int i = 1; i < kmer.length; i++) {
                kmer[i] = curr.kmer[i - 1];
            }
        }

        return kmer;*/
    }
}
