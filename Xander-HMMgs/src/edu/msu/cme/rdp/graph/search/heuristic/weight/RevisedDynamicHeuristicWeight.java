/*
 * Copyright (C) 2014 gilmanma
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

package edu.msu.cme.rdp.graph.search.heuristic.weight;

import edu.msu.cme.rdp.alignment.hmm.ProfileHMM;
import edu.msu.cme.rdp.graph.search.AStarNode;

/**
 *
 * @author gilmanma
 * 
 * Based on algorithm given in
 *  Thayer, J., Ruml W., and Kreis J. 2009. Using Distance Estimates in 
 *      Heuristic Search: A Re-evaluation. Symposium on Combinatorial Search.
 */
public class RevisedDynamicHeuristicWeight extends HeuristicWeight {
    
    public RevisedDynamicHeuristicWeight(double epsilon) {
        super(epsilon);
    }

    public RevisedDynamicHeuristicWeight(double epsilon, ProfileHMM hmm) {
        super(epsilon, hmm);
    }

    @Override
    public double w(AStarNode node) {
        double dn = hmm.M() - node.length;
        return Math.max(1.0, (1.0+epsilon) * dn / hmm.M());
    }
    
}
