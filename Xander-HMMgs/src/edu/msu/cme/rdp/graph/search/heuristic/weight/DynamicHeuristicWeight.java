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
 */
public class DynamicHeuristicWeight extends HeuristicWeight {
    
    public DynamicHeuristicWeight(double epsilon) {
        super(epsilon);
    }

    public DynamicHeuristicWeight(double epsilon, ProfileHMM hmm) {
        super(epsilon, hmm);
    }

    @Override
    public double w(AStarNode node) {
        return (1 + epsilon) * (1 - node.length/hmm.M());
    }
    
}
