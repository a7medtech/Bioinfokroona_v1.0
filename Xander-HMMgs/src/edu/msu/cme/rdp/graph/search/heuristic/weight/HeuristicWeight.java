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
 * The type of weighting to use when calculating the heuristic
 * for a weighted A* search algorithm.
 * 
 * @author gilmanma
 * 
 */
public abstract class HeuristicWeight {
    
    protected double epsilon;
    protected ProfileHMM hmm;
    
    /**
     * Class constructor.
     * 
     * @param epsilon   the "base" weight; algorithms should be designed such
     *                  that the solution is no worse than (1 + epsilon) times
     *                  the optimal one
     */
    public HeuristicWeight(double epsilon) {
        this(epsilon, null);
    }
    
    /**
     * 
     * Class constructor that allows the specification of a profile hmm for
     * more advanced functionality.
     * 
     * @param epsilon   the "base" weight; algorithms should be designed such
     *                  that the solution is no worse than (1 + epsilon) times
     *                  the optimal one
     * @param hmm       a ProfileHMM from RDPTool's AlignmentTools 
     */
    public HeuristicWeight(double epsilon, ProfileHMM hmm) {
        this.epsilon = epsilon;
        this.hmm = hmm;
    }
    
    /**
     * Allows for the specification of the profile hmm after the instance has
     * been created.
     * 
     * @param hmm   a ProfileHMM from RDPTool's AlignmentTools 
     */
    public void setHMM(ProfileHMM hmm) {
        this.hmm = hmm;
    }
    
    /**
     * The "meat" of the class, this calculates and returns the amount by
     * which to weight the heuristic based on the current node, epsilon, and
     * the profile hmm.
     * 
     * @param node  current AStarNode to find the heuristic cost of
     * @return      amount by which to multiply the heuristic cost
     */
    public abstract double w(AStarNode node);
    
}
