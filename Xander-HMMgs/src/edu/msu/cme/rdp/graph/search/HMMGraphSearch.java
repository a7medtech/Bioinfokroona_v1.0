/*
 * Copyright (C) 2012 Michigan State University <rdpstaff at msu.edu>
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
import edu.msu.cme.rdp.alignment.hmm.scoring.ForwardScorer;
import edu.msu.cme.rdp.alignment.hmm.scoring.HMMScorer;
import edu.msu.cme.rdp.graph.filter.CodonWalker;
import edu.msu.cme.rdp.graph.filter.InvalidDNABaseException;
import edu.msu.cme.rdp.graph.search.heuristic.weight.HeuristicWeight;
import edu.msu.cme.rdp.kmer.Kmer;
import edu.msu.cme.rdp.kmer.NuclKmer;
import edu.msu.cme.rdp.readseq.SequenceType;
import edu.msu.cme.rdp.readseq.utils.ProteinUtils;
import java.io.IOException;
import java.util.*;

/**
 *
 * @author fishjord (edited by gilmanma)
 */
public class HMMGraphSearch {

    public static class PartialResult {

        public String maxSeq;
        public String alignment;
        public double maxScore;
        public long searchTime;
    }

    public static class HackTerminateException extends RuntimeException {
    }
    public static final int PRUNE_NODE = 20;  // prune the search if the score does not improve after p number of nodes
    public static final int INT_SCALE = 10000; //This is the number of sigfigs in a HMMER3 model, so it works out quite nicely if you ask me
    private static final int upperBound = Integer.MIN_VALUE;
    private final int maxk;  // max number of shortest paths to search
    
    
    /**
     * Map of previously discovered terminal nodes.  Used to shorten search times
     * when algorithm encounters nodes from a previously cached path.  
     */
    private Map<AStarNode, AStarNode> termNodes = new HashMap<AStarNode, AStarNode>();
    
    /**
     * prune the search if the score does not improve after p number of nodes
     * Set to -1 to deactivate.
     */
    private int heuristicPruning = PRUNE_NODE;
    
    /**
     * Method to use when weighting heuristic
     */
    private HeuristicWeight hweight;

    public HMMGraphSearch(int maxk, int n_nodes) {
        this.maxk = maxk;
        this.heuristicPruning = n_nodes;
    }
    
    
    /**
     * Set which method should be used to weight the heuristic during f-score
     * calculation
     * 
     * @param hweight   weighting method
     */
    public void setHWeight(HeuristicWeight hweight) {
        this.hweight = hweight;
    }

    public List<SearchResult> search(SearchTarget target) throws InterruptedException {
        String framedKmer = target.getKmer();
        List<SearchResult> ret = new ArrayList();

        System.err.println("left " );
        hweight.setHMM(target.getReverseHmm());
        int lStartingState = target.getReverseHmm().M() - target.getStartState() - target.getKmer().length() / ((target.isProt()) ? 3 : 1);

        try {
            List<PartialResult> leftParts = kpathsSearch(target.getReverseHmm(), lStartingState, framedKmer, target.getFilter().new LeftCodonFacade(target.getKmer()), false);
            for (PartialResult r : leftParts) {
                String nuclSeq = r.maxSeq + framedKmer;
                String alignment = r.alignment + framedKmer.toUpperCase();
                String protSeq = null;
                String scoringSeq = nuclSeq;

                HMMScorer scorer = new ForwardScorer(target.getForwardHmm(), -1);

                if (target.isProt()) {
                    scoringSeq = protSeq = ProteinUtils.getInstance().translateToProtein(nuclSeq, true, 11);
                }

                for (char c : scoringSeq.toCharArray()) {
                    scorer.consume(c);
                }

                ret.add(new SearchResult(target, target.getKmer(), nuclSeq, alignment, protSeq, SearchResult.SearchDirection.left, lStartingState, r.maxScore, scorer.getMaxScore(), r.searchTime));
            }

            System.err.println("right");
            hweight.setHMM(target.getForwardHmm());
            List<PartialResult> rightParts = kpathsSearch(target.getForwardHmm(), target.getStartState(), framedKmer, target.getFilter().new RightCodonFacade(target.getKmer()), true);

            for (PartialResult r : rightParts) {
                String nuclSeq = framedKmer + r.maxSeq;
                String alignment = framedKmer.toUpperCase() + r.alignment;
                String protSeq = null;
                String scoringSeq = nuclSeq;

                HMMScorer scorer = new ForwardScorer(target.getForwardHmm(), -1);

                if (target.isProt()) {
                    scoringSeq = protSeq = ProteinUtils.getInstance().translateToProtein(nuclSeq, true, 11);
                }

                for (char c : scoringSeq.toCharArray()) {
                    scorer.consume(c);
                }

                ret.add(new SearchResult(target, target.getKmer(), nuclSeq, alignment, protSeq, SearchResult.SearchDirection.right, target.getStartState(), r.maxScore, scorer.getMaxScore(), r.searchTime));
            }
        } catch(InvalidDNABaseException ex){
            System.err.println("Warning: " + ex.getMessage());
        }

        return ret;
    }

    public List<AStarNode> searchGraph(SearchTarget target) throws InterruptedException {
        String framedKmer = target.getKmer();
        List<AStarNode> ret = new ArrayList();

        int lStartingState = target.getReverseHmm().M() - target.getStartState() - target.getKmer().length() / ((target.isProt()) ? 3 : 1);

        List<CandidatePath> leftParts = kpathsSearchGraph(target.getReverseHmm(), lStartingState, framedKmer, target.getFilter().new LeftCodonFacade(target.getKmer()), false);

        for (CandidatePath r : leftParts) {
            ret.add(r.get(r.length() - 1));
        }

        List<CandidatePath> rightParts = kpathsSearchGraph(target.getForwardHmm(), target.getStartState(), framedKmer, target.getFilter().new RightCodonFacade(target.getKmer()), true);

        for (CandidatePath r : rightParts) {
            ret.add(r.get(r.length() - 1));
        }

        return ret;
    }

    private List<PartialResult> kpathsSearch(ProfileHMM hmm, int startingState, String framedWord, CodonWalker walker, boolean forward) throws InterruptedException {
        List<PartialResult> results = new ArrayList();
        for (CandidatePath path : kpathsSearchGraph(hmm, startingState, framedWord, walker, forward)) {
            AStarNode goal = path.get(path.length() - 1);
            results.add(partialResultFromGoal(goal, forward, hmm.getAlphabet() == SequenceType.Protein, framedWord.length(), path.generationTime));
        }

        return results;
    }

    private List<CandidatePath> kpathsSearchGraph(ProfileHMM hmm, int startingState, String framedWord, CodonWalker walker, boolean forward) throws InterruptedException {
        List<CandidatePath> bestPaths = new ArrayList();  // best path to return, each contains at least one new kmer, different from bestBasePaths
        CandidatePath bestBasePath = null;   // best path to be used to find base
        PriorityQueue<CandidatePath> candidatePaths = new PriorityQueue<CandidatePath>();
        Map<AStarNode, Set<AStarNode>> shortestPathEdges = new HashMap();
        Set<Kmer> seenKmers = new HashSet(); // keep a set of kmers seen from the previous chosen best paths

        long kTime = System.currentTimeMillis();
        try {
            AStarNode goalNode = astarSearch(hmm, startingState, framedWord, walker, forward,  new HashSet(), true);

            CandidatePath bestPath = new CandidatePath(goalNode);
            bestPath.generationTime = (System.currentTimeMillis() - kTime);
            bestPaths.add(bestPath);
            bestBasePath = bestPath;
            // add seen kmers
            for ( int n = 0; n < bestPath.length(); n++){
                seenKmers.add(bestPath.get(n).kmer);
            }
                
            for(int index = 0;index < (bestPath.path.size()-1);index++) {
                AStarNode node = bestPath.path.get(index);
                termNodes.put(node, bestPath.path.get(index+1));
            }		            
            while (bestPaths.size() < maxk) {   //Where k is the current kth shortest path
                CandidatePath pathAk = bestBasePath;
                kTime = System.currentTimeMillis();
               
                /*
                 * Candidate generation
                 */
                for (int i = pathAk.i; i < pathAk.length() - 1; i++) { //Lawler's Observation
                    CandidatePath base = pathAk.subpath(i + 1);
                    AStarNode starting = base.get(i);
                    AStarNode ak_i_1 = pathAk.get(i + 1);
                    
                    //If edge between starting and ak-i-1 was previously blocked to find a spur
                    //This is to avoid the combinatorial of multiple paths
                    if (shortestPathEdges.containsKey(starting) && shortestPathEdges.get(starting).contains(ak_i_1)) {
                        break;
                    }

                    if (!shortestPathEdges.containsKey(starting)) {
                        shortestPathEdges.put(starting, new HashSet());
                    }
                    shortestPathEdges.get(starting).add(ak_i_1);
                    goalNode = astarSearch(hmm, starting, walker, shortestPathEdges.get(starting), false);
                    if(goalNode == null) {
                        break;
                    }

                    CandidatePath candidate = new CandidatePath(goalNode);
                    candidate.i = i;
                    candidate.k = bestPaths.size();

                    //System.err.println(bestPaths.size() + " find goalnode " + goalNode.kmer + " " + goalNode.stateNo + " " + goalNode.score + " realscore " +  goalNode.realScore);
                    //System.err.println("\nstarting i " + i + " " + starting.kmer + " " + starting.stateNo + " ak_i_1 " + ak_i_1.kmer + " " + ak_i_1.stateNo);
                    
                    if (!bestPaths.contains(candidate) && !Double.isInfinite(candidate.score) ) { 
                        candidatePaths.add(candidate);                        
                        /* This is very impportant: only save the edges starting from the ak_i_1 node in the path found */
                        
                        for( int startIndex = candidate.i+1 ;startIndex < candidate.path.size() -1;startIndex++) {
                            AStarNode node = candidate.path.get(startIndex);                            
                            termNodes.put(node,  candidate.path.get(startIndex+1));                           
                        }		
                    } else { // need to figure out if this is needed
                        // shortestPathEdges.get(starting).remove(ak_i_1);
                    }
                    
                }

                CandidatePath kthPath = candidatePaths.poll();
                if (kthPath == null || Double.isInfinite(kthPath.score)) {  
                    break;
                }
                bestBasePath = kthPath;  // we need added to best path the set to find the bases
                // but we only keep that path in the returned list if it contains at least one new kmer
                Set<Kmer> newKmers = new HashSet();
                for ( int n = 0; n < kthPath.length(); n++){
                    newKmers.add(kthPath.get(n).kmer);
                }
                newKmers.removeAll(seenKmers);
                if ( ! newKmers.isEmpty()) {
                    kthPath.generationTime = (System.currentTimeMillis() - kTime);
                    bestPaths.add(kthPath); 
                    seenKmers.addAll(newKmers);
                }  
            }
        } catch (HackTerminateException e) {
            System.err.println("Terminated on path " + bestPaths.size() + " (candidates=" + candidatePaths.size() + ")");            
            throw e;
        } catch (IOException ignore) {
            throw new RuntimeException(ignore);
        }
        return bestPaths;
    }

    /**
     * initializes the starting node and starts the main search 
     * (NOTE: Walker MUST be initialized to the passed starting kmer)
     *
     * @param hmm           ProfileHMM for the target gene
     * @param startingState 
     * @param framedWord
     * @param walker
     * @param forward
     * @param disallowedLinks
     * @return
     * @throws IOException
     */
    private AStarNode astarSearch(final ProfileHMM hmm,
            int startingState, String framedWord,
            CodonWalker walker,
            boolean forward,
            Set<AStarNode> disallowedLinks,
            boolean bestOnlySearch) throws IOException, InterruptedException {
        framedWord = framedWord.toLowerCase();

        char[] startingCodon = framedWord.substring(framedWord.length() - 3).toCharArray();

        String scoringWord = framedWord;
        if (!forward) {
            if (hmm.getAlphabet() == SequenceType.Protein) {
                StringBuilder tmp = new StringBuilder(ProteinUtils.getInstance().translateToProtein(framedWord, false, 11));
                scoringWord = tmp.reverse().toString();
            } //else {
               // scoringWord = framedWord;
            //}
            framedWord = new StringBuilder(framedWord).reverse().toString(); 
        } else if (hmm.getAlphabet() == SequenceType.Protein) {
            scoringWord = ProteinUtils.getInstance().translateToProtein(framedWord, false, 11);
        }

        Kmer kmer = new NuclKmer(framedWord.toCharArray());

        AStarNode startingNode;
        if (hmm.getAlphabet() == SequenceType.Protein) {
            startingNode = new AStarNode(null, kmer, walker.getFwdHash(), walker.getRcHash(), startingState + (framedWord.length() / 3), 'm');
            startingNode.length = framedWord.length() / 3;
        } else {
            startingNode = new AStarNode(null, kmer, walker.getFwdHash(), walker.getRcHash(), startingState, 'm');
            startingNode.length = framedWord.length();
        }

        startingNode.fval = 0;
        startingNode.score = scoreStart(hmm, scoringWord, startingState);
        startingNode.realScore = realScoreStart(hmm, scoringWord, startingState);

        return astarSearch(hmm, startingNode, walker,  disallowedLinks, bestOnlySearch);
    }

    private float scoreStart(ProfileHMM hmm, String startingKmer, int startingState) {
        //System.err.println("Starting kmer: " + startingKmer + ", state: " + startingState);
        float ret = 0;
        char[] residues = startingKmer.toCharArray();
        for (int index = 1; index <= residues.length; index++) {
            ret += hmm.msc(startingState + index, residues[index - 1]) + hmm.tsc(startingState + index, TSC.MM) - hmm.getMaxMatchEmission(startingState + index);
        }

        return ret;
    }

    private float realScoreStart(ProfileHMM hmm, String startingKmer, int startingState) {
        //System.err.println("Starting kmer: " + startingKmer + ", state: " + startingState);
        float ret = 0;
        char[] residues = startingKmer.toCharArray();
        for (int index = 1; index <= residues.length; index++) {
            ret += hmm.msc(startingState + index, residues[index - 1]) + hmm.tsc(startingState + index - 1, TSC.MM);
            //System.err.println(index+"\t"+hmm.msc(startingState + index, residues[index - 1])+"\t"+hmm.tsc(startingState + index - 1, TSC.MM));
            //System.err.println(index + "\t" + residues[index - 1] + "\t" + ret + "\t" + (ret + Math.log(2.0 / (index + startingState + 2)) * 2 - HMMScorer.getNull1(index + startingState)) / HMMScorer.ln2);
        }

        return ret;
    }
    // assume the longest protein we care is 3000 aa.
    private static double[] exitProbabilities = new double[3000];

    static {
        for (int index = 0; index < exitProbabilities.length; index++) {
            exitProbabilities[index] = Math.log(2.0 / (index + 2)) * 2;
        }
    }
    private static final double ln2 = Math.log(2);

    /**
     * The core algorithm of the search.  Finds the path from starting kmer to the end
     * of the hmm.
     * 
     * @param hmm               ProfileHMM for target gene
     * @param startingNode      node at which to start the search
     * @param walker            used to find next node(s) to open
     * @param disallowedLinks   node transitions that the algorithm should not take
     * @return final node in path
     * @throws IOException
     * @throws InterruptedException 
     */
    private AStarNode astarSearch(final ProfileHMM hmm,
            AStarNode startingNode,
            CodonWalker walker,           
            Set<AStarNode> disallowedLinks,
             boolean bestOnlySearch ) throws IOException, InterruptedException {

        // if the starting node is already at the end of the model, return it and exit
        if (startingNode.stateNo >= hmm.M()) {
            System.err.println(startingNode.kmer + "\t-\t-\t-\t-\t-\t" + false);
            return startingNode;
        }
        
        NodeEnumerator nodeEnumerator = new NodeEnumerator(hmm, hweight);
        PriorityQueue<AStarNode> open = new PriorityQueue<AStarNode>();
        Set<AStarNode> closed = new HashSet();
        AStarNode curr;
        // maximum number of insertions and deletions in the path
        int openedNodes = 1;
        
        Map<AStarNode, AStarNode> openHash = new HashMap<AStarNode, AStarNode>();
        
        int repeatedNodes = 0;
        int replacedNodes = 0;
        int prunedNodes = 0;

        //First step, enumerate all the nodes and remove any disallowed transitions
        //This way we only have to look at the set (disallowedLinks) once 
        // And if this is for the best path search, we only need to open the child node from the saved best edge
        if ( bestOnlySearch && termNodes.containsKey(startingNode)){            
            AStarNode childtoopen = termNodes.get(startingNode);
            for (AStarNode next : nodeEnumerator.enumerateNodes(startingNode, walker, childtoopen)) { 
                open.add(next);
            }
        }else {
            for (AStarNode next : nodeEnumerator.enumerateNodes(startingNode, walker)) {                
                if (!disallowedLinks.contains(next)) {                    
                    open.add(next);
                }
            }
        }
        
        if (open.isEmpty()) {
            return null;
        }
        //Decide the intermediate goal
        AStarNode interGoal = startingNode;

        //While we have more things to close
        while ((curr = open.poll()) != null) {
            if (closed.contains(curr)) { // we may need to examine the nodes in the closed set because of the bounded relaxation heuristic used
                continue;
            }

            if (curr.stateNo >= hmm.M()) { //We're at an "end" state
                curr.partial = false;

                System.err.println(startingNode.kmer + "\t" + openedNodes + "\t" + open.size() + "\t" + closed.size() + "\t" + repeatedNodes + "\t" + replacedNodes + "\t" + prunedNodes + "\t"+ false);
                // compare to the partial goal nodes to see which one has higher score                                   
                if ((curr.realScore + exitProbabilities[curr.length] - HMMScorer.getNull1(curr.length)) / ln2
                        > (interGoal.realScore + exitProbabilities[interGoal.length] - HMMScorer.getNull1(interGoal.length)) / ln2) {
                //if (curr.realScore > interGoal.realScore) {    
                    interGoal = curr;
                } 
                
                // find the highest scoring node along the path as the returning goal node  
                return getHighestScoreNode(interGoal);                   
            }

            closed.add(curr);            
            
            if ((curr.realScore + exitProbabilities[curr.length] - HMMScorer.getNull1(curr.length)) / ln2
                    > (interGoal.realScore + exitProbabilities[interGoal.length] - HMMScorer.getNull1(interGoal.length)) / ln2) {
                interGoal = curr;
            }           

            Set<AStarNode> tempNodesToOpen = null;
            //if the node is already on a previous shortest path, we just need to open the next child node from the saved edge            
            if(termNodes.containsKey(curr)) {  
                AStarNode childtoopen = termNodes.get(curr);
                tempNodesToOpen = nodeEnumerator.enumerateNodes(curr, walker, childtoopen);               
            }else {
                tempNodesToOpen = nodeEnumerator.enumerateNodes(curr, walker);
            }
            
            //Look at the adjacent nodes
            for (AStarNode next : tempNodesToOpen) {
                boolean openNode = false;
                if(heuristicPruning > 0) {
                    // Don't open nodes that fail to pass this heuristic
                    
                    if((next.length < 5 || next.negativeCount <= heuristicPruning) && next.realScore > 0.0) {
                        // check to see if this node has been opened previously
                        AStarNode old = openHash.get(next);
                        if(old != null) {
                            // if it has, only open it again if its path is better than the old copy's
                            repeatedNodes++;
                            if (old.compareTo(next) > 0) {
                                replacedNodes++;
                                // remove a node from heap is O(n) operation, we don't need to do this because it does not hurt to have multiple entries in the heap
                                //open.remove(next);
                                //open.remove(old);
                                openNode = true;
                            }
                        } else {
                            openNode = true;
                        }
                    } else {
                        prunedNodes++;
                    }
                } else {
                    AStarNode old = openHash.get(next);
                    if(old != null) {
                        repeatedNodes++;
                        if (old.compareTo(next) > 0) {
                            replacedNodes++;
                            //open.remove(next);
                            //open.remove(old);
                            openNode = true;
                        }
                    } else {
                        openNode = true;
                    }
                }
                
                if(openNode) {
                    openHash.put(next, next);
                    openedNodes++;
                    open.add(next);
                }
            }
        }
        
        // if the loop exits, then the search failed to reach the end state
        System.err.println(startingNode.kmer + "\t" + openedNodes + "\t" + open.size() + "\t" + closed.size() + "\t" + repeatedNodes + "\t" + replacedNodes + "\t" + prunedNodes + "\t"+  true);
        interGoal.partial = true;
        return getHighestScoreNode(interGoal);
    }
    
    
    private AStarNode getHighestScoreNode(AStarNode node){
        // need to find the highest scoring node along the path as the returning goal node
        // because sometimes the nodes after the highest scoring node are junk
        AStarNode tempGoal = node;
        AStarNode highestScoreNode = node;
        while ( tempGoal.discoveredFrom != null){
            tempGoal = tempGoal.discoveredFrom;
            if ( tempGoal.realScore > highestScoreNode.realScore){
                highestScoreNode = tempGoal;
            }
        }
        
        return highestScoreNode;                   
    }
    
    /**
     * Assemble the path, from the goal to the start
     */
    public static PartialResult partialResultFromGoal(AStarNode goal, boolean forward, boolean protSearch, int kmerLength, long searchTime) {
        StringBuilder nuclSeq = new StringBuilder();
        StringBuilder alignmentSeq = new StringBuilder();

        char[] gap = (protSearch) ? new char[]{'-', '-', '-'} : new char[]{'-'};

        PartialResult result = new PartialResult();
        if (goal != null) {
            result.maxScore = goal.score;

            while (goal.discoveredFrom != null) {
                char[] kmer = goal.kmer.toString().toCharArray();
                char[] emission;

                if (protSearch) {
                    if (forward) {
                        //If we're to the right the codon is at the end of the kmer
                        emission = new char[]{kmer[kmer.length - 3], kmer[kmer.length - 2], kmer[kmer.length - 1]};
                    } else {
                        //If we're going to the left the codon is still at the right
                        //end of this kmer BUT it is in reverse order (DIFFERENT
                        //than in the path returned by CodonWalker.getPathString())
                        emission = new char[]{kmer[kmer.length - 1], kmer[kmer.length - 2], kmer[kmer.length - 3]};
                    }
                } else {
                    //In the single emission case the last character in the kmer
                    //is always right
                    emission = new char[]{kmer[kmer.length - 1]};
                }

                if (goal.state == 'd') {
                    if (forward) {
                        alignmentSeq.insert(0, gap);
                    } else {
                        alignmentSeq.append(gap);
                    }
                } else if (goal.state == 'm') {
                    if (forward) {
                        alignmentSeq.insert(0, new String(emission).toUpperCase());
                    } else {
                        alignmentSeq.append(new String(emission).toUpperCase());
                    }
                } else if (goal.state == 'i') {
                    if (forward) {
                        alignmentSeq.insert(0, new String(emission).toLowerCase());
                    } else {
                        alignmentSeq.append(new String(emission).toLowerCase());
                    }
                }


                if (goal.state != 'd') { //No emission on delete states
                    if (forward) {
                        //prepend for forward
                        nuclSeq.insert(0, emission);
                    } else {
                        //append for reverse (we're building in the 'right' direction)
                        nuclSeq.append(emission);
                    }
                }

                goal = goal.discoveredFrom;
            }
        }
        result.maxSeq = nuclSeq.toString();
        result.alignment = alignmentSeq.toString();
        result.searchTime = searchTime;

        return result;
    }

    private static double getAverageScore(AStarNode node, int over) {
        double ret = 0;
        int index;
        for (index = 0; index < over && node != null; index++, node = node.discoveredFrom) {
            ret += node.thisNodeScore;
        }

        if (index != over) {
            throw new IllegalArgumentException("Path does not go back " + over + " nodes");
        }

        return ret / over;
    }

    public static double getMemRatio() {
        return (double) (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / Runtime.getRuntime().maxMemory();
    }
}
