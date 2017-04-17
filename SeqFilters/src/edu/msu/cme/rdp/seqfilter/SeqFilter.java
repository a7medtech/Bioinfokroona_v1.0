/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.msu.cme.rdp.seqfilter;

import edu.msu.cme.rdp.readseq.readers.Sequence;

/**
 *
 * @author fishjord
 */
public interface SeqFilter {
    /**
     * @param seq sequence to be filtered
     * @return sequence after filtering or null if the sequence fails the filter
     */
    public SeqFilterResult filterSequence(Sequence seq);
    
    public String getName();
}
