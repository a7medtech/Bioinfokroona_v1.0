/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.msu.cme.rdp.seqfilter.output;

import edu.msu.cme.rdp.readseq.readers.Sequence;

/**
 *
 * @author fishjord
 */
public interface SeqFilterOutput {
    public void appendSequence(Sequence s);
}
