/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.msu.cme.rdp.seqfilter.filters;

import edu.msu.cme.rdp.readseq.QSequence;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.seqfilter.SeqFilter;
import edu.msu.cme.rdp.seqfilter.SeqFilterResult;

/**
 *
 * @author fishjord
 */
public class NSeqFilter implements SeqFilter {
    private int maxNs;
    
    public NSeqFilter(int maxNs) {
        this.maxNs = maxNs;
    }

    public String getName() {
        return "N count > " + maxNs + " seq filter";
    }

    public SeqFilterResult filterSequence(Sequence s) {
        int numNs = s.getSeqString().replaceAll("[^Nn]", "").length();
        if(numNs > maxNs)
            return new SeqFilterResult(numNs + " is greater than the maximum allowed Ns " + maxNs);

        return new SeqFilterResult(s);
    }

    public SeqFilterResult filterSequence(QSequence s) {
        return filterSequence((Sequence)s);
    }
}
