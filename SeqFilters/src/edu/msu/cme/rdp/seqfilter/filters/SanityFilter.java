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
public class SanityFilter implements SeqFilter {

    public enum SanityFilterExpectedSeqs {

        QSeq, Seq, Both
    };
    private SanityFilterExpectedSeqs expectedSeqType;

    public SanityFilter(SanityFilterExpectedSeqs expectedSeqType) {
        this.expectedSeqType = expectedSeqType;
    }

    public SeqFilterResult filterSequence(Sequence seq) {
        if (seq instanceof QSequence) {
            if (expectedSeqType == SanityFilterExpectedSeqs.Seq) {
                return new SeqFilterResult("Expecting Seqs only");
            }

            if (((QSequence) seq).getQuality().length != seq.getSeqString().length()) {
                return new SeqFilterResult("Sequence length [" + seq.getSeqString().length() + "] different from quality length [" + ((QSequence) seq).getQuality().length + "]");
            }
        } else {
            if (expectedSeqType == SanityFilterExpectedSeqs.QSeq) {
                return new SeqFilterResult("Expecting QSeqs only");
            }
        }
        return new SeqFilterResult(seq);
    }

    public String getName() {
        throw new UnsupportedOperationException("Sanity Filter");
    }
}
