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
public class CaseTransformFilter implements SeqFilter {

    private boolean toLower = true;

    public CaseTransformFilter() {
    }

    public CaseTransformFilter(boolean toLowerCase) {
        this.toLower = toLowerCase;
    }

    public SeqFilterResult filterSequence(Sequence seq) {
        if (seq instanceof QSequence) {
            return new SeqFilterResult(new QSequence(seq.getSeqName(), seq.getDesc(), doFilter(seq.getSeqString().toLowerCase()), ((QSequence)seq).getQuality()));
        } else {
            return new SeqFilterResult(new Sequence(seq.getSeqName(), seq.getDesc(), doFilter(seq.getSeqString().toLowerCase())));
        }
    }

    private String doFilter(String seqString) {
        if (toLower) {
            return seqString.toLowerCase();
        } else {
            return seqString.toUpperCase();
        }
    }

    public String getName() {
        if (toLower) {
            return "Lower case transform";
        } else {
            return "Upper case transform";
        }
    }
}
