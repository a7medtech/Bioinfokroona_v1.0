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
public class SeqFilterResult {
    private Sequence resultSeq;
    private String errorMessage = "";

    public SeqFilterResult(Sequence resultSeq) {
        this.resultSeq = resultSeq;
    }

    public SeqFilterResult(String errorMessage) {
        this.errorMessage = errorMessage;
    }

    public boolean failed() {
        return resultSeq == null;
    }

    public Sequence getResultSeq() {
        return resultSeq;
    }

    public String getErrorMessage() {
        return errorMessage;
    }
}
