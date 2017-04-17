/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.msu.cme.rdp.seqfilter.output;

import edu.msu.cme.rdp.readseq.readers.Sequence;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author fishjord
 */
public class SeqFilterListOutput implements SeqFilterOutput {
    private List<Sequence> seqs = new ArrayList();

    public SeqFilterListOutput() {}

    public SeqFilterListOutput(List<Sequence> seqs) {
        this.seqs = seqs;
    }

    @Override
    public void appendSequence(Sequence s) {
        seqs.add(s);
    }

    public List<Sequence> getSeqs() {
        return seqs;
    }
}
