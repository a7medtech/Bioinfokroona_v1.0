/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.msu.cme.rdp.seqfilter.output;

import edu.msu.cme.rdp.readseq.writers.FastaWriter;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import java.io.File;
import java.io.IOException;

/**
 *
 * @author fishjord
 */
public class SeqFilterFileOutput implements SeqFilterOutput {
    private FastaWriter out;

    public SeqFilterFileOutput(File outFile) throws IOException {
        out = new FastaWriter(outFile);
    }

    @Override
    public void appendSequence(Sequence s) {
        out.writeSeq(s);
    }
}
