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
public class ExpQualityFilter implements SeqFilter {

    private static final double[] QUALITY_MAPPING = new double[100];

    static {
        for (int q = 0; q < QUALITY_MAPPING.length; q++) {
            QUALITY_MAPPING[q] = Math.pow(10, -((double)q / 10));
        }
    }
    private int minQuality;

    public ExpQualityFilter(int minQuality) {
        this.minQuality = minQuality;
    }

    public SeqFilterResult filterSequence(Sequence inseq) {
        if (!(inseq instanceof QSequence)) {
            return new SeqFilterResult(inseq);
        }
        QSequence seq = (QSequence) inseq;

        double avgEQual = 0;
        for (byte q : seq.getQuality()) {
            avgEQual += QUALITY_MAPPING[q];
        }

        double avgQual = -10 * Math.log10(avgEQual / seq.getQuality().length);

        if (avgQual < minQuality) {
            return new SeqFilterResult("Quality score=" + avgQual);
        }

        return new SeqFilterResult(seq);
    }

    public String getName() {
        return "Exponential Quality Filter";
    }
}
