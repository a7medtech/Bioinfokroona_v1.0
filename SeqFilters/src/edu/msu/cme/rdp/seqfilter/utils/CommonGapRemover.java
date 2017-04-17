package edu.msu.cme.rdp.seqfilter.utils;

import java.util.Set;

import edu.msu.cme.rdp.readseq.readers.Sequence;

public class CommonGapRemover {

    private Set<Integer> gapset = null;

    public String getName() {
        return "common gap remover";
    }

    public CommonGapRemover(Set<Integer> gapset) {
        this.gapset = gapset;
    }

    public String removeCommonGaps(Sequence seq) {
        return removeCommonGaps(seq.getSeqString());
    }

    public String removeCommonGaps(String seqstring) {
        StringBuffer ret = new StringBuffer();

        char[] bases = seqstring.toCharArray();
        for (int index = 0; index < bases.length; index++) {
            if (!gapset.contains(index)) {
                ret.append(bases[index]);
            }
        }
        return ret.toString();
    }
}
