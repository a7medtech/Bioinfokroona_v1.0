/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.seqfilter.filters;

import edu.msu.cme.rdp.readseq.QSequence;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.seqfilter.SeqFilter;
import edu.msu.cme.rdp.seqfilter.SeqFilterResult;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * @author fishjord
 */
public class ValidAlphabetFilter implements SeqFilter {

    public static enum ValidAlphabetFilterMode {

        DROP_SEQUENCE, DROP_BASE
    };
    private Set<Character> alphabet;
    private ValidAlphabetFilterMode mode;

    public ValidAlphabetFilter(ValidAlphabetFilterMode mode, Character... alphabet) {
        this(mode, new HashSet(Arrays.asList(alphabet)));
    }

    public ValidAlphabetFilter(ValidAlphabetFilterMode mode, Set<Character> alphabet) {
        this.alphabet = alphabet;
        this.mode = mode;
    }

    public SeqFilterResult filterSequence(Sequence seq) {
        if (seq instanceof QSequence) {
            return doFilter(seq.getSeqName(), seq.getDesc(), seq.getSeqString().toCharArray(), ((QSequence)seq).getQuality());
        } else {
            return doFilter(seq.getSeqName(), seq.getDesc(), seq.getSeqString().toCharArray(), null);
        }
    }

    private SeqFilterResult doFilter(String seqName, String desc, char[] seq, byte[] qual) {
        StringBuffer seqRet = new StringBuffer();
        List<Byte> qualRet = new ArrayList();
        Boolean drop = false;
        String invalidChars = "";

        for (int index = 0; index < seq.length; index++) {
            char c = seq[index];
            if (!alphabet.contains(c)) {
                if (mode == ValidAlphabetFilterMode.DROP_SEQUENCE) {
                    drop = true;
                    invalidChars += c;
                }
            } else {
                seqRet.append(c);
                if (qual != null && qual.length > index) {
                    qualRet.add(qual[index]);

                }
            }
        }

        if (drop) {
            return new SeqFilterResult("Sequence contains invalid characters \"" + invalidChars + "\"");
        }

        if (qual != null) {
            byte[] qualTmp = null;
            qualTmp = new byte[qualRet.size()];
            for (int index = 0; index < qualTmp.length; index++) {
                qualTmp[index] = qualRet.get(index);
            }
            return new SeqFilterResult(new QSequence(seqName, desc, seqRet.toString(), qualTmp));
        }

        return new SeqFilterResult(new Sequence(seqName, desc, seqRet.toString()));

    }

    public String getName() {
        return "Sequence Alphabet Filter";
    }
}
