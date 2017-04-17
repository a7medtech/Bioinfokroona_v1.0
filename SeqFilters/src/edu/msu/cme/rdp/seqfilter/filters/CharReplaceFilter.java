/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.msu.cme.rdp.seqfilter.filters;

import edu.msu.cme.rdp.readseq.QSequence;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.seqfilter.SeqFilter;
import edu.msu.cme.rdp.seqfilter.SeqFilterResult;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author fishjord
 */
public class CharReplaceFilter implements SeqFilter {
    private Map<Character, Character> replaceMap;

    public static final Map<Character, Character> rnaReplaceMap;
    public static final Map<Character, Character> dnaReplaceMap;
    static {
        Map<Character, Character> tmpMap = new HashMap();
        tmpMap.put('X', 'N');
        tmpMap.put('x', 'n');
        tmpMap.put('T', 'U');
        tmpMap.put('t', 'u');
        tmpMap.put('I', 'D');
        tmpMap.put('i', 'd');
        rnaReplaceMap = Collections.unmodifiableMap(tmpMap);

        tmpMap = new HashMap();
        tmpMap.put('X', 'N');
        tmpMap.put('x', 'n');
        tmpMap.put('U', 'T');
        tmpMap.put('u', 't');
        tmpMap.put('I', 'D');
        tmpMap.put('i', 'd');
        dnaReplaceMap = Collections.unmodifiableMap(tmpMap);
    }

    public CharReplaceFilter(Map<Character, Character> replaceMap) {
        this.replaceMap = replaceMap;
    }

    public SeqFilterResult filterSequence(Sequence seq) {
        if(seq instanceof QSequence) {
            return new SeqFilterResult(new QSequence(seq.getSeqName(), seq.getDesc(), doFilter(seq.getSeqString()), ((QSequence)seq).getQuality()));
        } else {
            return new SeqFilterResult(new Sequence(seq.getSeqName(), seq.getDesc(), doFilter(seq.getSeqString())));
        }
    }

    private String doFilter(String seqString) {
        for(Character k : replaceMap.keySet()) {
            seqString = seqString.replace(k, replaceMap.get(k));
        }

        return seqString;
    }

    public String getName() {
        return "Character Replacement Filter";
    }

}
