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
public class SeqLengthFilter implements SeqFilter {

    /**
     * Absolute - sequence length must be the specified length to pass filter
     * LessThan - sequence length must be less than specified length to pass filter
     * 
     */
    public enum LengthFilterType {

        Absolute, LessThan, GreaterThan, Range
    };
    private LengthFilterType filterType;
    private int length;
    private int slush;

    public SeqLengthFilter(int length, LengthFilterType type) {
        if (type == LengthFilterType.Range) {
            throw new IllegalArgumentException("To use the range filter please use the SeqLengthFilter(int length, int slush) constructor");
        }

        this.length = length;
        this.filterType = type;
    }

    public SeqLengthFilter(int length, int slush) {
        this.length = length;
        this.slush = slush;
        this.filterType = LengthFilterType.Range;
    }

    public String getName() {
        if(filterType == LengthFilterType.Range)
            return filterType.toString() + " " + length + "+/-" + slush + " sequence length filter";
        else
            return filterType.toString() + " " + length + " sequence length filter";
    }

    public SeqFilterResult filterSequence(Sequence seq) {
        int seqLength = seq.getSeqString().length();

        if (filterType == LengthFilterType.Absolute && seqLength != length) {
            return new SeqFilterResult("length [" + seqLength + "] doesn't match required length [" + length + "]");
        } else if (filterType == LengthFilterType.LessThan && seqLength > length) {
            return new SeqFilterResult("length [" + seqLength + "] is greater than maximum length [" + length + "]");
        } else if (filterType == LengthFilterType.GreaterThan && seqLength < length) {
            return new SeqFilterResult("length [" + seqLength + "] is less than minimum length [" + length + "]");
        } else if (filterType == LengthFilterType.Range && !(seqLength >= length - slush && seqLength <= length + slush)) {
            return new SeqFilterResult("length [" + seqLength + "] isn't in required range [" + (length - slush) + "-" + (length + slush));
        }

        return new SeqFilterResult(seq);
    }
}
