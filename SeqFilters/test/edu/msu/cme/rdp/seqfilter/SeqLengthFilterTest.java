/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.msu.cme.rdp.seqfilter;

import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.seqfilter.filters.SeqLengthFilter;
import edu.msu.cme.rdp.seqfilter.filters.SeqLengthFilter.LengthFilterType;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class SeqLengthFilterTest {

    private static final Sequence absBadSeq = new Sequence("Absolute_Length_Bad_Seq", "",  "AAAA");
    private static final Sequence absGoodSeq = new Sequence("Absolute_Length_Good_Seq", "",  "AAAAA");
    private static final Sequence lessThanBadSeq = new Sequence("Less_Than_Bad_Seq", "",  "AAAAAA");
    private static final Sequence lessThanGoodSeq = new Sequence("Less_Than_Good_Seq", "",  "AAAA");
    private static final Sequence greaterThanBadSeq = new Sequence("Greater_Than_Bad_Seq", "",  "AAAA");
    private static final Sequence greaterThanGoodSeq = new Sequence("Greater_Than_Good_Seq", "",  "AAAAAA");
    private static final Sequence rangeLowerSeq = new Sequence("Range_Low_Seq", "",   "AA");
    private static final Sequence rangeHigherSeq = new Sequence("Range_High_Seq", "", "AAAAAAAAA");
    private static final Sequence rangeGoodSeq = new Sequence("Range_Good_Seq", "", "AAAAAA");

    public SeqLengthFilterTest() {
    }

    /**
     * Test of filterSequence method, of class SeqLengthFilter.
     */
    @Test
    public void testAbsLength() {
        SeqLengthFilter filter = new SeqLengthFilter(5, LengthFilterType.Absolute);
        SeqFilterResult result;

        result = filter.filterSequence(absBadSeq);
        assertTrue(absBadSeq.getSeqName() + " passed when it shouldn't have", result.failed());

        result = filter.filterSequence(absGoodSeq);
        assertFalse(absGoodSeq.getSeqName() + " failed when it shouldn't have " + result.getErrorMessage(), result.failed());
    }

    @Test
    public void testLessThanLength() {
        SeqLengthFilter filter = new SeqLengthFilter(5, LengthFilterType.LessThan);
        SeqFilterResult result;

        result = filter.filterSequence(lessThanBadSeq);
        assertTrue(lessThanBadSeq.getSeqName() + " passed when it shouldn't have", result.failed());

        result = filter.filterSequence(lessThanGoodSeq);
        assertFalse(lessThanGoodSeq.getSeqName() + " failed when it shouldn't have " + result.getErrorMessage(), result.failed());
    }

    @Test
    public void testGreaterThanLength() {
        SeqLengthFilter filter = new SeqLengthFilter(5, LengthFilterType.GreaterThan);
        SeqFilterResult result;

        result = filter.filterSequence(greaterThanBadSeq);
        assertTrue(greaterThanBadSeq.getSeqName() + " passed when it shouldn't have", result.failed());

        result = filter.filterSequence(greaterThanGoodSeq);
        assertFalse(greaterThanGoodSeq.getSeqName() + " failed when it shouldn't have " + result.getErrorMessage(), result.failed());
    }

    @Test
    public void testRangeLength() {
        SeqLengthFilter filter = new SeqLengthFilter(5, 2);
        SeqFilterResult result;

        result = filter.filterSequence(rangeLowerSeq);
        assertTrue(rangeLowerSeq.getSeqName() + " passed when it shouldn't have", result.failed());

        result = filter.filterSequence(rangeHigherSeq);
        assertTrue(rangeHigherSeq.getSeqName() + " passed when it shouldn't have", result.failed());

        result = filter.filterSequence(rangeGoodSeq);
        assertFalse(rangeGoodSeq.getSeqName() + " failed when it shouldn't have " + result.getErrorMessage(), result.failed());
    }

}