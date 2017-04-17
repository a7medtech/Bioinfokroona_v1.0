/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.msu.cme.rdp.seqfilter.filters;

import edu.msu.cme.rdp.readseq.QSequence;
import edu.msu.cme.rdp.seqfilter.SeqFilterResult;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class ExpQualityFilterTest {

    private static final QSequence goodSeq = new QSequence("good", "", "abcdefg", new byte[] {25, 26, 25, 20, 30, 25, 24});
    private static final QSequence badSeq = new QSequence("bad", "", "abcdefg", new byte[] {20, 10, 10, 5, 30, 5, 10});

    public ExpQualityFilterTest() {
    }

    @Test
    public void test() {
        ExpQualityFilter filter = new ExpQualityFilter(20);
        SeqFilterResult result;

        result = filter.filterSequence(goodSeq);
        assertFalse("Good sequence failed <pout>", result.failed());
        result = filter.filterSequence(badSeq);
        assertTrue("Good heavens...a bad sequence passed", result.failed());
    }
}