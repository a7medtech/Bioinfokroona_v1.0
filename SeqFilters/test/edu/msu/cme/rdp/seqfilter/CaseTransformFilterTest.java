/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.msu.cme.rdp.seqfilter;

import edu.msu.cme.rdp.seqfilter.filters.*;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.seqfilter.SeqFilterResult;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class CaseTransformFilterTest {
    private static final Sequence testSeq = new Sequence("test", "", "ABCDEFabcdef");

    private static final String lowerExpected = "abcdefabcdef";
    private static final String upperExpected = "ABCDEFABCDEF";

    /**
     * Test of filterSequence method, of class CaseTransformFilter.
     */
    @Test
    public void testFilterSequence() {
        CaseTransformFilter filter = new CaseTransformFilter();

        SeqFilterResult result = filter.filterSequence(testSeq);
        assertFalse("Test failed when it shouldn't have...ever", result.failed());

        assertEquals(result.getResultSeq().getSeqString(), lowerExpected);

        filter = new CaseTransformFilter(false);

        result = filter.filterSequence(testSeq);
        assertFalse("Test failed when it shouldn't have...ever", result.failed());

        assertEquals(result.getResultSeq().getSeqString(), upperExpected);
    }

}