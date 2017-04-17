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
public class CharReplaceFilterTest {

    public CharReplaceFilterTest() {
    }
    
    private static final Sequence testSeq = new Sequence("test1", "test", "aAcCgGtTuUxXnN");
    private static final String expectedString = "aAcCgGuUuUnNnN";

    /**
     * Test of filterSequence method, of class CharReplaceFilter.
     */
    @Test
    public void testFilterSequence() {
        CharReplaceFilter filter = new CharReplaceFilter(CharReplaceFilter.rnaReplaceMap);

        SeqFilterResult result = filter.filterSequence(testSeq);
        assertFalse("Test failed when it shouldn't have...ever", result.failed());

        assertEquals(result.getResultSeq().getSeqString(), expectedString);
    }

}