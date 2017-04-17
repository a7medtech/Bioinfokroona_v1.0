/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.msu.cme.rdp.seqfilter;

import edu.msu.cme.rdp.seqfilter.filters.NSeqFilter;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class NSeqFilterTest {

    private static final Sequence veryBadSeq = new Sequence("Very_Bad_Sequence", "",  "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");
    private static final Sequence badSeq = new Sequence("Bad_Sequence", "",  "nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANN");
    private static final Sequence okSeq = new Sequence("Ok_Sequence", "",   "NAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAN");
    private static final Sequence goodSeq = new Sequence("Good_Sequence", "", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

    public NSeqFilterTest() {
    }

    /**
     * Test of filterSequence method, of class NSeqFilter.
     */
    @Test
    public void testFilterSequence() {
        NSeqFilter instance = new NSeqFilter(2);
        SeqFilterResult result;

        result = instance.filterSequence(veryBadSeq);
        assertTrue("Very bad sequence didn't fail like it should", result.failed());

        result = instance.filterSequence(badSeq);
        assertTrue("Bad sequence didn't fail like it should", result.failed());

        result = instance.filterSequence(okSeq);
        assertFalse("Ok sequence failed when it shouldn't have because " + result.getErrorMessage(), result.failed());

        result = instance.filterSequence(goodSeq);
        assertFalse("Good sequence failed when it shouldn't have because " + result.getErrorMessage(), result.failed());
    }

}