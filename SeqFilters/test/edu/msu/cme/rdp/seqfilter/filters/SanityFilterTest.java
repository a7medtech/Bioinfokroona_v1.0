/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.msu.cme.rdp.seqfilter.filters;

import edu.msu.cme.rdp.readseq.QSequence;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.seqfilter.SeqFilterResult;
import edu.msu.cme.rdp.seqfilter.filters.SanityFilter.SanityFilterExpectedSeqs;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class SanityFilterTest {

    private Sequence goodSeq = new Sequence("good", "", "abcdefg");
    private QSequence goodQSeq = new QSequence("goodq", "", "abcdefg", new byte[]{1, 2, 3, 4, 5, 6, 7});
    private QSequence badQSeq = new QSequence("bad", "", "abcdefg", new byte[]{1, 2, 3, 4, 5, 6});

    public SanityFilterTest() {
    }

    @Test
    public void test() {
        SanityFilter filter = new SanityFilter(SanityFilterExpectedSeqs.Both);
        SeqFilterResult result;

        result = filter.filterSequence(goodSeq);
        assertFalse("Sequence failed when it shouldn't have", result.failed());

        result = filter.filterSequence(goodQSeq);
        assertFalse("Sequence failed when it shouldn't have", result.failed());

        result = filter.filterSequence(badQSeq);
        assertTrue("QSequence passed when it shouldn't have", result.failed());

        filter = new SanityFilter(SanityFilterExpectedSeqs.QSeq);

        result = filter.filterSequence(goodSeq);
        assertTrue("Sequence passed when it shouldn't have", result.failed());

        result = filter.filterSequence(goodQSeq);
        assertFalse("Sequence failed when it shouldn't have", result.failed());

        filter = new SanityFilter(SanityFilterExpectedSeqs.Seq);

        result = filter.filterSequence(goodSeq);
        assertFalse("Sequence failed when it shouldn't have", result.failed());

        result = filter.filterSequence(goodQSeq);
        assertTrue("Sequence passed when it shouldn't have", result.failed());
    }
}