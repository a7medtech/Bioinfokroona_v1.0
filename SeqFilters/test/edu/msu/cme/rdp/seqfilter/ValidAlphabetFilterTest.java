/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.msu.cme.rdp.seqfilter;

import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import edu.msu.cme.rdp.readseq.QSequence;
import edu.msu.cme.rdp.seqfilter.filters.*;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.seqfilter.filters.ValidAlphabetFilter.ValidAlphabetFilterMode;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class ValidAlphabetFilterTest {

    private static final Sequence rnaTest = new QSequence("test_rna", "", "aAcCgGouUrRyYmMkKsxSwWbBdDhHvVnNq", new byte[]{1,   2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33});
    //a   A   c   C   g   G   o   u   U   r   R   y   Y   m   M   k   K   s   x   S   w   W   b   B   d   D   h   H   v   V   n   N   q
    //a   A   c   C   g   G       u   U   r   R   y   Y   m   M   k   K   s       S   w   W   b   B   d   D   h   H   v   V   n   N
    byte[] expectedRNAQual = {1,   2,  3,  4,  5,  6,      8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18,     20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
    private static final Sequence protTest = new Sequence("test_protein", "", "aAcCdDeEfFgGhHi2IkKlL;mMnNpPq4QrRsStTvVwWyYzZbBxX");

    private static final String expectedRNA = "aAcCgGuUrRyYmMkKsSwWbBdDhHvVnN";
    private static final String expectedProtein = "aAcCdDeEfFgGhHiIkKlLmMnNpPqQrRsStTvVwWyYzZbBxX";
    /**
     * Test of filterSequence method, of class ValidAlphabetFilter.
     */
    @Test
    public void testFilterSequence() {
        ValidAlphabetFilter filter = new ValidAlphabetFilter(ValidAlphabetFilterMode.DROP_BASE, SeqUtils.RNAAlphabet);

        SeqFilterResult result = filter.filterSequence(rnaTest);
        assertFalse(result.getErrorMessage(), result.failed());
        assertEquals("Expected QSeq", QSequence.class, result.getResultSeq().getClass());
        byte[] actualQual = ((QSequence)result.getResultSeq()).getQuality();
        assertEquals("Expected qual length != actual qual length", expectedRNAQual.length, actualQual.length);
        for(int index = 0;index < actualQual.length;index++) {
            assertEquals("Qual[" + index + "]", expectedRNAQual[index], actualQual[index]);
        }

        assertEquals(result.getResultSeq().getSeqString(), expectedRNA);

        filter = new ValidAlphabetFilter(ValidAlphabetFilterMode.DROP_SEQUENCE, SeqUtils.RNAAlphabet);

        result = filter.filterSequence(rnaTest);
        assertTrue("Sequence should have been dropped!", result.failed());

        filter = new ValidAlphabetFilter(ValidAlphabetFilterMode.DROP_BASE, SeqUtils.proteinAlphabet);

        result = filter.filterSequence(protTest);
        assertFalse(result.getErrorMessage(), result.failed());

        assertEquals(result.getResultSeq().getSeqString(), expectedProtein);

        filter = new ValidAlphabetFilter(ValidAlphabetFilterMode.DROP_SEQUENCE, SeqUtils.proteinAlphabet);

        result = filter.filterSequence(protTest);
        assertTrue("Sequence should have been dropped!", result.failed());
    }
}