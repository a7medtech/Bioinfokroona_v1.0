/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.seqfilter.utils;

import edu.msu.cme.rdp.seqfilter.utils.TrimHelper.TrimResult;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class TrimHelperTest {

    private static final String[] primers1 = new String[]{"AACGTURYMKSWBDHVN"};
    private static final String[] primers2 = new String[]{"CMNBDNVNNNNNNNNNN"};
    private static final String[] primers3 = new String[]{"CMNBDNMNNNNNNNNNN"};
    private static final String testSeq = "ACGTURYMKSWBDHVNAAAAAAAAAAAAAAAAAAAACGTUCGT";
    private static final String expectedTrimmedSeq = "AAAAAAAAAAAAAAAAAAAACGTUCGT";
    private static final int expectedPrimer = 1;
    private static final int expectedStop = 16;
    private static final int expectedEdit = 16;
    
    private static final String expectedMisTrimmedSeq = "DHVNAAAAAAAAAAAAAAAAAAAACGTUCGT";
    private static final int expectedMisPrimer = 1;
    private static final int expectedMisStop = 12;
    private static final int expectedMisEdit = 12;

    @Test
    public void testPartialMatchExactPrimer() {
        TrimResult result = TrimHelper.partialMatch(testSeq, primers1);

        assertEquals("Incorrect trimmed seq", expectedTrimmedSeq, result.getTrimmedSeq());
        assertEquals("Incorrect primer distance", expectedPrimer, result.getPrimer());
        assertEquals("Incorrect primer stop distance", expectedStop, result.getPrimerStopIndex());
        assertEquals("Incorrect edit distance", expectedEdit, result.getEditDistance());

        result = TrimHelper.partialMatch(testSeq, primers2);

        assertEquals("Incorrect trimmed seq", expectedTrimmedSeq, result.getTrimmedSeq());
        assertEquals("Incorrect primer distance", expectedPrimer, result.getPrimer());
        assertEquals("Incorrect primer stop distance", expectedStop, result.getPrimerStopIndex());
        assertEquals("Incorrect edit distance", expectedEdit, result.getEditDistance());

        result = TrimHelper.partialMatch(testSeq, primers3);

        assertEquals("Incorrect trimmed seq", expectedMisTrimmedSeq, result.getTrimmedSeq());
        assertEquals("Incorrect primer distance", expectedMisPrimer, result.getPrimer());
        assertEquals("Incorrect primer stop distance", expectedMisStop, result.getPrimerStopIndex());
        assertEquals("Incorrect edit distance", expectedMisEdit, result.getEditDistance());
    }
}
