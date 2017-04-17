/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.msu.cme.rdp.seqfilter;

import edu.msu.cme.rdp.readseq.QSequence;
import edu.msu.cme.rdp.seqfilter.filters.PrimerFilter;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import java.util.Arrays;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class PrimerFilterTest {

    private static final String[] fPrimer = {"AYTGGGYDTAAAGNG", "TYTGGGYDTAAAGNG"};

    private static final String rPrimer[] = {"TACCRGGGTHTCTAATCC", "TACCAGAGTATCTAATTC", "CTACDSRGGTMTCTAATC", "TACNVGGGTATCTAATCC"};
    private static final String rPrimer_mod[] = {"TACCRGGGTHTCTAATCC", "TACCAGAGTATCTAATTC", "TACDSRGGTMTCTAATCN", "TACNVGGGTATCTAATCC"};
    //1 GGATTAGADACCCYGGTA  2 GAATTAGATACTCTGGTA  3 GATTAGAKACCYSHGTAG  4 GGATTAGATACCCBNGTA
    //1 GGATTAGA[A,T]ACCC[C,T]GGTA  2 GAATTAGATACTCTGGTA  3 GATTAGA[T,G]ACC[C,T][C,G][A,T,C]GTAG  4 GGATTAGATACCC[C,T,G][A,C,T,G]GTA
    private static final int maxEditDist = 2;

    //match fprimer, match exact rPrimer[0]
    private static final QSequence goodSeq1 = new QSequence(
            "EQKL81L01D8Y7X",
            "length=244",
            "ATCTACTGGGTGTAAAGGGTGCGTAGGCGGGTTTTTAAGTCAGGGGTGAAATCCTGGAGCTCAACTCCAGAACTGCCTTTGATACTGAAGATCTTGAGTCCGGGAGAGGTGAGTGGAACTGCGAGTGTAGAGGTGAAATTCGTAGATATTCGCAAGAACACCAGTGGCGAAGGCGGCTCACTGGCCCGGTACTGACGCCGAGGCACGAAAGCGTGGGGAGCAAACAGGATTAGATACCCCGGTA",
            new byte[] {25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 27, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 30, 25, 25, 25, 25, 25, 25, 25, 30, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25}
    );
    //match fprimer, 1 mismatch to rPrimer[1]
    private static final Sequence goodSeq2 = new Sequence("EQKL81L01DVMLI", "length=244", "ATCTATTGGGTGTAAAGCGCTCGTAGGCGGCTTGTCGCGTCGTCCGTGAAAACCTGGGGCTCAACTCCAGGCTTGCGGTCGATACGGGCAGGCTCGAGTTCGGTAGGGGAGACTGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGGTGGCGAAGGCGGGCCTCTGGGCCGATACTGACGCTGAGGAGCGAAAGCGTGGGGAGCGAACAGATTAGAGACCYAHGTAG");
    //match fprimer, 1 mismatch to rPrimer[3], has 1 n
    //private static final Sequence badSeq1 = new Sequence("EQKL81L01AS83V", "length=190", "ATCTACTGGGTGTAAAGGGCTCGTAGGCGGCCGACTAAGTCAGTTGTGAAATCCCTCGGCTTAACCGGGGAATTGCGTCTGATACTGGACGGCTTGAGTTTGGGAGAGGGATGCGGAATTCCAGGTGTAGCGGTGAAATGCGTAGATATCTGGAGGAACACCGGTGGCGAAGGCGGCANTCTCGACACACGGATTAGATACCCBNGAA");
    //match fprimer, doesn't match rPrimer
    private static final Sequence badSeq2 = new Sequence("EQKL81L01DEOFU", "length=194", "ATCTACTGGGTATAAAGGGTGCGTAGGTGGTTGTTTAAGTCTGTTGTGAAAGCCCTGGGCTCAACCTGGGAACTGCAGTGGAAACTGGACGACTAGAGTGTGGTGGAGGGTAGCGGAATTCCTGGTGTAGCAGTGAAATGCGTAGAGATCAGGAGGAACATCCATGGCGAAGGCAGCTACCTGGACCAACACTA");
    //match fprimer, match partial primer rPrimer[3]
    private static final Sequence partial = new Sequence("partialend", "length=236", "TCATATTGGGTGTAAAGGGAGCGCAGGCGGAATGATAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAACTGCATCGGAAACTGTCATTCTTGAGTGCAGAAGAGGAGAGTGGAATTCCATGTGTAGCGGTGGAATGCGTAGATATATGGAAGAACACCAGTGGCGAAGGCGGCTCTCTGGTCTGCAACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGATA");
    //reverse primer is too short
    private static final Sequence shortSeq = new Sequence("shortend", "length=236", "TCATATTGGGTGTAAAGGGAGCGCAGGCGGAATGATAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAACTGCATCGGAAACTGTCATTCTTGAGTGCAGAAGAGGAGAGTGGAATTCCATGTGTAGCGGTGGAATGCGTAGATATATGGAAGAACACCAGTGGCGAAGGCGGCTCTCTGGTCTGCAACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATT");
    // match modfPrimer[1]
    private static final Sequence matchf2_goodSeq1 = new Sequence("matchfprimer2", "length=244", "ATCTTCTGGGTGTAAAGGGTGCGTAGGCGGGTTTTTAAGTCAGGGGTGAAATCCTGGAGCTCAACTCCAGAACTGCCTTTGATACTGAAGATCTTGAGTCCGGGAGAGGTGAGTGGAACTGCGAGTGTAGAGGTGAAATTCGTAGATATTCGCAAGAACACCAGTGGCGAAGGCGGCTCACTGGCCCGGTACTGACGCCGAGGCACGAAAGCGTGGGGAGCAAACAGGATTAGATACCCCGGTA");

    private static final Sequence partialAmbigSeq = new Sequence("F3SGFHA01C0UH2", "", "AGAGAGATTGGGTGTAAAAGGGCGCGTAGGCGGTTTGTTAAGTGTGAAGTGAAATGCCTGGGCTCAACCTGGGACGTGCTTTGCATACTGATGAACTTGAGTCCAAGAGGGGGTGGTGGAATTCCTGGTGTAGGGGTGAAATCCGTAGATATCAGGAGGAACACCGTTGGCGAAGGCGGCCACCTGATTGGTACTGACGCTGAGGCGCGAAAGCGTGGGGAGCGAGCAGGATTAGATACCCGTG");

    private static final Sequence parital_badSeq = new Sequence("parital_badSeq", "", "CTCTGATCACTGGGACGTAAAAGCGAGTGCAGGTCGGACTACGAATAAGTCTAGATGTAGAAACGCCTTCGCTCAACCGGAGAATTGACAATACAAGAAACTGTCGAGCTTGAGTACAAGAAGAGAGAGTGGAACTCCATGTGTAGCGGTAGAAATGCGTAAGATATAGTGGAAACGAACCACCGGTGACGAAGCGCT");
    private static final String good1Expected = "TGCGTAGGCGGGTTTTTAAGTCAGGGGTGAAATCCTGGAGCTCAACTCCAGAACTGCCTTTGATACTGAAGATCTTGAGTCCGGGAGAGGTGAGTGGAACTGCGAGTGTAGAGGTGAAATTCGTAGATATTCGCAAGAACACCAGTGGCGAAGGCGGCTCACTGGCCCGGTACTGACGCCGAGGCACGAAAGCGTGGGGAGCAAACA";
    private static final byte[] good1QualExpected = new byte[] {27, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 30};
    private static final String good2Expected = "CTCGTAGGCGGCTTGTCGCGTCGTCCGTGAAAACCTGGGGCTCAACTCCAGGCTTGCGGTCGATACGGGCAGGCTCGAGTTCGGTAGGGGAGACTGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGGTGGCGAAGGCGGGCCTCTGGGCCGATACTGACGCTGAGGAGCGAAAGCGTGGGGAGCGAACA";
    private static final String partialExpected = "AGCGCAGGCGGAATGATAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAACTGCATCGGAAACTGTCATTCTTGAGTGCAGAAGAGGAGAGTGGAATTCCATGTGTAGCGGTGGAATGCGTAGATATATGGAAGAACACCAGTGGCGAAGGCGGCTCTCTGGTCTGCAACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACA";
    private static final String secondForwardExpected = "TGCGTAGGCGGGTTTTTAAGTCAGGGGTGAAATCCTGGAGCTCAACTCCAGAACTGCCTTTGATACTGAAGATCTTGAGTCCGGGAGAGGTGAGTGGAACTGCGAGTGTAGAGGTGAAATTCGTAGATATTCGCAAGAACACCAGTGGCGAAGGCGGCTCACTGGCCCGGTACTGACGCCGAGGCACGAAAGCGTGGGGAGCAAACA";
    private static final String partialAmbigExpected = "CGCGTAGGCGGTTTGTTAAGTGTGAAGTGAAATGCCTGGGCTCAACCTGGGACGTGCTTTGCATACTGATGAACTTGAGTCCAAGAGGGGGTGGTGGAATTCCTGGTGTAGGGGTGAAATCCGTAGATATCAGGAGGAACACCGTTGGCGAAGGCGGCCACCTGATTGGTACTGACGCTGAGGCGCGAAAGCGTGGGGAGCGAGCA";

    /**
     * Test of filterSequence method, of class PrimerFilter.
     */
    @Test
    public void testFilterSequence() {
        PrimerFilter filter = new PrimerFilter(Arrays.asList(fPrimer), Arrays.asList(rPrimer), maxEditDist, maxEditDist, true, false);
        SeqFilterResult result;
        Sequence testSeq;

        testSeq = goodSeq1;
        result = filter.filterSequence(testSeq);
        assertFalse(testSeq.getSeqName() + " failed primer trimming...it shouldn't have! " + result.getErrorMessage(), result.failed());
        assertEquals(testSeq.getSeqName() + " didn't produce the expected trimming result", result.getResultSeq().getSeqString(), good1Expected);
        assertEquals("Expected to get a qseq back from goodSeq1", QSequence.class, result.getResultSeq().getClass());

        byte[] resultQual = ((QSequence)result.getResultSeq()).getQuality();
        assertEquals("Expected qual seq length doesn't match actual qual seq length",good1QualExpected.length, resultQual.length);
        for(int index = 0;index < good1QualExpected.length;index++) {
            assertEquals(good1QualExpected[index], resultQual[index]);
        }

        testSeq = goodSeq2;
        result = filter.filterSequence(testSeq);
        assertFalse(testSeq.getSeqName() + " failed primer trimming...it shouldn't have! " + result.getErrorMessage(), result.failed());
        assertEquals(testSeq.getSeqName() + " didn't produce the expected trimming result", result.getResultSeq().getSeqString(), good2Expected);

        testSeq = partial;
        result = filter.filterSequence(testSeq);
        assertFalse(testSeq.getSeqName() + " failed primer trimming...it shouldn't have! " + result.getErrorMessage(), result.failed());
        assertEquals(testSeq.getSeqName() + " didn't produce the expected trimming result", result.getResultSeq().getSeqString(), partialExpected);

        testSeq = partialAmbigSeq;
        result = filter.filterSequence(testSeq);
        assertFalse(testSeq.getSeqName() + " failed primer trimming...it shouldn't have! " + result.getErrorMessage(), result.failed());
        assertEquals(testSeq.getSeqName() + " didn't produce the expected trimming result", result.getResultSeq().getSeqString(), partialAmbigExpected);

        testSeq = matchf2_goodSeq1;
        result = filter.filterSequence(testSeq);
        assertFalse(testSeq.getSeqName() + " failed primer trimming...it shouldn't have! " + result.getErrorMessage(), result.failed());
        assertEquals(testSeq.getSeqName() + " didn't produce the expected trimming result", result.getResultSeq().getSeqString(), secondForwardExpected);

        /*testSeq = badSeq1;
        result = filter.filterSequence(testSeq);
        assertTrue(testSeq.getSeqName() + " passed when it shouldn't have", result.failed());*/

        testSeq = badSeq2;
        result = filter.filterSequence(testSeq);
        assertTrue(testSeq.getSeqName() + " passed when it shouldn't have", result.failed());

        testSeq = shortSeq;
        result = filter.filterSequence(testSeq);
        assertTrue(testSeq.getSeqName() + " passed when it shouldn't have", result.failed());

        testSeq = parital_badSeq;
        result = filter.filterSequence(testSeq);
        assertTrue(testSeq.getSeqName() + " passed when it shouldn't have", result.failed());

        // change the filter
        int maxRevDist = 0;
        filter = new PrimerFilter(Arrays.asList(fPrimer), Arrays.asList(rPrimer_mod), maxEditDist, maxRevDist, true, true);

        testSeq = parital_badSeq;
        result = filter.filterSequence(testSeq);
        assertTrue(testSeq.getSeqName() + " passed when it shouldn't have", result.failed());

    }
}