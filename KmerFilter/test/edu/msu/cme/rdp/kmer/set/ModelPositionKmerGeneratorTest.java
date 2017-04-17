/*
 * Copyright (C) 2012 Jordan Fish <fishjord at msu.edu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package edu.msu.cme.rdp.kmer.set;

import edu.msu.cme.rdp.kmer.trie.*;
import edu.msu.cme.rdp.readseq.SequenceType;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class ModelPositionKmerGeneratorTest {

    private static final String testSeq1 = "-BCDeFG.Hijk";
    private static final int k = 3;
    private static final List<char[]> expected1 = Arrays.asList(
            new char[]{'b', 'c', 'd'},
            new char[]{'f', 'g', 'h'});
    private static final int[] expectedModelPos1 = new int[]{2, 5};
    private static final String testSeq2 = ".ABCdE.F.G.HijkLMXAPQ";
    private static final List<char[]> expected2 = Arrays.asList(
            new char[]{'a', 'b', 'c'},
            new char[]{'e', 'f', 'g'},
            new char[]{'f', 'g', 'h'},
            new char[]{'a', 'p', 'q'});
    private static final int[] expectedModelPos2 = new int[]{1, 4, 5, 11};
    private static final String testSeq3 = ".AGCtT.C.G.TijkTCCATG";
    private static final List<char[]> expected3 = Arrays.asList(
            new char[]{'a', 'g', 'c'},
            new char[]{'t', 'c', 'g'},
            new char[]{'c', 'g', 't'},
            new char[]{'t', 'c', 'c'},
            new char[]{'c', 'c', 'a'},
            new char[]{'c', 'a', 't'},
            new char[]{'a', 't', 'g'});
    private static final int[] expectedModelPos3 = new int[]{1, 4, 5, 8, 9, 10, 11, 12, 13, 14};

    @Test
    public void testGetKmers() {
        List<char[]> result = new ArrayList();
        ProtKmerGenerator kmer = new ProtKmerGenerator(testSeq1, k, true);
        while(kmer.hasNext()) {
            result.add(kmer.next().toString().toCharArray());
        }

        for (int index = 0; index < result.size(); index++) {
            assertArrayEquals("Testing index " + index + " of expected " + new String(expected1.get(index)) + " compared to " + new String(result.get(index)), expected1.get(index), result.get(index));
        }

        assertEquals(expected1.size(), result.size());

        kmer = new ProtKmerGenerator(testSeq2, k, true);
        result.clear();
        while(kmer.hasNext()) {
            result.add(kmer.next().toString().toCharArray());
        }

        for (int index = 0; index < result.size(); index++) {
            assertArrayEquals("Testing index " + index + " of expected " + new String(expected2.get(index)) + " compared to " + new String(result.get(index)), expected2.get(index), result.get(index));
        }

        assertEquals(expected2.size(), result.size());
    }

    @Test
    public void testIterator() {
        ProtKmerGenerator kmerIt = new ProtKmerGenerator(testSeq1, k, true);
        int index = 0;

        while (kmerIt.hasNext()) {
            char[] kmer = kmerIt.next().toString().toCharArray();
            assertArrayEquals(new String(kmer) + " testing index " + index + " of expected " + new String(expected1.get(index)) + " compared to " + new String(kmer), expected1.get(index), kmer);
            assertEquals("Model position check", expectedModelPos1[index], kmerIt.getPosition());
            index++;
        }

        assertEquals(expected1.size(), index);

        kmerIt = new ProtKmerGenerator(testSeq2, k, true);
        index = 0;

        while (kmerIt.hasNext()) {
            char[] kmer = kmerIt.next().toString().toCharArray();
            assertArrayEquals("Testing index " + index + " of expected " + new String(expected2.get(index)) + " compared to " + new String(kmer), expected2.get(index), kmer);
            assertEquals("Model position check", expectedModelPos2[index], kmerIt.getPosition());
            index++;
        }

        assertEquals(expected2.size(), index);
    }

    @Test
    public void testNuclIterator() {
        ProtKmerGenerator kmerIt = new ProtKmerGenerator(testSeq3, k, true);
        int index = 0;

        while (kmerIt.hasNext()) {
            char[] kmer = kmerIt.next().toString().toCharArray();
            assertArrayEquals(new String(kmer) + " testing index " + index + " of expected " + new String(expected3.get(index)) + " compared to " + new String(kmer), expected3.get(index), kmer);
            assertEquals("Model position check", expectedModelPos3[index], kmerIt.getPosition());
            index++;
        }

        assertEquals(expected3.size(), index);
    }
}
