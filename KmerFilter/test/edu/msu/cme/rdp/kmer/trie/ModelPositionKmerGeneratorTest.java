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
package edu.msu.cme.rdp.kmer.trie;

import edu.msu.cme.rdp.readseq.SequenceType;
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
            new char[]{'B', 'C', 'D'},
            new char[]{'F', 'G', 'H'});
    private static final int[] expectedModelPos1 = new int[]{2, 5};
    private static final String testSeq2 = ".ABCdE.F.G.HijkLMXOPQ";
    private static final List<char[]> expected2 = Arrays.asList(
            new char[]{'A', 'B', 'C'},
            new char[]{'E', 'F', 'G'},
            new char[]{'F', 'G', 'H'},
            new char[]{'O', 'P', 'Q'});
    private static final int[] expectedModelPos2 = new int[]{1, 4, 5, 11};

    @Test
    public void testGetKmers() {
        List<char[]> result = ModelPositionKmerGenerator.getKmers(testSeq1, k, SequenceType.Protein);

        for (int index = 0; index < result.size(); index++) {
            assertArrayEquals("Testing index " + index + " of expected " + new String(expected1.get(index)) + " compared to " + new String(result.get(index)), expected1.get(index), result.get(index));
        }

        assertEquals(expected1.size(), result.size());

        result = ModelPositionKmerGenerator.getKmers(testSeq2, k, SequenceType.Protein);

        for (int index = 0; index < result.size(); index++) {
            assertArrayEquals("Testing index " + index + " of expected " + new String(expected2.get(index)) + " compared to " + new String(result.get(index)), expected2.get(index), result.get(index));
        }

        assertEquals(expected2.size(), result.size());
    }

    @Test
    public void testIterator() {
        ModelPositionKmerGenerator instance = new ModelPositionKmerGenerator(testSeq1, k, SequenceType.Protein);
        int index = 0;

        while (instance.hasNext()) {
            char[] kmer = instance.next();
            assertArrayEquals("Testing index " + index + " of expected " + new String(expected1.get(index)) + " compared to " + new String(kmer), expected1.get(index), kmer);
            assertEquals("Model position check", expectedModelPos1[index], instance.getModelPosition());
            index++;
        }

        assertEquals(expected1.size(), index);

        instance = new ModelPositionKmerGenerator(testSeq2, k, SequenceType.Protein);
        index = 0;

        while (instance.hasNext()) {
            char[] kmer = instance.next();
            assertArrayEquals("Testing index " + index + " of expected " + new String(expected2.get(index)) + " compared to " + new String(kmer), expected2.get(index), kmer);
            assertEquals("Model position check", expectedModelPos2[index], instance.getModelPosition());
            index++;
        }

        assertEquals(expected2.size(), index);
    }

    @Test
    public void testIterable() {
        ModelPositionKmerGenerator instance = new ModelPositionKmerGenerator(testSeq1, k, SequenceType.Protein);
        int index = 0;

        for (char[] kmer : instance) {
            assertArrayEquals("Testing index " + index + " of expected " + new String(expected1.get(index)) + " compared to " + new String(kmer), expected1.get(index++), kmer);
        }

        assertEquals(expected1.size(), index);

        instance = new ModelPositionKmerGenerator(testSeq2, k, SequenceType.Protein);
        index = 0;

        for (char[] kmer : instance) {
            assertArrayEquals("Testing index " + index + " of expected " + new String(expected2.get(index)) + " compared to " + new String(kmer), expected2.get(index++), kmer);
        }

        assertEquals(expected2.size(), index);
    }
    
    
    /*
    @Test
    public void testReal() {
        String seq = "........................................................---MRKVAIYGKG.GIGKSTTTQNTVAGLA.E.M.G.R.K.IMVVGCDPKADSTRLLLGGLAQKSVLDT..LRE.E.....G..........-....-...E..D......V.......E.LDDIRKPG...FGG............TWCVESGGPEPGVGCAGRGIITSINMLES.L.G..A.Y...EeseG....L..DYAFY.DVLGDVVCGGFAMPIRDGKAEEIYIVCSGEMMAMYAANNICKGIMKYAE..S..GGVR..LGGLICNSRNT.D.READLITELANKLGTQMIYFVPRDNDVQRAEINRKTVIE.WDGTVPQ.ADQYRGLAKAIDGN.T.-MFTVPTPL...E.IEDLEQLLLDYGIME----..-..--..-----------a........................................................";
        String modelSeq = "---MRKVAIYGKGIGKSTTTQNTVAGLAMVVGCDPKADSTRLLLGGLAQKSVLDTRE--DDIRKPGGGWCVESGGPEPGVGCAGRGIITSINMLESGYAFYVLGDVVCGGFAMPIRDGKAEEIYIVCSGEMMAMYAANNICKGIMKYAEGVRGGLICNSRNTEADLITELANKLGTQMIYFVPRDNDVQRAEINRKTVIEDGTVPQDQYRGLAKAIDGN-MFTVPTPLEDLEQLLLDYGIME------------------";
        System.out.println(modelSeq);
        
        ModelPositionKmerGenerator instance = new ModelPositionKmerGenerator(seq, 10);
        int index = 0;

        while (instance.hasNext()) {
            String next = new String(instance.next());
            System.out.println(instance.getModelPosition() + " " + next);
        }
        
        seq = "........................................................---MRQIAIYGKG.GIGKSTTTQNLTASLS.T.M.G.N.K.IMLVGCDPKADSTRMLLGGLNQKTVLDT..LRS.E.....G..........D....-...E..G......V.......D.LDVVMQRG...FGD............IKCVESGGPEPGVGCAGRGIITSIGLLEN.L.G..A.Y...T..dD....L..DYVFY.DVLGDVVCGGFAMPIREGKAKEIYIVASGELMAIYAANNICKGLAKFAK..-..GGAR..LGGIICNSRNV.D.GERELLDAFAKKLGSQLIHFIPRDNIVQRAEINRKTVID.FDPESNQ.AKEYLTLAHNVQNN.D.-KLVVPTPL...P.MEELESMMVEFGIVD----..-..--..-----------m........................................................";
        instance = new ModelPositionKmerGenerator(seq, 10);
        System.out.println();

        while (instance.hasNext()) {
            String next = new String(instance.next());
            System.out.println(instance.getModelPosition() + " " + next);
        }
    }*/
}
