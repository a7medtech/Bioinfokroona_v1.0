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

import edu.msu.cme.rdp.kmer.Kmer;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class NuclKmerGeneratorTest {

    @Test
    public void testGenerator() {
        String seq = "agtcgctacatgaactgactacttaggttaacgtcatgcctaagcttacatacg";
        NuclKmerGenerator kmer = new NuclKmerGenerator(seq, 30);
        long[] expected = new long[]{
            205411726394520252L,
            821646905578081008L,
            980744613098630081L,
            464213938573979398L,
            703934249689070619L,
            509893989542588525L,
            886654453563507124L,
            87853300433487571L,
            351413201733950286L,
            252731302328954169L,
            1010925209315816677L,
            584936323442725783L,
            33902284557209180L,
            135609138228836720L,
            542436552915346882L,
            1016824707054540553L,
            608534314397621287L,
            128294248376791199L,
            513176993507164796L,
            899786469421812209L,
            140381363866707908L,
            561525455466831635L,
            1093180317260479564L,
            913956755221377329L,
            197062507064968390L
        };

        String[] expectedKmers = new String[]{
            "agtcgctacatgaactgactacttaggtta",
            "gtcgctacatgaactgactacttaggttaa",
            "tcgctacatgaactgactacttaggttaac",
            "cgctacatgaactgactacttaggttaacg",
            "gctacatgaactgactacttaggttaacgt",
            "ctacatgaactgactacttaggttaacgtc",
            "tacatgaactgactacttaggttaacgtca",
            "acatgaactgactacttaggttaacgtcat",
            "catgaactgactacttaggttaacgtcatg",
            "atgaactgactacttaggttaacgtcatgc",
            "tgaactgactacttaggttaacgtcatgcc",
            "gaactgactacttaggttaacgtcatgcct",
            "aactgactacttaggttaacgtcatgccta",
            "actgactacttaggttaacgtcatgcctaa",
            "ctgactacttaggttaacgtcatgcctaag",
            "tgactacttaggttaacgtcatgcctaagc",
            "gactacttaggttaacgtcatgcctaagct",
            "actacttaggttaacgtcatgcctaagctt",
            "ctacttaggttaacgtcatgcctaagctta",
            "tacttaggttaacgtcatgcctaagcttac",
            "acttaggttaacgtcatgcctaagcttaca",
            "cttaggttaacgtcatgcctaagcttacat",
            "ttaggttaacgtcatgcctaagcttacata",
            "taggttaacgtcatgcctaagcttacatac",
            "aggttaacgtcatgcctaagcttacatacg"
        };

        int index = 0;

        while (kmer.hasNext()) {
            Kmer temp = kmer.next();
            long val = temp.getPart(0);
            assertEquals(expected[index], val);
            assertEquals(expectedKmers[index], temp.toString());
            assertEquals(index + 1, kmer.getPosition());
            index++;
        }
    }

    @Test
    public void testKmerwithN() {
        String seq = "agtcgctacatgaactgactacttaggttaacgNctacttaggttaacgtcatgcctaagcttacatacg";
        NuclKmerGenerator kmer = new NuclKmerGenerator(seq, 30);
        long[] expected = new long[]{
            205411726394520252L,
            821646905578081008L,
            980744613098630081L,
            464213938573979398L,
            
            513176993507164796L,
            899786469421812209L,
            140381363866707908L,
            561525455466831635L,
            1093180317260479564L,
            913956755221377329L,
            197062507064968390L
        };

        String[] expectedKmers = new String[]{
            "agtcgctacatgaactgactacttaggtta",
            "gtcgctacatgaactgactacttaggttaa",
            "tcgctacatgaactgactacttaggttaac",
            "cgctacatgaactgactacttaggttaacg",
            
            "ctacttaggttaacgtcatgcctaagctta",
            "tacttaggttaacgtcatgcctaagcttac",
            "acttaggttaacgtcatgcctaagcttaca",
            "cttaggttaacgtcatgcctaagcttacat",
            "ttaggttaacgtcatgcctaagcttacata",
            "taggttaacgtcatgcctaagcttacatac",
            "aggttaacgtcatgcctaagcttacatacg"
        };

        int index = 0;

        while (kmer.hasNext()) {
            Kmer temp = kmer.next();
            long val = temp.getPart(0);
            assertEquals(expected[index], val);
            assertEquals(expectedKmers[index], temp.toString());
            assertEquals(index + 1, kmer.getPosition());
            index++;
            if ( index == 4) break;
        }
        
        // after N
        while (kmer.hasNext()) {
            Kmer temp = kmer.next();
            long val = temp.getPart(0);
            assertEquals(expected[index], val);
            assertEquals(expectedKmers[index], temp.toString());
            assertEquals(index + 31, kmer.getPosition());
            index++;
        }
        
        
        String seq2 = "agtcgctacatgaactgactacttaggttaacgnct";
        kmer = new NuclKmerGenerator(seq2, 30);
        index = 0;

        while (kmer.hasNext()) {
            Kmer temp = kmer.next();
            long val = temp.getPart(0);
            assertEquals(expected[index], val);
            assertEquals(expectedKmers[index], temp.toString());
            assertEquals(index + 1, kmer.getPosition());
            index++;
        }
        assertEquals(index, 4);
    }
    
    
    @Test
    public void testInvalidKmer() {
        String str = "acgty";
        try {
            new NuclKmerGenerator(str, 32);
            fail("Kmer size of 32 should have triggered an exception");
        } catch(IllegalArgumentException e) {

        }

        try {
            new NuclKmerGenerator(str, 10);
            fail("Smaller sequence than kmer should have triggered an exception");
        } catch(IllegalArgumentException e) {

        }

        try {
            NuclKmerGenerator kmer = new NuclKmerGenerator(str, 2);
            assertEquals(1L, (long)kmer.next().getPart(0));
            assertEquals(1, kmer.getPosition());
            assertEquals(6L, (long)kmer.next().getPart(0));
            assertEquals(2, kmer.getPosition());
            assertEquals(11L, (long)kmer.next().getPart(0));
            assertEquals(3, kmer.getPosition());
            kmer.next();
            fail("Next should've triggered an exception for invalid character");
        } catch(IllegalArgumentException e) {

        }
    }
}
