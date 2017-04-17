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
public class KmerIteratorTest {

    @Test
    public void testGenerator() {
        String seq = "agtcgctacatgaactgactacttaggttaacgtcatgcctaagcttacatacgcgcgccgg";
        KmerIterator kmer = new KmerIterator(seq, 30);

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
            "aggttaacgtcatgcctaagcttacatacg",
            "ggttaacgtcatgcctaagcttacatacgc",
            "gttaacgtcatgcctaagcttacatacgcg",
            "ttaacgtcatgcctaagcttacatacgcgc",
            "taacgtcatgcctaagcttacatacgcgcg",
            "aacgtcatgcctaagcttacatacgcgcgc",
            "acgtcatgcctaagcttacatacgcgcgcc",
            "cgtcatgcctaagcttacatacgcgcgccg",
            "gtcatgcctaagcttacatacgcgcgccgg"
        };

        int index = 0;

        while (kmer.hasNext()) {
            Kmer val = kmer.next();
            assertEquals(expectedKmers[index], val.toString());
            index++;
        }
        
        // test size 60
        kmer = new KmerIterator(seq, 60);

        String[] expected60Kmers = new String[]{
            "agtcgctacatgaactgactacttaggttaacgtcatgcctaagcttacatacgcgcgcc",
            "gtcgctacatgaactgactacttaggttaacgtcatgcctaagcttacatacgcgcgccg",
            "tcgctacatgaactgactacttaggttaacgtcatgcctaagcttacatacgcgcgccgg"
        };
        
        index = 0;
        while (kmer.hasNext()) {
            Kmer val = kmer.next();
            assertEquals(expected60Kmers[index], val.toString());
            index++;
        }
        
    }
}
