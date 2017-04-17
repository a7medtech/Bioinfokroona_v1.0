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
package edu.msu.cme.rdp.kmer;

import edu.msu.cme.rdp.readseq.utils.NuclBinMapping;
import java.util.Random;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class KmerTest {

    @Test
    public void testKmer() {

        Random rand = new Random();

        for (int length = 1; length < 15; length++) {
           for (int iteration = 0; iteration < 10; iteration++) {
                StringBuilder strKmer = new StringBuilder();
                for (int index = 0; index < length; index++) {
                    strKmer.append(NuclBinMapping.intToChar[rand.nextInt(4)]);
                }
       
                NuclKmer kmer = new NuclKmer(strKmer.toString().toCharArray());
                assertEquals(strKmer.toString(), kmer.toString());

                StringBuilder rightShift = new StringBuilder(strKmer.toString());
                
                char shift = NuclBinMapping.intToChar[rand.nextInt(4)];
                rightShift.insert(0, shift);
                rightShift.deleteCharAt(rightShift.length() - 1);

                Kmer newKmer = kmer.shiftRight(shift);
                assertEquals("Shift Right Original kmer: " + strKmer + ", shift=" + shift,rightShift.toString(), newKmer.toString());
                                
                // test decodedLong()
                String decodedStr = newKmer.decodeLong(newKmer.getLongKmers());
                assertEquals(rightShift.toString(), decodedStr);
                
                StringBuilder leftShift = new StringBuilder(strKmer.toString());
                shift = NuclBinMapping.intToChar[rand.nextInt(4)];

                leftShift.append(shift);
                leftShift.deleteCharAt(0);

                newKmer = kmer.shiftLeft(shift);
                assertEquals("Shift Left Original kmer: " + strKmer + ", length: " + length + ", shift=" + shift, leftShift.toString(), newKmer.toString());
           }  
        
        }
    }
    
    /**
     * Test of shiftRight method, of class ProtKmer.
     */
    @Test
    
    public void testShift() {
        System.out.println("shift");
        String seq = "agtcgctacatgaactgactacttaggttaacgtcatgcctaagcttacatacgcgcgccgg";
        String[] expectedKmers = new String[]{
            "agtcgctacatgaactgactacttaggttaacgtcatgcctaagcttacatacgcgcgcc",
            "gtcgctacatgaactgactacttaggttaacgtcatgcctaagcttacatacgcgcgccg",
            "tcgctacatgaactgactacttaggttaacgtcatgcctaagcttacatacgcgcgccgg"
        };
        int size = 60;
        int position = size;
        Kmer kmer = new NuclKmer(seq.substring(0, size).toCharArray());
        assertEquals(expectedKmers[position - size], kmer.toString());
        // shift left
        while ( position < seq.length() ) {
            kmer = kmer.shiftLeft(seq.charAt(position));   
            assertEquals(expectedKmers[position - size +1], kmer.toString());
            position++;
        }
        
        // shift right
        while ( (position -size -1) >= 0 ) {
            kmer = kmer.shiftRight(seq.charAt(position -size -1));   
            assertEquals(expectedKmers[position - size -1], kmer.toString());
            position--;
        }
        
    }
    
}
