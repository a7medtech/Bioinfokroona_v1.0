/*
 * Copyright (C) 2014 wangqion
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

import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author wangqion
 */
public class ProtKmerTest {
    
    public ProtKmerTest() {
    }


    /**
     * Test of shiftRight method, of class ProtKmer.
     */
    @Test
    public void testShift() {
        System.out.println("shift");
        String seq = "tmamrqcalygkggigkstttqnlvaa";
        String[] expectedKmers = new String[]{
            "tmamrqcalygkggigkstt",
            "mamrqcalygkggigksttt",
            "amrqcalygkggigkstttq",
            "mrqcalygkggigkstttqn",
            "rqcalygkggigkstttqnl",
            "qcalygkggigkstttqnlv",
            "calygkggigkstttqnlva",
            "alygkggigkstttqnlvaa"
        };
        int size = 20;
        int position = size;
        Kmer kmer = new ProtKmer(seq.substring(0, size).toCharArray());
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
