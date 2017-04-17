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

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class KmerGeneratorTest {
    
    private static final String testSeq = "abcdefghijk";
    private static final int k = 3;
    private static final List<char[]> expected = Arrays.asList(
            new char[] { 'a', 'b', 'c' },
            new char[] { 'b', 'c', 'd' },
            new char[] { 'c', 'd', 'e' },
            new char[] { 'd', 'e', 'f' },
            new char[] { 'e', 'f', 'g' },
            new char[] { 'f', 'g', 'h' },
            new char[] { 'g', 'h', 'i' },
            new char[] { 'h', 'i', 'j' },
            new char[] { 'i', 'j', 'k' }
            );
    
    /**
     * Test of getKmers method, of class KmerGenerator.
     */
    @Test
    public void testGetKmers() {
        List<char[]> result = KmerGenerator.getKmers(testSeq, k);

        for(int index = 0;index < result.size();index++) {
            assertArrayEquals(expected.get(index), result.get(index));
        }
        
        assertEquals(expected.size(), result.size());
    }

    /**
     * Test of iterator method, of class KmerGenerator.
     */
    @Test
    public void testIterator() {
        KmerGenerator instance = new KmerGenerator(testSeq, k);
        int index = 0;
                
        Iterator<char[]> result = instance.iterator();
        while(result.hasNext()) {
            assertArrayEquals(expected.get(index++), result.next());
        }
        
        assertEquals(expected.size(), index);
    }
    /**
     * Test of iterator method, of class KmerGenerator.
     */
    @Test
    public void testIterable() {
        KmerGenerator instance = new KmerGenerator(testSeq, k);
        int index = 0;
        
        for(char[] kmer : instance) {
            assertArrayEquals(expected.get(index++), kmer);
        }
        
        assertEquals(expected.size(), index);
    }
}
