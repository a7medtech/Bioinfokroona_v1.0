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
package edu.msu.cme.rdp.graph.hash;

import java.util.Arrays;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author wangqion
 */
public class CyclicHashTest {

    public CyclicHashTest() {
    }

    /**
     * Test of getInitialHashvalue method, of class CyclicHash.
     */
    @Test
    public void testGetInitialHashvalue() {
        System.out.println("getInitialHashvalue");
        String s = "GGCGC";
        CyclicHash instance = new CyclicHash(5);
        long expResult = 0L;
        long result = instance.getInitialHashvalue(s);
        //System.err.println("initial value=" + result + " " + Long.toBinaryString(result));
        // The hash value changes each time because a random value was chosen, can't assert
        //assertEquals(expResult, result);

    }

    /**
     * Test of updateRight method, of class CyclicHash.
     */
    @Test
    public void testUpdateRight() {
        System.out.println("updateRight");
        int wordlength = 2;
        CyclicHash instance = new CyclicHash(wordlength);
        String s = "GGCGCAGACG";
        long expectedHash = instance.getInitialHashvalue("CA");
        long reversehash = instance.getInitialHashvalue("GG");
        assertFalse((expectedHash == reversehash));
        reversehash = instance.updateRight(reversehash, 2, 1);  // GC
        reversehash = instance.updateRight(reversehash, 2, 0);  // CA
        assertTrue((expectedHash == reversehash));

        wordlength = 5;
        instance = new CyclicHash(wordlength);
        long hashvalue = instance.getInitialHashvalue(s);
        //System.err.println("original GGCGCAGA=" + hashvalue + " " + Long.toBinaryString(hashvalue));

        byte outchar = 2;
        byte inchar = 0;

        long expResult = 0L;
        long updateRight_result = instance.updateRight(hashvalue, outchar, inchar);
        //System.err.println("updateRight=" + updateRight_result + " " + Long.toBinaryString(updateRight_result));


        //should be the same as updateright_hashvalue
        long updateRight_hashvalue = instance.getInitialHashvalue("GCGCA");
        //System.err.println("GCGCA=" + updateRight_hashvalue + " " + Long.toBinaryString(updateRight_hashvalue));

        //should be the same as the original hashvalue
        long updateLeft_result = instance.updateLeft(updateRight_result, inchar, outchar);
        //System.err.println("updateLeft=" + updateLeft_result + " " + Long.toBinaryString(updateLeft_result));


        assertEquals(updateRight_result, updateRight_hashvalue);

         assertEquals(updateLeft_result, hashvalue);

    }

    private String format(long l) {
        String ret = Long.toBinaryString(l);

        if (ret.length() < 64) {
            char[] pad = new char[64 - ret.length()];
            Arrays.fill(pad, '0');
            ret = new String(pad) + ret;
        }

        return ret;
    }

    /**
     * Test of updateLeft method, of class CyclicHash.
     */
    @Test
    public void testUpdateLeft() {
        System.out.println("updateLeft");
        int wordlength = 2;
        CyclicHash instance = new CyclicHash(wordlength);
        long expectedHash = instance.getInitialHashvalue("GG");
        long reversehash = instance.getInitialHashvalue("CA");
        assertFalse((expectedHash == reversehash));
        reversehash = instance.updateLeft(reversehash, 0, 2);  // GC
        reversehash = instance.updateLeft(reversehash, 1, 2);  // GG
        //System.err.println("reversehash =" + reversehash + " " + format(reversehash));
        //System.err.println("expectedHash GG=" + expectedHash + " " + format(expectedHash));

        assertTrue((expectedHash == reversehash));

        instance = new CyclicHash(5);
        String s = "GGCGCAGACG";
        long hashvalue = instance.getInitialHashvalue(s);
        //System.err.println("original GGCGC=" + hashvalue + " " + format(hashvalue));
        byte outchar = 1;
        byte inchar = 0;

        long expResult = 0L;

        long updateLeft_result = instance.updateLeft(hashvalue, inchar, outchar);
        //System.err.println("updateLeft=" + updateLeft_result + " " + format(updateLeft_result));

        // should be the same as updateLeft_result
        long hashvalue_2 = instance.getInitialHashvalue("AGGCG");
        //System.err.println("AGGCG=" + hashvalue_2 + " " + format(hashvalue_2));

    }
}