/*
 * Copyright (C) 2013 Jordan Fish <fishjord at msu.edu>
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
package edu.msu.cme.rdp.graph.utils;

import edu.msu.cme.rdp.readseq.utils.NuclBinMapping;
import edu.msu.cme.rdp.readseq.utils.ProteinUtils;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Jordan Fish <fishjord at msu.edu>
 */
public class BackTranslationIteratorTest {

    @Test
    public void testBackTranslationIterator() {
        runTest("vh");
        runTest("ag");
        runTest("vhalelv");
        runTest("VH");
    }

    private void runTest(String protKmerStr) {
        char[] protKmer = protKmerStr.toCharArray();
        long expectedTranslations = BackTranslationIterator.countCombinations(protKmer);
        BackTranslationIterator it = new BackTranslationIterator(protKmer);

        byte[] mer;
        String s;
        String prot;
        while ((mer = it.next()) != null) {
            s = NuclBinMapping.fromByteArray(mer);
            prot = ProteinUtils.getInstance().translateToProtein(s, true, 11);
            assertEquals(protKmerStr.toLowerCase(), prot.toLowerCase());
        }

        assertEquals(expectedTranslations, it.getGenerated());
    }
}