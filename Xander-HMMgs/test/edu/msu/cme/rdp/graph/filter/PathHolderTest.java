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
package edu.msu.cme.rdp.graph.filter;

import edu.msu.cme.rdp.kmer.Kmer;
import edu.msu.cme.rdp.kmer.NuclKmer;
import edu.msu.cme.rdp.readseq.utils.NuclBinMapping;
import java.util.Random;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class PathHolderTest {

    @Test
    public void testPathHolder() {

        Random r = new Random();

        for (int i = 0; i < 100; i++) {
            StringBuilder expected = new StringBuilder();
            int l = r.nextInt(200) + 150;

            PathHolder holder = new PathHolder();
            for (int index = 0; index < l; index++) {
                byte add = (byte) r.nextInt(4);
                holder.push(add);
                expected.append(NuclBinMapping.intToChar[add]);
            }

            assertEquals(expected.toString(), new String(holder.toCharArray()));
            StringBuilder buf = new StringBuilder();
            for (int index = 0; index < holder.size(); index++) {
                buf.append(NuclBinMapping.intToChar[holder.get(index)]);
            }
            assertEquals(expected.toString(), buf.toString());
        }
    }

    @Test
    public void specificTest() {
        PathHolder path = new PathHolder();
        char[] expected = "aacgttugta".toCharArray();
        String expectedString = new String(expected).replace('u', 't');
        String expectedReverse = new StringBuilder(expectedString).reverse().toString();

        path.init(new NuclKmer(expectedReverse.toCharArray()));
        assertEquals(expectedReverse, new String(path.toCharArray()));
        assertEquals(expectedReverse, pathToStringViaGet(path));

        path.init(new NuclKmer(expected));
        assertEquals(expectedString, new String(path.toCharArray()));
        assertEquals(expectedString, pathToStringViaGet(path));

        int substr = expected.length - 1;
        while(path.size() > 0) {
            char exp = expectedString.charAt(substr);
            String e = expectedString.substring(0, substr--);
            char removed = NuclBinMapping.intToChar[path.remove()];

            assertEquals(exp, removed);
            assertEquals(e, new String(path.toCharArray()));
            assertEquals(e, pathToStringViaGet(path));
        }
    }

    @Test
    public void pathFromKmerTest() {
        Kmer kmer;
        PathHolder path = new PathHolder();

        kmer = NuclKmer.randomKmer(30);
        path.init(kmer);
        assertEquals(kmer.toString(), new String(path.toCharArray()));

        kmer = NuclKmer.randomKmer(60);
        path.init(kmer);
        assertEquals(kmer.toString(), new String(path.toCharArray()));
    }

    private String pathToStringViaGet(PathHolder path) {
        StringBuilder ret = new StringBuilder();
        for(int index = 0;index < path.size();index++) {
            ret.append(NuclBinMapping.intToChar[path.get(index)]);
        }

        return ret.toString();
    }
}
