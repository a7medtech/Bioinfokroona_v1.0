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
package edu.msu.cme.rdp.graph.search;

import edu.msu.cme.rdp.kmer.NuclKmer;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class AStarSearchNodeTest {

    public AStarSearchNodeTest() {
    }

    @Test
    public void testComparison() {
        //AStarSearchNode discoveredFrom, long kmer, long fwdHash, long rcHash, int stateNo, char state
        AStarNode node1 = new AStarNode(null, new NuclKmer(new char[]{'a', 'a', 'a'}), 0, 0, 15, 'm');
        AStarNode node2 = new AStarNode(null, new NuclKmer(new char[]{'a', 'a', 'a'}), 0, 0, 15, 'd');
        AStarNode node3 = new AStarNode(null, new NuclKmer(new char[]{'a', 'a', 'a'}), 0, 0, 15, 'i');
        AStarNode node4 = new AStarNode(null, new NuclKmer(new char[]{'a', 'c', 'g'}), 0, 0, 15, 'm');
        AStarNode node5 = new AStarNode(null, new NuclKmer(new char[]{'a', 'a', 'a'}), 0, 0, 15, 'm');
        AStarNode node6 = new AStarNode(null, new NuclKmer(new char[]{'c', 'g', 'a'}), 0, 0, 16, 'm');
        AStarNode node7 = new AStarNode(null, new NuclKmer(new char[]{'t', 'a', 'c'}), 0, 0, 50, 'm');

        assertTrue(node1.equals(node1));
        assertTrue(node1.equals(node5));
        assertTrue(node5.equals(node1));

        assertFalse(node1.equals(node2));
        assertFalse(node1.equals(node3));
        assertFalse(node1.equals(node4));
        assertFalse(node1.equals(node6));
        assertFalse(node1.equals(node7));

        assertTrue(node1.compareTo(node2) < 0);
        assertTrue(node1.compareTo(node1) == 0);
        assertTrue(node1.compareTo(node7) < 0);
        assertTrue(node7.compareTo(node1) > 0);

        node1.fval = -1;
        node6.fval = -6;

        node2.fval = -10;

        assertTrue(node2.compareTo(node6) > 0);
        assertTrue(node6.compareTo(node2) < 0);

        assertTrue(node1.compareTo(node6) < 0);
        assertTrue(node6.compareTo(node1) > 0);
    }
}
