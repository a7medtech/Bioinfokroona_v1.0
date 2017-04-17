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
package edu.msu.cme.rdp.graph.search.method.astar;

import edu.msu.cme.rdp.graph.search.AStarNode;
import edu.msu.cme.rdp.graph.search.CandidatePath;
import edu.msu.cme.rdp.kmer.NuclKmer;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class YenTest {

    @Test
    public void testAStarNode() {
        AStarNode node = new AStarNode(null, new NuclKmer(new char[]{ 'a', 'a', 'a' }), 0, 0, 0, 'm');
        AStarNode good = new AStarNode(null, new NuclKmer(new char[]{ 'a', 'a', 'a' }), 0, 0, 0, 'm');
        AStarNode bad = new AStarNode(null, new NuclKmer(new char[]{ 'a', 'a', 't' }), 0, 0, 0, 'd');

        assertEquals(node, good);
        assertFalse(node.equals(bad));
    }

    @Test
    public void testPath() {
        AStarNode node11 = new AStarNode(null, new NuclKmer(new char[]{ 'a', 'a', 'a' }), 0, 0, 0, 'm');
        AStarNode node12 = new AStarNode(node11, new NuclKmer(new char[]{ 'a', 'a', 't' }), 0, 0, 1, 'm');
        AStarNode node13 = new AStarNode(node12, new NuclKmer(new char[]{ 'a', 't', 'a' }), 0, 0, 1, 'i');
        AStarNode node14 = new AStarNode(node13, new NuclKmer(new char[]{ 'a', 't', 'a' }), 0, 0, 2, 'd');

        CandidatePath path1 = new CandidatePath(node14);

        AStarNode node21 = new AStarNode(null, new NuclKmer(new char[]{ 'a', 'a', 'a' }), 0, 0, 0, 'm');
        AStarNode node22 = new AStarNode(node21, new NuclKmer(new char[]{ 'a', 'a', 't' }), 0, 0, 1, 'm');
        AStarNode node23 = new AStarNode(node22, new NuclKmer(new char[]{ 'a', 't', 'a' }), 0, 0, 1, 'i');
        CandidatePath path2 = new CandidatePath(node23);

        AStarNode node31 = new AStarNode(null, new NuclKmer(new char[]{ 'a', 'a', 'a' }), 0, 0, 0, 'm');
        AStarNode node32 = new AStarNode(node31, new NuclKmer(new char[]{ 'a', 'a', 'g' }), 0, 0, 1, 'm');
        AStarNode node33 = new AStarNode(node32, new NuclKmer(new char[]{ 'a', 't', 'a' }), 0, 0, 1, 'd');
        CandidatePath path3 = new CandidatePath(node33);

        AStarNode node41 = new AStarNode(null, new NuclKmer(new char[]{ 'a', 'a', 'g' }), 0, 0, 0, 'm');
        AStarNode node42 = new AStarNode(node41, new NuclKmer(new char[]{ 'a', 'a', 't' }), 0, 0, 1, 'm');
        AStarNode node43 = new AStarNode(node42, new NuclKmer(new char[]{ 'a', 't', 'a' }), 0, 0, 1, 'i');
        CandidatePath path4 = new CandidatePath(node43);

        assertEquals("Path 1's length incorrect", 4, path1.length());
        assertEquals("Path 2's length incorrect", 3, path2.length());
        assertEquals("Path 3's length incorrect", 3, path3.length());

        assertFalse("Path 1 equals path 2", path1.equals(path2));
        assertFalse("Path 1 equals path 3", path1.equals(path3));
        assertFalse("Path 2 equals path 3", path2.equals(path3));

        assertEquals(node11, path1.get(0));
        assertEquals(node12, path1.get(1));
        assertEquals(node13, path1.get(2));
        assertEquals(node14, path1.get(3));

        assertEquals(node21, path2.get(0));
        assertEquals(node22, path2.get(1));
        assertEquals(node23, path2.get(2));

        assertEquals(node31, path3.get(0));
        assertEquals(node32, path3.get(1));
        assertEquals(node33, path3.get(2));

        CandidatePath subPath = path1.subpath(3);

        assertEquals(3, subPath.length());
        assertEquals(subPath, path2);
	assertTrue(subPath.sameBase(path1, 2));
        assertFalse(subPath.equals(path1));
        assertFalse(subPath.equals(path3));

        Map<AStarNode, Set<AStarNode>> disallowedLinks = new HashMap();

        path1.disallowTransition(disallowedLinks, 1);
        path2.disallowTransition(disallowedLinks, 1);
        path3.disallowTransition(disallowedLinks, 1);
        path4.disallowTransition(disallowedLinks, 1);

        assertEquals(2, disallowedLinks.size());

        assertEquals(2, disallowedLinks.get(node12).size());
        assertEquals(2, disallowedLinks.get(node22).size());

        assertEquals(1, disallowedLinks.get(node32).size());
    }
}

