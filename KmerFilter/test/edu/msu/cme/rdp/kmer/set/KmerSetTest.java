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

import java.util.Arrays;
import java.util.Set;
import java.util.HashMap;
import java.util.Map;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class KmerSetTest {
    
    @Test
    public void testKmerSet() {
        long[] values = new long[]{
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
            197062507064968390L,
            197062507064968390L
        };
        
        
        KmerSet<Integer> set = new KmerSet();
        
        for(int index = 0;index < values.length;index +=2) {
            long[] temp = new long[2];
            temp[0] = values[index];
            temp[1] = values[index +1];
            set.add(temp, index);
        }
        
        assertEquals(values.length/2, set.size());
        
        for(int index = 0;index < values.length;index +=2) {
            long[] temp = new long[2];
            temp[0] = values[index];
            temp[1] = values[index +1];
            assertEquals(index, (int)set.get(temp));
        }
        
    }
}
