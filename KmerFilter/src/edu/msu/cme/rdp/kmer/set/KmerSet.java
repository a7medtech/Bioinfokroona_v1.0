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

import java.io.Serializable;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author fishjord
 */
public class KmerSet<E> implements Serializable {

    private static class Item<E> implements Serializable {

        private final long[] key;
        private E value;
        private Item<E> next;

        public Item(long[] key, E value) {
            this.key = key;
            this.value = value;
        }
    }
    private final int size = 1299827;
    private Item<E>[] items = new Item[size];
    private int elems = 0;

    public void add(long[] key, E val) {
        int bucket = Math.abs( (int) (key[0] % size));

        if (items[bucket] == null) {
            items[bucket] = new Item<E>(key, val);
            elems++;
        } else {
            Item<E> item = items[bucket];
            Item<E> last = null;
            boolean found = false;

            while (item != null) {
                if (Arrays.equals(item.key, key) ) {
                    item.value = val;
                    found = true;
                    break;
                }

                last = item;
                item = item.next;
            }

            if (!found) {
                last.next = new Item<E>(key, val);
                elems++;
            }
        }

    }
    
    public int size() {
        return elems;
    }

    public E get(long[] key) {
        int bucket =  Math.abs( (int) (key[0] % size) );
        Item<E> item = items[bucket];

        while (item != null) {
            if (Arrays.equals(item.key, key)) {
                return item.value;
            }
            item = item.next;
        }

        return null;
    }

    public double getLoad() {
        int load = 0;
        for(int index = 0;index < size;index++) {
            if(items[index] != null) {
                load++;
            }
        }
        
        return (double)load / size;
    }
    
    public void printStats() {
        int collisions = 0;
        int tail = 0;
        int load = 0;
        
        for(int index = 0;index < size;index++) {
            if(items[index] != null) {
                load++;
                
                if(items[index].next != null) {
                    collisions++;
                    Item<E> item = items[index];
                    
                    while(item.next != null) {
                        tail++;
                        item = item.next;
                    }
                }
            }
        }
        
        System.err.println("Load:       " + ((double)load / size));
        System.err.println("Collisions: " + collisions);
        System.err.println("Average tail: " + ((double)tail / collisions));
    }
    
    public boolean containsKey(long[] key) {
        return get(key) != null;
    }

    public Set<long[]> getKeys() {
        Set<long[]> keys = new HashSet();
        Item<E> item = null;

        for (int index = 0; index < items.length; index++) {
            item = items[index];
            while (item != null) {
                keys.add(item.key);
                item = item.next;
            }
        }

        return keys;
    }
}
