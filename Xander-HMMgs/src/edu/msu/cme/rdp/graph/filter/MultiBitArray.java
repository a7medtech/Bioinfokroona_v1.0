/*
 * Copyright (C) 2014 gilmanma
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

import java.util.BitSet;
import java.util.List;
import java.util.ArrayList;
import java.io.Serializable;

/**
 * A class similar to a BitSet that allows multiple bits per bucket
 * 
 * @author gilmanma
 * 
 */
public class MultiBitArray implements Serializable {
    
    /**
     * Log2 of maximum number of BitSets
     */
    private static final int MAX_HASH_SIZE_LOG2 = 64;
    
    /**
     * Log2 of maximum size of individual BitSet
     */
    private static final int MAX_BITSET_SIZE_LOG2 = 30;

    public static final long serialVersionUID = -1919209030378954154L;
    /**
     * Size of an individual BitSet
     */
    private final int bitSetSize;
    private final long bitSetMask;
    private final int bitSetSizeLog2;

    private List<BitSet> bitSetList;

    /**
     * Number of bits per bucket
     */
    private int numBits;

    public MultiBitArray(int hashSizeLog2, int bitSetSizeLog2) {
        this(hashSizeLog2, bitSetSizeLog2, 1);
    }

    public MultiBitArray(int hashSizeLog2, int bitSetSizeLog2, int numBits) {
        if(bitSetSizeLog2 > MAX_BITSET_SIZE_LOG2) {
            throw new IllegalArgumentException("Can't have bitset larger than 2^" + MAX_BITSET_SIZE_LOG2);
        }
        if(hashSizeLog2 > MAX_HASH_SIZE_LOG2) {
            throw new IllegalArgumentException("Can't have filter larger than 2^" + MAX_HASH_SIZE_LOG2);
        }

        this.bitSetSizeLog2 = bitSetSizeLog2;
        this.numBits = numBits;

        this.bitSetSize = (1 << this.bitSetSizeLog2);
        this.bitSetMask = this.bitSetSize - 1;

        int bitSetCount = 1;
        if(hashSizeLog2 > bitSetSizeLog2) {
            bitSetCount = (1 << (hashSizeLog2-bitSetSizeLog2));
        }
        bitSetCount *= numBits;

        bitSetList = new ArrayList<BitSet>(bitSetCount);
        for(int i = 0; i < bitSetCount; i++) {
            bitSetList.add(i, new BitSet(bitSetSize));
        }
    }

    /** 
     * @param bit   bit for which to find the offset of
     * 
     * @return position of bit within the BitSet
     */
    protected int getOffset(long bit) {
        return (int) (bit & bitSetMask);
    }

    /**
     * @param bit   bit for which to find the BitSet of
     * 
     * @return number of BitSet that contains the bit
     */
    protected int getSetNum(long bit) {             
        return (int) (bit >>> bitSetSizeLog2) * numBits;
    }

    
     /**
     * Increments count at given position
     * 
     * @param bit   which bit to increment
     * @return      whether the bit was previously set
     */
    public boolean setBit(long bit) {
        int bitSetOffset = getOffset(bit);
        int bitSetNum = getSetNum(bit);
        boolean wasSet = false;
        
        //  we can speed up count 1 and 2 more than 40% compared to the for loop below
       if ( numBits == 1){ 
            boolean val = bitSetList.get(bitSetNum ).get(bitSetOffset);
            wasSet |= val;
            bitSetList.get(bitSetNum).set(bitSetOffset, true);
        }else if (numBits == 2){ //this is faster than the for loop below 
            boolean val = bitSetList.get(bitSetNum ).get(bitSetOffset);
            wasSet |= val;
            if ( !val){
                bitSetList.get(bitSetNum).set(bitSetOffset, true);
            }else{                
                val = bitSetList.get(bitSetNum +1 ).get(bitSetOffset);
                wasSet |= val;
                if ( !val ){
                    bitSetList.get(bitSetNum).set(bitSetOffset, false);
                    bitSetList.get(bitSetNum+1).set(bitSetOffset, true);
                }
            }
        }else { // this loop works for all the counts but slow
        // if the count already reached the max, don't do anything
            boolean isBitSet = true;
            for(int i = 0; i < numBits; ++i) {
                isBitSet &= bitSetList.get(bitSetNum + i).get(bitSetOffset);
                if ( !isBitSet){
                    break;
                }
            }
            if ( isBitSet){
                 return wasSet;
            }

        // carry = whether the previous bit rolled over into this one
            // for first bit, this is just adding one to it
            boolean carry = true;

            for(int i = 0; i < numBits; ++i) {
                boolean val = bitSetList.get(bitSetNum + i).get(bitSetOffset);
                wasSet |= val;
                bitSetList.get(bitSetNum + i).set(bitSetOffset, (val ^ carry));
                carry &= val;
            }
        }
        return !wasSet;
    }
    
    
    /**
     * 
     * @param bit   which bit to check
     * @return      whether the bit is set
     */
    public boolean isSet(long bit) {
        int bitSetOffset = getOffset(bit);
        int bitSetNum = getSetNum(bit);
        
        boolean isSet = false;
        for(int i = 0; i < numBits; ++i) {
            isSet |= bitSetList.get(bitSetNum + i).get(bitSetOffset);
        }
        return isSet;
    }

    /**
     * @param bit   bit to get the count of
     * @return      count at the given bit
     */
    public int getCount(long bit) {
        int bitSetOffset = getOffset(bit);
        int bitSetNum = getSetNum(bit);

        int count = 0;
        for(int i = 0; i < numBits; ++i) {
            count += (bitSetList.get(bitSetNum + i).get(bitSetOffset) ? 1 << i : 0);
        }
        return count;
    }

    /** 
     * "Collapses" a MultiBitArray with multiple bits per bucket to one with one
     * bit per bucket
     * 
     * @param cutoff    how many bits must be in the bucket for the resulting
     *                  MultiBitArray to have a 1 at that position
     */
    public void collapse(int cutoff) {
        if(numBits == 1) {
            return;
        }

        List<BitSet> newList = new ArrayList<BitSet>(bitSetList.size()/numBits);

        for(int i = 0; i < bitSetList.size()/numBits; ++i) {
            BitSet baseBitSet = bitSetList.get(i*numBits);
            long pos = (long) i * bitSetSize;
            for(int j = 0; j < bitSetSize; ++j) {
                baseBitSet.set(j, getCount(pos+j) >= cutoff);
            }
            newList.add(baseBitSet);
        }
        numBits = 1;
        bitSetList = newList;
    }
    
    public int getBitSetSize() {
        return bitSetSize;
    }
    
    public long getBitSetMask() {
        return bitSetMask;
    }
    
    /**
     * Returns the effective number of BitSets 
     * (i.e. How many there are in the list divided by the number of buckets,
     * since there is one BitSet per bucket)
     * 
     * @return  number of BitSets
     */
    public int getNumBitSets() {
        return bitSetList.size()/numBits;
    }
}
