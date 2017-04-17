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

import edu.msu.cme.rdp.readseq.SequenceType;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author fishjord
 */
public class KmerTrie implements Serializable {

    private static final byte[] alphaMap = new byte[127];
    private static final int rnaAlphaSize = 4;
    private static final int proteinAlphaSize;

    static {
        Arrays.fill(alphaMap, (byte) -1);
        alphaMap['a'] = alphaMap['A'] = 0;
        alphaMap['c'] = alphaMap['C'] = 1;
        alphaMap['g'] = alphaMap['G'] = 2;
        alphaMap['t'] = alphaMap['T'] = alphaMap['u'] = alphaMap['U'] = 3;

        byte nextIndex = 4;
        for (Character c : SeqUtils.proteinAlphabet) {
            c = Character.toLowerCase(c);
            if (alphaMap[c] == -1) {
                alphaMap[c] = alphaMap[Character.toUpperCase(c)] = nextIndex++;
            }
        }

        proteinAlphaSize = nextIndex;
    }

    public static class RefPos {
        public Integer modelPos;
        public String seqid;
    }

    private abstract static class TrieNode {
    }

    private static class TrieInteriorNode extends TrieNode {

        TrieNode[] children;
    }

    public static class TrieLeaf extends TrieNode {

        private int frame = 0;
        private int count = 0;
        private int queryCount = 0;
        private Map<Integer, Set<RefPos>> refSetToModelStarts = new HashMap();
        private Set<String> seqids = new HashSet();

        public TrieLeaf(int frame) {
            this.frame = frame;
        }

        public int getFrame() {
            return frame;
        }

        public int getCount() {
            return count;
        }

        public Set<Integer> getRefSets() {
            return Collections.unmodifiableSet(refSetToModelStarts.keySet());
        }

        public Set<RefPos> getModelStarts(int refSet) {
            return Collections.unmodifiableSet(refSetToModelStarts.get(refSet));
        }

        public synchronized void incQueryCount() {
            queryCount++;
        }

        public int getQueryCount() {
            return queryCount;
        }

        void inc() {
            count++;
        }
    }
    private final TrieNode root;
    private final boolean isProtein;
    private final int k;
    private int seqCount;

    public KmerTrie(int k, boolean isProtein) {
        this.k = k;
        this.isProtein = isProtein;

        TrieInteriorNode rootNode = new TrieInteriorNode();
        rootNode.children = new TrieNode[isProtein ? proteinAlphaSize : rnaAlphaSize];
        root = rootNode;
    }

    public int getWordSize() {
        return k;
    }

    public static KmerTrie buildTrie(SeqReader reader, int k) throws IOException {
        Sequence seq = reader.readNextSequence();

        if (seq == null) {
            throw new IOException("Sequence file contains no sequences");
        }
        KmerTrie trie = new KmerTrie(k, SeqUtils.guessSequenceType(seq) == SequenceType.Protein);

        do {
            trie.addSequence(seq);
        } while ((seq = reader.readNextSequence()) != null);

        return trie;
    }

    public static KmerTrie buildTrieFromAligned(SeqReader reader, final int k) throws IOException {
        Sequence seq = reader.readNextSequence();

        if (seq == null) {
            throw new IOException("Sequence file contains no sequences");
        }
        KmerTrie trie = new KmerTrie(k, SeqUtils.guessSequenceType(seq) == SequenceType.Protein);

        do {
            if (seq.getSeqName().startsWith("#")) {
                continue;
            }

            trie.addModelSequence(seq);
        } while ((seq = reader.readNextSequence()) != null);

        return trie;
    }

    private void addKmer(char[] kmer, int frame, int modelPos, int refSet, String seqid) {
        if (kmer.length != k) {
            throw new IllegalArgumentException(new String(kmer) + "'s length doesn't match expected (" + k + ")");
        }
        addKmer(kmer, 0, frame, modelPos, refSet, seqid);
    }

    private void addKmer(char[] kmer, int offset, int frame, int modelPos, int refSet, String seqid) {
        if (offset + k > kmer.length) {
            throw new IllegalArgumentException("Array offset doesn't leave enough room for a full kmer");
        }

        TrieNode curr = root;
        TrieInteriorNode node;
        int end = offset + k;

        for (int index = offset; index < end; index++) {
            char c = kmer[index];
            node = (TrieInteriorNode) curr;

            if(alphaMap[c] == -1) {
                throw new IllegalArgumentException("kmer contains unmappable character " + c);
            }

            if (node.children[alphaMap[c]] == null) {
                if (index + 1 == end) {
                    curr = new TrieLeaf(frame);
                    node.children[alphaMap[c]] = curr;
                } else {
                    TrieInteriorNode interNode = new TrieInteriorNode();
                    interNode.children = new TrieNode[isProtein ? proteinAlphaSize : rnaAlphaSize];
                    node.children[alphaMap[c]] = interNode;

                    curr = interNode;
                }
            } else {
                curr = node.children[alphaMap[c]];
            }
        }

        TrieLeaf leaf = (TrieLeaf) curr;
        leaf.count++;
        if(!leaf.refSetToModelStarts.containsKey(refSet)) {
            leaf.refSetToModelStarts.put(refSet, new HashSet());
        }
        RefPos pos = new RefPos();
        pos.modelPos = modelPos;
        pos.seqid = seqid;
        leaf.refSetToModelStarts.get(refSet).add(pos);
    }

    public void addSequence(Sequence seq) {
        addSequence(seq, 0);
    }

    public void addSequence(Sequence seq, int refSet) {

        int frame = 0;
        char[] bases = seq.getSeqString().toCharArray();

        for (int index = 0; index <= bases.length - k; index++) {
            if (isProtein) {
                addKmer(bases, index, -1, -1, refSet, seq.getSeqName());
            } else {
                addKmer(bases, index, frame, -1, refSet, seq.getSeqName());
                frame = (frame > 1) ? 0 : frame++;
            }
        }
        seqCount++;
    }

    public void addModelSequence(Sequence seq) {
        addModelSequence(seq, 0);
    }

    public void addModelSequence(Sequence seq, int refSet) {

        int frame = 0;
        ModelPositionKmerGenerator kmers = new ModelPositionKmerGenerator(seq.getSeqString(), k, getTreeSeqType());

        for (char[] kmer : kmers) {
            addKmer(kmer, -1, kmers.getModelPosition(), refSet, seq.getSeqName());

        }
        seqCount++;
    }

    public TrieLeaf contains(char[] kmer) {
        if (kmer.length != k) {
            throw new IllegalArgumentException(new String(kmer) + "'s length doesn't match expected (" + k + ")");
        }

        return contains(kmer, 0);
    }

    public TrieLeaf contains(char[] kmer, int offset) {
        return containsInternal(kmer, offset);
    }

    private TrieLeaf containsInternal(char[] kmer, int offset) {
        if (offset + k > kmer.length) {
            throw new IllegalArgumentException("Array offset doesn't leave enough room for a full kmer " + offset + ", " + k + " " + kmer.length);
        }
        TrieNode curr = root;
        TrieInteriorNode node;

        for (int index = offset; index < offset + k; index++) {
            char c = kmer[index];
            node = (TrieInteriorNode) curr;

            if (node.children[alphaMap[c]] == null) {
                return null;
            }
            curr = node.children[alphaMap[c]];
        }

        return ((TrieLeaf) curr);
    }

    public SequenceType getTreeSeqType() {
        return isProtein ? SequenceType.Protein : SequenceType.Nucleotide;
    }

    public int getSeqCount() {
        return seqCount;
    }

    public int uniqueWords() {
        return countUniqueWords(root);
    }

    private int countUniqueWords(TrieNode node) {
        int ret = 0;

        if (node instanceof TrieLeaf) {
            return 1;
        } else {
            TrieInteriorNode curr = (TrieInteriorNode) node;
            TrieNode child = null;
            for (int index = 0; index < curr.children.length; index++) {
                child = curr.children[index];

                if (child != null) {
                    ret += countUniqueWords(child);
                }
            }
        }

        return ret;
    }

    public int countNodes() {
        return countNodes(root);
    }

    private int countNodes(TrieNode node) {
        int ret = 0;

        if (node instanceof TrieLeaf) {
            return 1;
        } else {
            TrieInteriorNode curr = (TrieInteriorNode) node;
            TrieNode child = null;
            for (int index = 0; index < curr.children.length; index++) {
                child = curr.children[index];

                if (child != null) {
                    ret += countNodes(child);
                }
            }
        }

        return ret + 1;
    }

    public void printWordHistogram(PrintWriter out) {
        printWordHistogram(root, new char[k], 0, out);
    }

    private int printWordHistogram(TrieNode node, char[] currWord, int depth, PrintWriter out) {
        int ret = 0;

        if (node instanceof TrieLeaf) {
            out.print(currWord);
            out.println("\t" + ((TrieLeaf) node).count);
        } else {
            TrieInteriorNode curr = (TrieInteriorNode) node;
            TrieNode child = null;
            for (int index = 0; index < curr.children.length; index++) {
                child = curr.children[index];

                if (child != null) {
                    char c = 'n';

                    for (int alphaIndex = 0; alphaIndex < proteinAlphaSize; alphaIndex++) {
                        if (alphaMap[alphaIndex] == index) {
                            c = (char) alphaIndex;
                        }
                    }

                    currWord[depth] = c;
                    printWordHistogram(child, currWord, depth + 1, out);
                }
            }
        }

        return ret + 1;
    }
}
