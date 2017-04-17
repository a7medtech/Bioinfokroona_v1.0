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

package edu.msu.cme.rdp.kmer;

import edu.msu.cme.rdp.kmer.trie.KmerTrie;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import java.io.File;
import java.io.PrintWriter;

/**
 *
 * @author fishjord
 */
public class TrieInfo {

    public static void main(String[] args) throws Exception {
        if(args.length != 2 && args.length != 3) {
            System.err.println("USAGE: TrieInfo <ref_file> <word_size> [word histogram]");
            System.exit(1);
        }

        int k = Integer.valueOf(args[1]);
        System.err.println("Word size is " + k);

        System.err.println("Building trie...");
        long startTime = System.currentTimeMillis();
        KmerTrie kmerTrie = KmerTrie.buildTrie(new SequenceReader(new File(args[0])), k);
        long stop = System.currentTimeMillis();
        System.err.println("Trie built in " + (System.currentTimeMillis() - startTime) + " ms");

        System.err.println("Trie contains " + kmerTrie.getTreeSeqType() + " kmers");
        System.err.println("Trie contains " + kmerTrie.countNodes() + " nodes");
        System.err.println("Trie contains " + kmerTrie.uniqueWords() + " unique words");

        System.out.println(k + "\t" + kmerTrie.countNodes() + "\t" + kmerTrie.uniqueWords() + "\t" + (stop - startTime));

        if(args.length == 3) {
            PrintWriter out = new PrintWriter(args[2]);
            kmerTrie.printWordHistogram(out);
            out.close();
        }
    }
}
