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

package edu.msu.cme.rdp.kmer.cli;

import edu.msu.cme.rdp.kmer.KmerFilter;
import edu.msu.cme.rdp.kmer.KmerSearch;
import edu.msu.cme.rdp.kmer.ReadGenerator;
import edu.msu.cme.rdp.kmer.TrieInfo;
import java.util.Arrays;

/**
 *
 * @author fishjord
 */
public class Main {
    public static void main(String[] args) throws Exception {
        String usage = "Commands: kmer_search, kmer_filter, fast_kmer_filter, kmer_coverage, shared_kmer, read_generator, trie_info";
        if(args.length == 0) {
            System.err.println(usage);
            return;
        }

        String cmd = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);

        if(cmd.equals("kmer_search")) {
            KmerSearch.main(newArgs);
        } else if(cmd.equals("kmer_filter")) {
            KmerFilter.main(newArgs);
        } else if(cmd.equals("fast_kmer_filter")) {
            FastKmerFilter.main(newArgs);
        } else if(cmd.equals("kmer_coverage")) {
            KmerCoverage.main(newArgs);
        } else if(cmd.equals("read_generator")) {
            ReadGenerator.main(newArgs);
        } else if(cmd.equals("trie_info")) {
            TrieInfo.main(newArgs);
        } else if(cmd.equals("shared_kmer")) {
            SharedKmers.main(newArgs);
        } else {
            System.err.println(usage);
        }
    }
}
