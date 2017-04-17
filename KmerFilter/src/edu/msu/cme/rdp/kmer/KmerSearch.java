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
import edu.msu.cme.rdp.kmer.trie.KmerGenerator;
import edu.msu.cme.rdp.readseq.SequenceType;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.utils.IUBUtilities;
import edu.msu.cme.rdp.readseq.utils.ProteinUtils;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import edu.msu.cme.rdp.readseq.writers.FastaWriter;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

/**
 *
 * @author fishjord
 */
public class KmerSearch {

    private static final Options options = new Options();

    static {
        options.addOption("o", "out", true, "Redirect output to file");
        options.addOption("t", "transl-table", true, "Translation table to use when translating nucleotide to protein sequences");
        options.addOption("c", "correct", false, "Assert sequences are in the correct orientation and reading frame");
    }

    public static void main(String[] args) throws IOException {
        KmerTrie kmerTrie = null;
        SeqReader queryReader = null;
        SequenceType querySeqType = SequenceType.Unknown;
        FastaWriter out = null;
        boolean exhaustive = true;
        boolean translQuery = false;
        int wordSize = -1;
        int translTable = 11;

        try {
            CommandLine cmdLine = new PosixParser().parse(options, args);
            args = cmdLine.getArgs();

            if (args.length != 3) {
                throw new Exception("Unexpected number of arguments");
            }

            if (cmdLine.hasOption("out")) {
                out = new FastaWriter(cmdLine.getOptionValue("out"));
            } else {
                out = new FastaWriter(System.out);
            }

            if (cmdLine.hasOption("correct")) {
                exhaustive = false;
            } else {
                exhaustive = true;
            }

            if (cmdLine.hasOption("transl-table")) {
                translTable = Integer.valueOf(cmdLine.getOptionValue("transl-table"));
            }

            File trainingFile = new File(args[0]);
            wordSize = Integer.valueOf(args[1]);
            File queryFile = new File(args[2]);

            querySeqType = SeqUtils.guessSequenceType(queryFile);
            queryReader = new SequenceReader(new File(args[2]));

            kmerTrie = KmerTrie.buildTrie(new SequenceReader(trainingFile), wordSize);

            if (querySeqType == SequenceType.Protein && kmerTrie.getTreeSeqType() == SequenceType.Nucleotide) {
                throw new Exception("Trie is made of nucleotide sequences but the query sequences are protein");
            }

            if (querySeqType == SequenceType.Nucleotide && kmerTrie.getTreeSeqType() == SequenceType.Protein) {
                translQuery = true;
                System.err.println("Query sequences are nucleotide but trie is in protein space, query sequences will be translated");
            }

            if(querySeqType == SequenceType.Protein && exhaustive) {
                System.err.println("Cannot do an exaustive search with protein sequences, disabling");
                exhaustive = false;
            }

        } catch (Exception e) {
            new HelpFormatter().printHelp("KmerSearch <ref_file> <word_size> <query_file>", options);
            System.err.println(e.getMessage());
            System.exit(1);
        }

        long startTime = System.currentTimeMillis();
        long seqCount = 0, passedCount = 0;

        Sequence querySeq;

        while ((querySeq = queryReader.readNextSequence()) != null) {
            seqCount++;

            List<Sequence> testSequences;

            if (!exhaustive) {
                if (translQuery) {
                    testSequences = Arrays.asList(new Sequence(querySeq.getSeqName(), "", ProteinUtils.getInstance().translateToProtein(querySeq.getSeqString(), true, translTable)));
                } else {
                    testSequences = Arrays.asList(querySeq);
                }
            } else {
                if (translQuery) {
                    testSequences = ProteinUtils.getInstance().allTranslate(querySeq);
                } else {
                    testSequences = Arrays.asList(querySeq,
                            new Sequence(querySeq.getSeqName(), "", IUBUtilities.reverseComplement(querySeq.getSeqString())));
                }
            }

            boolean passed = false;
            for (Sequence seq : testSequences) {
                for (char[] kmer : KmerGenerator.getKmers(seq.getSeqString(), wordSize)) {
                    if (kmerTrie.contains(kmer) != null) {
                        passed = true;
                        break;
                    }
                }
                
                if (passed) {
                    out.writeSeq(seq);
                    passedCount++;
                    break;
                }
            }
        }
        System.err.println("Processed: " + seqCount);
        System.err.println("Passed: " + passedCount);
        System.err.println("Failed: " + (seqCount - passedCount));
        System.err.println("Time: " + (System.currentTimeMillis() - startTime) + " ms");
    }
}
