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

import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.utils.gbregion.RegionParser;
import edu.msu.cme.rdp.readseq.utils.gbregion.SingleSeqRegion;
import edu.msu.cme.rdp.readseq.writers.FastaWriter;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

/**
 *
 * @author fishjord
 */
public class ReadGenerator {

    public static Map<String, Map<SingleSeqRegion, String>> readRegionMap(String file) throws IOException {
        Map<String, Map<SingleSeqRegion, String>> regionMap = new HashMap();
        BufferedReader reader = new BufferedReader(new FileReader(file));
        String line;

        while ((line = reader.readLine()) != null) {
            String[] lexemes = line.split("\t");

            if (lexemes.length != 3) {
                continue;
            }

            try {
                if (!regionMap.containsKey(lexemes[0])) {
                    regionMap.put(lexemes[0], new HashMap());
                }

                regionMap.get(lexemes[0]).put(RegionParser.parse(lexemes[2]), lexemes[1]);
            } catch (Exception e) {
                throw new RuntimeException("Failed to parse region " + lexemes[2], e);
            }
        }

        reader.close();

        return regionMap;
    }

    public static void main(String[] args) throws Exception {
        if (args.length < 3) {
            System.err.println("USAGE: ReadGenerator <read_size> <coverage> <region_mapping> infiles...");
            System.exit(1);
        }

        final int minReadSize, maxReadSize;
        if (args[0].contains("-")) {
            String[] lexemes = args[0].split("-");
            minReadSize = Integer.valueOf(lexemes[0]);
            maxReadSize = Integer.valueOf(lexemes[1]);
        } else {
            minReadSize = maxReadSize = Integer.valueOf(args[0]);
        }

        final int coverage = Integer.parseInt(args[1]);
        final int averageReadSize = (minReadSize + maxReadSize) / 2;
        final int readSizeDiff = maxReadSize - minReadSize;
        int readSize = maxReadSize;

        Map<String, Map<SingleSeqRegion, String>> regionMap = readRegionMap(args[2]);

        Random rand = new Random();
        Sequence seq;

        System.err.println("Begining read generation, size range = [" + minReadSize + "-" + maxReadSize + "], size diff=" + readSizeDiff + ", coverage = " + coverage);

        FastaWriter out = new FastaWriter(System.out);

        for (int index = 3; index < args.length; index++) {
            SeqReader reader = new SequenceReader(new File(args[index]));

            while ((seq = reader.readNextSequence()) != null) {
                int seqSize = seq.getSeqString().length();
                int requiredReads = (seqSize * coverage) / averageReadSize;
                int maxReadStart = seqSize - maxReadSize;
                Map<SingleSeqRegion, String> regions = regionMap.get(seq.getSeqName());
                if (regions == null) {
                    regions = new HashMap();
                }

                if (seqSize <= maxReadSize) {
                    System.err.println(seq.getSeqName() + "\t" + seqSize + "\t" + maxReadStart + "\ttoo short to generate reads");
                    continue;
                }

                for (int read = 0; read < requiredReads; read++) {
                    if (readSizeDiff != 0) {
                        readSize = rand.nextInt(readSizeDiff) + minReadSize;
                        maxReadStart = seqSize - readSize;
                    }

                    int readStart = rand.nextInt(maxReadStart);
                    String readName = seq.getSeqName() + "_read_" + read;
                    String readSeq = seq.getSeqString().substring(readStart, readStart + readSize);
                    out.writeSeq(readName, readSeq);

                    String gene = "none";

                    for (SingleSeqRegion r : regions.keySet()) {
                        int regionSize = (int) (Math.min(readSize, r.getSeqStop() - r.getSeqStart()) * .5);
                        if (readStart + readSize > r.getSeqStart() && readStart < r.getSeqStop()) {
                            gene = regions.get(r);
                            break;
                        }
                    }

                    System.err.println(seq.getSeqName() + "\t" + readName + "\t" + readStart + "\t" + readSize + "\t" + gene);
                }
            }
        }
        out.close();
    }
}
