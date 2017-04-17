/*
 * Copyright (C) 2012 Michigan State University <rdpstaff at msu.edu>
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
package org.apache.hadoop.mapred;

import edu.msu.cme.rdp.hadoop.distance.mapred.keys.DistanceAndComparison;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.DataInputBuffer;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.compress.GzipCodec;
import org.apache.hadoop.util.GenericOptionsParser;

/**
 *
 * @author Jordan Fish <fishjord at msu.edu>
 */
public class SummarizeAttempts {

    private static void guessSplitRanges() throws IOException {

        long numSeqs = 2334222;
        int inMemSeqs = 2000;
        long rangesPerSplit = 100;
        // figure out number of ranges per split to use
        long numRanges = (numSeqs / inMemSeqs) * (numSeqs / inMemSeqs) / 2;
        long numSplits = Math.max(numRanges / rangesPerSplit, 1);

        int splits = 0;
        int i = 0;
        PrintStream out = new PrintStream("/home/fishjord/tasks_to_splits.txt");
        for (long x = 0; x < numSeqs; x += inMemSeqs) {
            for (long y = x; y < numSeqs; y += inMemSeqs) {
                if (splits >= rangesPerSplit) {

                    i++;
                    splits = 0;
                }
                out.println(i + "\t" + splits + "\t" + x + "\t" + y + "\t" + inMemSeqs);
                splits++;
            }
        }
        out.close();
    }

    public static void main(String[] args) throws Exception {
        //guessSplitRanges();
        Configuration conf = new Configuration();
        args = new GenericOptionsParser(conf, args).getRemainingArgs();

        if (args.length != 1) {
            System.err.println("USAGE: ReadPartialResultFile <job_dir>");
            System.exit(1);
        }

        File jobDir = new File(args[0]);
        String jobName = jobDir.getName();
        if (!jobName.startsWith("job_")) {
            throw new IOException("Expected job dir name to start with 'job_'");
        }

        FileSystem fs = FileSystem.getLocal(conf);

        DataInputBuffer keyBuf = new DataInputBuffer();
        DataInputBuffer valBuf = new DataInputBuffer();

        System.out.println("path\tattempt\tnum_dists\tsorted?\tminX\tmaxX\tminY\tmaxY\tavg_dist\ttime(s)");
        GzipCodec codec = new GzipCodec();
        codec.setConf(conf);

        for (File attemptDir : jobDir.listFiles()) {
            if (!attemptDir.isDirectory() || !attemptDir.getName().startsWith("attempt")) {
                continue;
            }

            String attemptNum = attemptDir.getName().split("_")[4];

            int numDirs = 0;
            for (File f : attemptDir.listFiles()) {
                if (f.isDirectory()) {
                    numDirs++;
                }
            }

            if (numDirs != 1) {
                System.err.println(attemptDir + " does not have exactly one sub directory");
                continue;
            }

            File gzFileOut = new File(attemptDir, "output/file.out.gz");
            File fileOut = new File(attemptDir, "output/file.out");
            boolean zipped;
            Path path;

            if (fileOut.exists()) {
                zipped = false;
                path = new Path(fileOut.getAbsolutePath());
            } else if (gzFileOut.exists()) {
                zipped = true;
                path = new Path(gzFileOut.getAbsolutePath());
            } else {
                System.err.println("Couldn't find output directory in " + attemptDir);
                continue;
            }

            try {
                IFile.Reader<DistanceAndComparison, NullWritable> reader = new IFile.Reader<DistanceAndComparison, NullWritable>(conf, fs, path, (zipped) ? codec : null);

                int numDists = 0;
                int minX = Integer.MAX_VALUE;
                int maxX = 0;
                int minY = Integer.MAX_VALUE;
                int maxY = 0;
                long totalDist = 0;
                DistanceAndComparison ds = new DistanceAndComparison();

                int lastDist = -1;
                boolean sorted = true;

                long startTime = System.currentTimeMillis();
                while (reader.next(keyBuf, valBuf)) {
                    ds.readFields(keyBuf);
                    numDists++;

                    if (ds.first < minX) {
                        minX = ds.first;
                    }

                    if (ds.first > maxX) {
                        maxX = ds.first;
                    }

                    if (ds.second < minY) {
                        minY = ds.second;
                    }

                    if (ds.second > maxY) {
                        maxY = ds.second;
                    }

                    if (lastDist > ds.distance) {
                        sorted = false;
                    }

                    lastDist = ds.distance;
                    totalDist += ds.distance;
                }

                reader.close();
                System.out.println(path + "\t" + attemptNum + "\t" + numDists + "\t" + sorted + "\t" + minX + "\t" + maxX + "\t" + minY + "\t" + maxY + "\t" + ((float) totalDist / numDists) + "\t" + (System.currentTimeMillis() - startTime) / 1000.f);
            } catch (Exception e) {
                System.err.println(path + "\t" + e.getLocalizedMessage());
            }
        }
    }
}
