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

import edu.msu.cme.pyro.cluster.dist.ThinEdge;
import edu.msu.cme.rdp.hadoop.distance.mapred.keys.DistanceAndComparison;
import edu.msu.cme.rdp.hadoop.utils.HDFSEdgeReader;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.util.GenericOptionsParser;

/**
 *
 * @author Jordan Fish <fishjord at msu.edu>
 */
public class ResortPartialResult {

    private static class ExpectedRange {

        final int minX, maxX;
        final int minY, maxY;
        final int taskid;

        IFile.Writer<DistanceAndComparison, NullWritable> out;

        public ExpectedRange(Configuration conf, FileSystem fs, String[] lexemes) throws IOException {
            taskid = Integer.valueOf(lexemes[0]);
            minX = Integer.valueOf(lexemes[1]);
            maxX = Integer.valueOf(lexemes[2]);
            minY = Integer.valueOf(lexemes[3]);
            maxY = Integer.valueOf(lexemes[4]);

            new File(lexemes[0]).mkdir();
            out = new IFile.Writer<DistanceAndComparison, NullWritable>(conf, fs, new Path(taskid + "/file.out"), DistanceAndComparison.class, NullWritable.class, null);
        }
    }

    private static List<ExpectedRange> loadRanges(Configuration conf, FileSystem fs, File f) throws IOException {
        List<ExpectedRange> ret = new ArrayList();

        BufferedReader reader = new BufferedReader(new FileReader(f));
        String line;
        while((line = reader.readLine()) != null) {
            String[] lexemes = line.split("\\s+");
            if(lexemes.length != 5) {
                continue;
            }

            ret.add(new ExpectedRange(conf, fs, lexemes));
        }

        return ret;
    }

    private static int compare(ExpectedRange range, DistanceAndComparison item) {
	if(item.first >= range.maxX) {
	    return 1;
	}

	if(item.first < range.minX) {
	    return -1;
	}

	if(item.second >= range.maxY) {
	    return 1;
	}

	if(item.second < range.minY) {
	    return -1;
	}

	return 0;
    }

    private static ExpectedRange binSearch(List<ExpectedRange> ranges, DistanceAndComparison item) {
	int start = 0;
	int end = ranges.size() - 1;
	int mid, c;
	ExpectedRange range;

	while(start <= end) {
	    mid = (start + end) / 2;

	    c = compare(ranges.get(mid), item);
	    //System.err.println("start = " + start + ", end = " + end + ", mid = " + mid + " minX = " + ranges.get(mid).minX +" maxX = " + ranges.get(mid).maxX +" minY = " + ranges.get(mid).minY +" maxY = " + ranges.get(mid).maxY + " item.first=" + item.first + " item.second = " + item.second + " c= " + c);

	    if(c == 0) {
		return ranges.get(mid);
	    } else if(c > 0) {
		start = mid + 1;
	    } else {
	        end = mid - 1;
	    }
	}

	return null;
    }

    public static void main(String[] args) throws Exception {
        //guessSplitRanges();
        Configuration conf = new Configuration();
        args = new GenericOptionsParser(conf, args).getRemainingArgs();

        if (args.length <2) {
            System.err.println("USAGE: ReadPartialResultFile <expected_splits> <matrix part...>");
            System.exit(1);
        }

        FileSystem fs = FileSystem.getLocal(conf);
        List<ExpectedRange> ranges = loadRanges(conf, fs, new File(args[0]));
        System.err.println("Loaded ranges");

        List<Path> matrixParts = new ArrayList();
        for(int index = 1;index < args.length;index++) {
            matrixParts.add(new Path(args[index]));
        }
        HDFSEdgeReader reader = new HDFSEdgeReader(conf, matrixParts);
        System.err.println("Created reader");
        DistanceAndComparison dist = new DistanceAndComparison();
        ThinEdge thinEdge;
        NullWritable nw = NullWritable.get();

	int edges = 0;
        while((thinEdge = reader.nextThinEdge()) != null) {
            dist.distance = thinEdge.getDist();
            dist.first = thinEdge.getSeqi();
            dist.second = thinEdge.getSeqj();
	    edges++;

	    ExpectedRange range = binSearch(ranges, dist);

            if(range == null) {
                System.err.println("Couldn't map " + dist.first + " " + dist.second);
            } else {
		range.out.append(dist, nw);
	    }

	    if(edges % 1000000 == 0) {
		System.err.println("Sorted " + edges + " edges");
	    }
        }

        for(ExpectedRange range : ranges) {
            range.out.close();
        }
    }
}
