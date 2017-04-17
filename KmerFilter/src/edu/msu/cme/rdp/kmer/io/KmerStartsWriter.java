/*
 * Copyright (C) 2013 Jordan Fish <fishjord at msu.edu>
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
package edu.msu.cme.rdp.kmer.io;

import java.io.*;

/**
 *
 * @author Jordan Fish <fishjord at msu.edu>
 */
public class KmerStartsWriter {

    private PrintStream out;

    public KmerStartsWriter(OutputStream os) throws IOException {
        out = new PrintStream(new BufferedOutputStream(os));
        out.println("#gene name\tquery id\trefid\tnucl kmer\tis prot?\tstarting_frame\tprot kmer\tmodel pos");
    }

    public KmerStartsWriter(File f) throws IOException {
        this(new FileOutputStream(f));
    }

    public KmerStartsWriter(String file) throws IOException {
        this(new File(file));
    }

    public void write(KmerStart start) throws IOException {
        out.print(start.getGeneName() + "\t"
                + start.getQueryId() + "\t"
                + start.getRefId() + "\t"
                + start.getNuclKmer() + "\t"
                + start.isProt() + "\t");

        if(start.isProt()) {
            out.print(start.getFrame() + "\t" + start.getProtKmer() + "\t");
        } else {
            out.print("-\t-\t");
        }

        out.println(start.getMpos());
    }

    public void close() {
        out.close();
    }
}
