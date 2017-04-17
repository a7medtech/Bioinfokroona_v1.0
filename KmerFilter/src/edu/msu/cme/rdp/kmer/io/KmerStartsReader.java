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
public class KmerStartsReader {

    private BufferedReader reader;

    public KmerStartsReader(InputStream is) throws IOException {
        this.reader = new BufferedReader(new InputStreamReader(is));
    }

    public KmerStartsReader(File f) throws IOException {
        this(new FileInputStream(f));
    }

    public KmerStartsReader(String file) throws IOException {
        this(new File(file));
    }

    public KmerStart readNext() throws IOException {
        String line;
        do {
            line = reader.readLine();
        } while (line != null && (line.isEmpty() || line.charAt(0) == '#'));

        if(line == null) {
            return null;
        }

        String[] lexemes = line.split("\t");

        if (lexemes.length != 8) {
            throw new IOException("Malformed line '" + line + "'");
        }

        boolean isprot = Boolean.valueOf(lexemes[4]);
        String protKmer;
        int frame;
        if (isprot) {
            frame = Integer.valueOf(lexemes[5]);
            protKmer = lexemes[6];
        } else {
            frame = 0;
            protKmer = null;
        }

        return new KmerStart(lexemes[0], lexemes[1], lexemes[2], lexemes[3], frame, Integer.valueOf(lexemes[7]), true, protKmer);
    }

    public void close() throws IOException {
        reader.close();
    }
}
