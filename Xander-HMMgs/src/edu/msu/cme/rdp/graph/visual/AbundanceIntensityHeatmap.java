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
package edu.msu.cme.rdp.graph.visual;

import edu.msu.cme.rdp.readseq.readers.IndexedSeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import javax.imageio.ImageIO;

/**
 *
 * @author fishjord
 */
public class AbundanceIntensityHeatmap {

    private static Color A = Color.green;
    private static Color C = Color.blue;
    private static Color G = Color.yellow;
    private static Color T = Color.red;
    private static final int a = 0, c = 1, g = 2, t = 3;

    public static class AbundanceSeq {

        private char[] seq;
        private int[] occurences;
    }
    private final int maxOccur;
    private final int alnCols;
    private final int[][] baseOccur;
    private final Color[][] baseColorRange;
    private final int numSeqs;
    private final byte[] mask;

    public AbundanceIntensityHeatmap(List<AbundanceSeq> seqs, int alnCols) {
        this.alnCols = alnCols;
        int maxAbund = 0;
        numSeqs = seqs.size();
        mask = new byte[alnCols];

        baseOccur = new int[4][alnCols];

        for (int col = 0; col < alnCols; col++) {
            for (AbundanceSeq seq : seqs) {
                if (seq.seq[col] == '-' || Character.isUpperCase(seq.seq[col])) {
                    mask[col] = 0;
                } else {
                    mask[col] = 1;
                }

                switch (seq.seq[col]) {
                    case 'A':
                    case 'a':
                        baseOccur[a][col] = Math.max(baseOccur[a][col], seq.occurences[col]);
                        break;
                    case 'C':
                    case 'c':
                        baseOccur[c][col] = Math.max(baseOccur[c][col], seq.occurences[col]);
                        break;
                    case 'G':
                    case 'g':
                        baseOccur[g][col] = Math.max(baseOccur[g][col], seq.occurences[col]);
                        break;
                    case 'T':
                    case 't':
                    case 'U':
                    case 'u':
                        baseOccur[t][col] = Math.max(baseOccur[t][col], seq.occurences[col]);
                        break;
                }
            }

            maxAbund = Math.max(maxAbund, Math.max(baseOccur[a][col], Math.max(baseOccur[c][col], Math.max(baseOccur[g][col], baseOccur[t][col]))));
        }

        System.out.println(maxAbund);

        this.maxOccur = maxAbund;
        baseColorRange = new Color[4][maxOccur + 1];

        initColors();
    }

    private void initColors() {

        System.out.println(maxOccur);

        float step = 255.0f / maxOccur;

        int lastVal = -1;
        for (int index = 0; index < maxOccur + 1; index++) {
            int val = (int) (index * step);

            if (val > lastVal) {
                baseColorRange[a][index] = new Color(255 - val, 255, 255 - val);
                baseColorRange[g][index] = new Color(255 - val, 255 - val, 255);
                baseColorRange[c][index] = new Color(255, 255 - val, 255);
                baseColorRange[t][index] = new Color(255, 255 - val, 255 - val);
                lastVal = val;
            } else {
                baseColorRange[a][index] = baseColorRange[a][index - 1];
                baseColorRange[c][index] = baseColorRange[c][index - 1];
                baseColorRange[g][index] = baseColorRange[g][index - 1];
                baseColorRange[t][index] = baseColorRange[t][index - 1];
            }
        }
    }

    public BufferedImage drawImage(int width, int height, int startCol, int endCol) {
        int drawCols = endCol - startCol;

        if (startCol < 0 || endCol > alnCols || endCol < startCol) {
            throw new IllegalArgumentException("Invalid starting and/or ending columns " + startCol + "/" + endCol);
        }

        if (height < 20 || width < drawCols * 5) {
            throw new IllegalArgumentException("Can't fit " + drawCols + " alignment columns in to a" + width + " by " + height + " image");
        }

        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        Graphics2D g = image.createGraphics();
        g.setColor(Color.white);
        g.fillRect(0, 0, width, height);

        int cellWidth = ((width - drawCols) / drawCols); //Black lines around each cell, and the borders
        int cellHeight = (height / 4);  //Black lines around each cell and the border

        int y = 0;
        for (int base = 0; base <= t; base++) {
            Color[] colorRange = baseColorRange[base];
            int[] occur = baseOccur[base];

            int x = 0;
            for (int col = startCol; col < endCol; col++) {
                g.setColor(colorRange[occur[col]]);
                g.fillRect(x, y, cellWidth, cellHeight);

                x += cellWidth + 1;
            }
            y += cellHeight + 1;
        }

        return image;
    }

    public BufferedImage getColorKey(int width, int height) {
        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        Graphics2D g = image.createGraphics();

        float segmentWidth = (float)width / maxOccur;
        float segmentHeight = (float)height / 4;
        if (segmentWidth < 1) {
            throw new IllegalArgumentException("Can't fit " + maxOccur + " segments in to " + height + " color key");
        }

        float x = 0;
        for (int index = 0; index <= maxOccur; index++) {
            float y = 0;
            for (int base = a; base <= t; base++) {
                g.setColor(baseColorRange[base][index]);
                g.fillRect((int)x, (int)y, (int)segmentWidth, (int)segmentHeight);

                y += segmentHeight;
            }

            g.setColor(Color.black);
            g.drawString(index + "", (int)x + 1, 50);
            x += segmentWidth;
        }

        g.setColor(Color.black);
        x = 0;
        for (int index = 0; index <= maxOccur; index++) {
            g.drawLine((int)x, 0, (int)x, height);
            x += segmentWidth;
        }

        return image;
    }

    public String getCons(int startCol, int endCol) {

        if (startCol < 0 || endCol > alnCols || endCol < startCol) {
            throw new IllegalArgumentException("Invalid starting and/or ending columns " + startCol + "/" + endCol);
        }
        StringBuilder cons = new StringBuilder();

        for (int col = startCol; col < endCol; col++) {
            int appendId = 0;
            int maxOccur = 0;
            int totalOccur = 0;

            for (int base = a; base <= t; base++) {
                if (baseOccur[base][col] > maxOccur) {
                    appendId = base;
                    maxOccur = baseOccur[base][col];
                }

                totalOccur += baseOccur[base][col];
            }

            char append;

            if (appendId == a) {
                append = 'a';
            } else if (appendId == c) {
                append = 'c';
            } else if (appendId == g) {
                append = 'g';
            } else if (appendId == t) {
                append = 't';
            } else {
                throw new IllegalArgumentException("I hate you.");
            }

            float perc = (float) maxOccur / totalOccur;
            if (perc < .25) {
                if (mask[col] == 1) {
                    cons.append(".");
                } else {
                    cons.append("-");
                }
            } else if (perc > .75) {
                append = Character.toUpperCase(append);
            }

            cons.append(append);
        }

        return cons.toString();
    }

    public static void main(String[] args) throws IOException {
        if (args.length != 3 && args.length != 6) {
            System.err.println("USAGE: IntensityHeatmap <occurence_file> <alignment_file> <output_stem> [cols_per_img img_width img_height]");
            System.exit(1);
        }

        File occurenceFile = new File(args[0]);
        File alignmentFile = new File(args[1]);
        File resultStem = new File(args[2]);

        int alnCols = 0;
        int colsPerImage = 200;
        int width = 1000;
        int height = 120;

        if (args.length > 3) {
            colsPerImage = Integer.valueOf(args[3]);
            width = Integer.valueOf(args[4]);
            height = Integer.valueOf(args[5]);
        }

        System.err.println("Generating Intensity Heatmap");
        System.err.println("*  Input file:         " + alignmentFile);
        System.err.println("*  Occurence file:     " + occurenceFile);
        System.err.println("*  Result stem:        " + resultStem);
        System.err.println("*  Columns per image:  " + colsPerImage);
        System.err.println("*  Per image width:    " + width);
        System.err.println("*  Per image height:   " + height);

        BufferedReader reader = new BufferedReader(new FileReader(occurenceFile));
        String line;

        IndexedSeqReader seqReader = new IndexedSeqReader(alignmentFile);
        List<AbundanceSeq> seqs = new ArrayList();

        while ((line = reader.readLine()) != null) {
            line = line.trim();
            if ("".equals(line)) {
                continue;
            }

            String lexemes[] = line.split("\\s+");
            String seqid = lexemes[0];
            int[] occur = new int[lexemes.length - 1];
            for (int index = 1; index < lexemes.length; index++) {
                occur[index - 1] = Integer.valueOf(lexemes[index]);
            }

            Sequence seq;
            try {
                seq = seqReader.readSeq(seqid);

                char[] seqString = seq.getSeqString().toCharArray();
                if (alnCols == 0) {
                    alnCols = seqString.length;
                } else if (seqString.length < alnCols) {
                    throw new IOException("Not enough columns in sequence " + seqid);
                }

                AbundanceSeq abundSeq = new AbundanceSeq();
                abundSeq.seq = seqString;
                abundSeq.occurences = new int[seqString.length];
                int occurIndex = 0;
                for (int index = 0; index < alnCols; index++) {
                    if (seqString[index] == '.' || seqString[index] == '-') {
                        abundSeq.occurences[index] = 0;
                    } else {
                        abundSeq.occurences[index] = occur[occurIndex++];
                    }
                }

                seqs.add(abundSeq);

            } catch (IOException e) {
                System.err.println("Occurence information for " + seqid + " doesn't have a corrosponding seqid");
                System.exit(1);
            }
        }

        AbundanceIntensityHeatmap heatmap = new AbundanceIntensityHeatmap(seqs, alnCols);

        int index;
        for (index = 0; index <= alnCols - colsPerImage; index += colsPerImage) {
            BufferedImage image = heatmap.drawImage(width, height, index, index + colsPerImage);
            ImageIO.write(image, "png", new File(String.format(resultStem.getAbsolutePath() + "_%05d-%05d.png", (index + 1), (index + colsPerImage))));
        }

        if (index != alnCols - colsPerImage) {
            int remainingCols = alnCols - index;
            int tmpWidth = (int) (((float) remainingCols / colsPerImage) * width);

            BufferedImage image = heatmap.drawImage(tmpWidth, height, index, index + remainingCols);
            ImageIO.write(image, "png", new File(String.format(resultStem.getAbsolutePath() + "_%05d-%05d.png", (index + 1), (index + remainingCols))));
        }

        ImageIO.write(heatmap.getColorKey(heatmap.maxOccur * 30, 100), "png", new File(resultStem.getAbsolutePath() + "_key.png"));

        System.out.println(heatmap.getCons(0, alnCols));
        for (index = 0; index < heatmap.mask.length; index++) {
            System.out.print((heatmap.mask[index] == 0) ? 'x' : '.');
        }

        System.out.println();
    }
}
