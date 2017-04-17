/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.seqfilter.utils;

import edu.msu.cme.rdp.probematch.myers99.BitVector64;
import edu.msu.cme.rdp.probematch.myers99.BitVector64Result.BitVector64Match;
import edu.msu.cme.rdp.probematch.myers99.PatternBitMask64;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;

/**
 *
 * @author fishjord
 */
public class TrimHelper {

    private static final int MATCH_SCORE = 1;
    private static final int MISMATCH_SCORE = -6;
    private static final int GAP_SCORE = -6;

    public static class TrimResult {

        private String trimmedSeq;
        private int primer;
        private int editDistance;
        private int primerStopIndex;

        public TrimResult(String trimmedSeq, int primer, int editDistance, int primerStopIndex) {
            this.trimmedSeq = trimmedSeq;
            this.primer = primer;
            this.editDistance = editDistance;
            this.primerStopIndex = primerStopIndex;
        }

        public int getEditDistance() {
            return editDistance;
        }

        public int getPrimer() {
            return primer;
        }

        public int getPrimerStopIndex() {
            return primerStopIndex;
        }

        public String getTrimmedSeq() {
            return trimmedSeq;
        }
    }

    private static int[] partialMatch(String seqString, String primer) {
        int[][] dpmatrix = new int[primer.length() + 1][seqString.length() + 1];
        char[] primerChars = primer.toCharArray();
        char[] seqChars = seqString.toCharArray();

        for (int i = 0; i <= primer.length(); ++i) {
            dpmatrix[i][0] = 0;
        }
        for (int i = 1; i <= seqString.length(); ++i) {
            dpmatrix[0][i] = Integer.MIN_VALUE / 2;
        }

        for (int i = 1; i <= primer.length(); ++i) {
            for (int j = 1; j <= seqString.length(); ++j) {
                int v1 = SeqUtils.IUPAC[primerChars[i - 1]];
                int v2 = SeqUtils.IUPAC[seqChars[j - 1]];

                if (v1 > v2) {
                    int tmp = v1;
                    v1 = v2;
                    v2 = tmp;
                }

                boolean match = ((v1 & v2) == v1);

                int down = dpmatrix[i - 1][j] + GAP_SCORE;
                int across = dpmatrix[i][j - 1] + GAP_SCORE;
                int diag = dpmatrix[i - 1][j - 1]
                        + ((match) ? MATCH_SCORE : MISMATCH_SCORE);
                dpmatrix[i][j] = Math.max(diag, Math.max(down, across));
            }
        }

        int bestPos = -1;
        int bestScore = 0;
        for (int position = 0; position <= seqString.length(); ++position) {
            int currentScore = dpmatrix[primer.length()][position];
            if (currentScore > bestScore) {
                bestScore = currentScore;
                bestPos = position;
            }
        }

        return new int[]{bestPos, bestScore};
    }

    public static TrimResult partialMatch(String seqString, String[] primers) {
        int bestScore = Integer.MIN_VALUE;
        int bestPos = 0;
        int bestPrimer = 0;

        for (int i = 0; i < primers.length; i++) {
            if (seqString.length() > primers[i].length()) {

                int[] tuple = partialMatch(seqString, primers[i]);
                int pos = tuple[0];
                int score = tuple[1];

                if (pos > -1 && score > bestScore) {
                    bestScore = score;
                    bestPos = pos;
                    bestPrimer = i + 1;
                }
            }
        }

        return new TrimResult(seqString.substring(bestPos), bestPrimer, bestScore, bestPos);
    }

    public static TrimResult trimSequence(String seqString, PatternBitMask64[] primers) {
        int bestScore = Integer.MAX_VALUE;
        int bestPrimer = 0;
        int bestPos = 0;

        for (int primerIndex = 0; primerIndex < primers.length; primerIndex++) {
            BitVector64Match result = BitVector64.process(seqString, primers[primerIndex]).getBestResult();

            if (result.getScore() < bestScore) {
                bestScore = result.getScore();
                bestPos = result.getPosition();
                bestPrimer = primerIndex + 1;
            }
        }

        return new TrimResult(seqString.substring(bestPos), bestPrimer, bestScore, bestPos);
    }
}
