/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.seqfilter.filters;

import edu.msu.cme.rdp.probematch.myers99.PatternBitMask64;
import edu.msu.cme.rdp.readseq.QSequence;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.utils.IUBUtilities;
import edu.msu.cme.rdp.readseq.utils.orientation.OrientationChecker;
import edu.msu.cme.rdp.seqfilter.SeqFilter;
import edu.msu.cme.rdp.seqfilter.SeqFilterResult;
import edu.msu.cme.rdp.seqfilter.utils.TrimHelper;
import edu.msu.cme.rdp.seqfilter.utils.TrimHelper.TrimResult;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author fishjord
 */
public class PrimerFilter implements SeqFilter {

    public static interface PrimerFilterListener {

        public void sequencePassed(String seqName, int seqStart, int seqStop, int forwardPrimer, int reversePrimer, int forwardScore, int reverseScore, int partialScore, int ns, int length, boolean reversed);
    }
    private static int MIN_DPSCORE = 6;
    private PatternBitMask64[] forwardPrimers;
    private PatternBitMask64[] reversePrimers;
    private String[] forwardPrimerStrs;
    private String[] reversePrimerStrs;
    private int maxForwardDist = 0;
    private int maxReverseDist = 0;
    private boolean reversible = true;
    private boolean keepPrimers = false;
    private PrimerFilterListener listener;

    public PrimerFilter(List<String> forwardPrimerStrs, List<String> reversePrimerStrs, int maxForwardDist, int maxReverseDist, boolean reversible, boolean keepPrimers) {
        this(forwardPrimerStrs, reversePrimerStrs, maxForwardDist, maxReverseDist, reversible, keepPrimers, null);
    }

    public PrimerFilter(List<String> forwardPrimerStrs, List<String> reversePrimerStrs, int maxForwardDist, int maxReverseDist, boolean reversible, boolean keepPrimers, PrimerFilterListener listener) {
        this.maxForwardDist = maxForwardDist;
        this.maxReverseDist = maxReverseDist;

        this.keepPrimers = keepPrimers;
        this.reversible = reversible;

        if (!forwardPrimerStrs.isEmpty()) {
            forwardPrimers = new PatternBitMask64[forwardPrimerStrs.size()];
            this.forwardPrimerStrs = forwardPrimerStrs.toArray(new String[forwardPrimerStrs.size()]);
            for (int primer = 0; primer < forwardPrimerStrs.size(); primer++) {
                forwardPrimers[primer] = new PatternBitMask64(forwardPrimerStrs.get(primer), true);
            }
        }

        if (!reversePrimerStrs.isEmpty()) {
            this.reversePrimerStrs = reversePrimerStrs.toArray(new String[reversePrimerStrs.size()]);

            this.reversePrimers = new PatternBitMask64[reversePrimerStrs.size()];
            for (int primer = 0; primer < reversePrimerStrs.size(); primer++) {
                reversePrimers[primer] = new PatternBitMask64(reversePrimerStrs.get(primer), true);
            }
        }

        this.listener = listener;
    }

    public PrimerFilter(String[] forwardPrimerStrs, int maxForwardDist, boolean reversible) {
        this.maxForwardDist = maxForwardDist;

        this.reversible = reversible;

        forwardPrimers = new PatternBitMask64[forwardPrimerStrs.length];
        this.forwardPrimerStrs = Arrays.copyOf(forwardPrimerStrs, forwardPrimerStrs.length);
        for (int primer = 0; primer < forwardPrimerStrs.length; primer++) {
            if (forwardPrimerStrs[primer].equals("")) {
                throw new IllegalArgumentException("Foward primer " + primer + " is empty");
            }

            forwardPrimers[primer] = new PatternBitMask64(forwardPrimerStrs[primer], true);
        }
    }

    public String getName() {
        return "Primer trimmer";
    }

    public SeqFilterResult filterSequence(Sequence seq) {
        if (seq instanceof QSequence) {
            return doFilter(seq.getSeqName(), seq.getDesc(), seq.getSeqString(), ((QSequence) seq).getQuality());
        } else {
            return doFilter(seq.getSeqName(), seq.getDesc(), seq.getSeqString(), null);
        }
    }

    public SeqFilterResult doFilter(String seqName, String desc, String origSeqString, byte[] qualSeq) {
        String newSeqString = origSeqString;
        int forwardPrimer = -1;
        int forwardScore = -1;
        int reversePrimer = -1;
        int reverseScore = -1;
        int partialMatch = -1;
        int forwardStop = 0;
        boolean reversed = false;

        if (forwardPrimers != null) {
            TrimResult forwardTrimming = TrimHelper.trimSequence(newSeqString, forwardPrimers);

            if (forwardTrimming.getEditDistance() > maxForwardDist) {
                return new SeqFilterResult("Forward primer didn't hit, best match is with primer " + forwardTrimming.getPrimer() + " with a distance of " + forwardTrimming.getEditDistance());
            }

            forwardPrimer = forwardTrimming.getPrimer();
            forwardScore = forwardTrimming.getEditDistance();
            forwardStop = forwardTrimming.getPrimerStopIndex();

            if (keepPrimers) {
                forwardStop = forwardTrimming.getPrimerStopIndex() - forwardPrimerStrs[forwardTrimming.getPrimer() - 1].length();
                if (forwardStop > 0) {
                    newSeqString = newSeqString.substring(forwardStop);
                } else {
                    forwardStop = 0;
                }
            } else {
                newSeqString = forwardTrimming.getTrimmedSeq();
            }
        }

        if (reversePrimers != null) {
            String reversePrimerSeq = IUBUtilities.reverseComplement(newSeqString);
            TrimResult reverseTrimming = TrimHelper.trimSequence(reversePrimerSeq, reversePrimers);

            if (reverseTrimming.getEditDistance() > maxReverseDist) {
                TrimResult partialTrimming = TrimHelper.partialMatch(reversePrimerSeq, reversePrimerStrs);

                if (partialTrimming.getEditDistance() < MIN_DPSCORE) {
                    return new SeqFilterResult("Reverse primer didn't hit"
                            + ((partialTrimming.getEditDistance() == Integer.MIN_VALUE) ? "" : ", best partial match is with primer " + partialTrimming.getPrimer() + " with a distance of " + partialTrimming.getEditDistance()));
                } else {
                    if (maxReverseDist == 0 && partialTrimming.getEditDistance() != partialTrimming.getPrimerStopIndex()) {
                        return new SeqFilterResult("Reverse primer didn't hit"
                                + ((partialTrimming.getEditDistance() == Integer.MIN_VALUE) ? "" : ", best partial match is with primer " + partialTrimming.getPrimer() + " with a distance of " + partialTrimming.getEditDistance()));
                    }
                }

                reversePrimer = partialTrimming.getPrimer();
                partialMatch = partialTrimming.getEditDistance();

                if (keepPrimers) {
                    newSeqString = reversePrimerSeq;
                } else {
                    newSeqString = partialTrimming.getTrimmedSeq();
                }
            } else {
                reversePrimer = reverseTrimming.getPrimer();
                reverseScore = reverseTrimming.getEditDistance();

                if (keepPrimers) {
                    if (reverseTrimming.getPrimerStopIndex() - reversePrimerStrs[reverseTrimming.getPrimer() - 1].length() > 0) {
                        newSeqString = reversePrimerSeq.substring(reverseTrimming.getPrimerStopIndex() - reversePrimerStrs[reverseTrimming.getPrimer() - 1].length());
                    } else {
                        newSeqString = reversePrimerSeq;
                    }
                } else {
                    newSeqString = reverseTrimming.getTrimmedSeq();
                }
            }


            newSeqString = IUBUtilities.reverseComplement(newSeqString);
        }

        if (reversible && OrientationChecker.getChecker().isSeqReversed(origSeqString)) {
            newSeqString = IUBUtilities.reverseComplement(newSeqString);
            reversed = true;
        }

        byte[] newQualSeq = null;
        if (qualSeq != null && qualSeq.length >= forwardStop + newSeqString.length()) {
            newQualSeq = Arrays.copyOfRange(qualSeq, forwardStop, forwardStop + newSeqString.length());

            if (reversed) {
                int left = 0;
                int right = newQualSeq.length - 1;

                while (left < right) {
                    byte tmp = newQualSeq[left];
                    newQualSeq[left] = newQualSeq[right];
                    newQualSeq[right] = tmp;

                    left++;
                    right--;
                }
            }
        }

        if (listener != null) {
            listener.sequencePassed(seqName,
                    forwardStop,
                    forwardStop + newSeqString.length(),
                    forwardPrimer,
                    reversePrimer,
                    forwardScore,
                    reverseScore,
                    partialMatch,
                    newSeqString.replaceAll("[nN]", "").length(),
                    newSeqString.length(),
                    reversed);
        }

        if (newQualSeq != null) {
            return new SeqFilterResult(new QSequence(seqName, "length=" + newSeqString.length(), newSeqString, newQualSeq));
        }

        return new SeqFilterResult(new Sequence(seqName, "length=" + newSeqString.length(), newSeqString));
    }
}
