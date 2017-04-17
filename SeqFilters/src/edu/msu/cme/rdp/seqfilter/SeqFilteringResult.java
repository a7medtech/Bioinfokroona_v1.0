/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.msu.cme.rdp.seqfilter;

import java.util.Map;

/**
 *
 * @author fishjord
 */
public class SeqFilteringResult {
    private Map<SeqFilter, Map<String, SeqFilterResult>> filteredSeqsMap;
    private int totalSeqs;
    private int filteredSeqs;

    public SeqFilteringResult(int totalSeqs, int filteredSeqs, Map<SeqFilter, Map<String, SeqFilterResult>> filteredSeqsMap) {
        this.filteredSeqsMap = filteredSeqsMap;
        this.totalSeqs = totalSeqs;
        this.filteredSeqs = filteredSeqs;
    }

    public int getFilteredSeqs() {
        return filteredSeqs;
    }

    public Map<SeqFilter, Map<String, SeqFilterResult>> getFilteredSeqsMap() {
        return filteredSeqsMap;
    }

    public int getTotalSeqs() {
        return totalSeqs;
    }
}
