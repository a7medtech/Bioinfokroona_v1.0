/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.seqfilter;

import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.SequenceParsingException;
import edu.msu.cme.rdp.readseq.readers.QSeqReader;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.seqfilter.output.SeqFilterOutput;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author fishjord
 */
public class SeqFilterChain {

    private List<SeqFilter> filters;

    public SeqFilterChain(SeqFilter... filters) {
        this.filters = Arrays.asList(filters);
    }

    public SeqFilteringResult filterSeqs(SeqReader reader, SeqFilterOutput out) throws IOException, SequenceParsingException {

        Sequence seq;
        Map<SeqFilter, Map<String, SeqFilterResult>> filteredSeqsMap = new LinkedHashMap();
        int totalSeqs = 0;
        int filteredSeqs = 0;

        for (SeqFilter filter : filters) {
            filteredSeqsMap.put(filter, new HashMap());
        }

        while ((seq = reader.readNextSequence()) != null) {
            SeqFilterResult result = runFilterChain(seq, filteredSeqsMap);
            if (result.failed()) {
                filteredSeqs++;
            } else {
                out.appendSequence(result.getResultSeq());

            }
            totalSeqs++;
        }

        return new SeqFilteringResult(totalSeqs, filteredSeqs, filteredSeqsMap);
    }

    public SeqFilteringResult filterSeqs(File inFile, SeqFilterOutput out) throws IOException, SequenceParsingException {
        SequenceReader seqReader = new SequenceReader(inFile);
        try {
            return filterSeqs(seqReader, out);
        } finally {
            seqReader.close();
        }
    }

    public SeqFilteringResult filterSeqs(File inFile, File qualityFile, SeqFilterOutput out) throws IOException, SequenceParsingException {
        QSeqReader seqReader = new QSeqReader(inFile, qualityFile);

        Sequence seq;
        Map<SeqFilter, Map<String, SeqFilterResult>> filteredSeqsMap = new LinkedHashMap();
        int totalSeqs = 0;
        int filteredSeqs = 0;

        for (SeqFilter filter : filters) {
            filteredSeqsMap.put(filter, new HashMap());
        }

        while ((seq = seqReader.readNextSequence()) != null) {
            SeqFilterResult result = runFilterChain(seq, filteredSeqsMap);
            if (result.failed()) {
                filteredSeqs++;
            } else {
                out.appendSequence(result.getResultSeq());

            }
            totalSeqs++;
        }

        seqReader.close();

        return new SeqFilteringResult(totalSeqs, filteredSeqs, filteredSeqsMap);
    }

    private SeqFilterResult runFilterChain(Sequence seq, Map<SeqFilter, Map<String, SeqFilterResult>> filteredSeqsMap) {

        SeqFilterResult result = new SeqFilterResult("No filters specified");

        for (SeqFilter filter : filters) {
            result = filter.filterSequence(seq);

            if (result.failed()) {
                filteredSeqsMap.get(filter).put(seq.getSeqName(), result);
                break;
            } else {
                seq = result.getResultSeq();
            }
        }

        return result;
    }
}
