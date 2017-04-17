package edu.msu.cme.rdp.initprocess;

import edu.msu.cme.pyro.PipelineGene;
import java.io.File;
import java.io.IOException;

import edu.msu.cme.rdp.readseq.SequenceParsingException;
import edu.msu.cme.rdp.readseq.readers.MultiFileSeqReader;
import edu.msu.cme.rdp.readseq.utils.BarcodeSorter;
import edu.msu.cme.rdp.readseq.utils.BarcodeUtils.BarcodeInvalidException;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import edu.msu.cme.rdp.seqfilter.SeqFilterChain;
import edu.msu.cme.rdp.seqfilter.SeqFilteringResult;
import edu.msu.cme.rdp.seqfilter.filters.CharReplaceFilter;
import edu.msu.cme.rdp.seqfilter.filters.ExpQualityFilter;
import edu.msu.cme.rdp.seqfilter.filters.NSeqFilter;
import edu.msu.cme.rdp.seqfilter.filters.PrimerFilter;
import edu.msu.cme.rdp.seqfilter.filters.PrimerFilter.PrimerFilterListener;
import edu.msu.cme.rdp.seqfilter.filters.SeqLengthFilter;
import edu.msu.cme.rdp.seqfilter.filters.SeqLengthFilter.LengthFilterType;
import edu.msu.cme.rdp.seqfilter.filters.ValidAlphabetFilter;
import edu.msu.cme.rdp.seqfilter.filters.ValidAlphabetFilter.ValidAlphabetFilterMode;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class InitialProcessor {

    private static InitialProcessOptions cloneOptions(InitialProcessOptions o) {
        InitialProcessOptions ret = new InitialProcessOptions();

        ret.bestScoreOutfile = o.bestScoreOutfile;
        ret.fPrimer = o.fPrimer;
        ret.forwardMaxEditDist = o.forwardMaxEditDist;
        ret.genename = o.genename;
        ret.keepPrimers = o.keepPrimers;
        ret.minExpQualScore = o.minExpQualScore;
        ret.minSeqLength = o.minSeqLength;
        ret.noofns = o.noofns;
        ret.processNoTag = o.processNoTag;
        ret.qualInfile = o.qualInfile;
        ret.rPrimer = o.rPrimer;
        ret.reverseMaxEditDist = o.reverseMaxEditDist;
        ret.seqInfile = o.seqInfile;
        ret.trimQualOutfile = o.trimQualOutfile;
        ret.trimSeqOutfile = o.trimSeqOutfile;

        return ret;
    }

    private static Map<String, InitialProcessOptions> loadCustomOpts(File tagfile, InitialProcessOptions defaultOpts) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(tagfile));
        String line;
        Map<String, InitialProcessOptions> ret = new HashMap();

        while ((line = reader.readLine()) != null) {
            String[] lexemes = line.split("\\s+");

            if (lexemes.length == 3) {
                String sampleName = lexemes[1];
                String[] optsList = lexemes[2].split(",");

                InitialProcessOptions opts = cloneOptions(defaultOpts);
                for (String opt : optsList) {
                    if (!opt.contains("=")) {
                        System.err.println("Couldn't process opt " + opt);
                        continue;
                    }

                    String[] optLexemes = opt.split("=");
                    String name = optLexemes[0];
                    String val = "";
                    if (optLexemes.length > 1) {
                        val = optLexemes[1];
                    }

                    if (name.equals("fprimer")) {
                        opts.fPrimer = Arrays.asList(val.split("\\|"));
                    } else if (name.equals("rprimer")) {
                        opts.rPrimer = Arrays.asList(val.split("\\|"));
                    } else if (name.equals("fedit")) {
                        opts.forwardMaxEditDist = Integer.valueOf(val);
                    } else if (name.equals("redit")) {
                        opts.reverseMaxEditDist = Integer.valueOf(val);
                    } else if (name.equals("gene")) {
                        if("RRNA16S".equals(val)) {
                            opts.genename = PipelineGene.RRNA_16S_BACTERIA;
                        } else if("RRNA28S".equals(val)) {
                            opts.genename = PipelineGene.RRNA_28S;
                        } else {
                            opts.genename = PipelineGene.valueOf(val);
                        }
                    } else if (name.equals("keep_primer")) {
                        opts.keepPrimers = Boolean.valueOf(val);
                    } else if (name.equals("min_qual")) {
                        opts.minExpQualScore = Integer.valueOf(val);
                    } else if (name.equals("min_length")) {
                        opts.minSeqLength = Integer.valueOf(val);
                    } else if (name.equals("max_ns")) {
                        opts.noofns = Integer.valueOf(val);
                    } else {
                        System.err.println("Unprocessed option " + opt);
                    }
                }

                ret.put(sampleName, opts);
            }
        }

        return ret;
    }

    private static void printParams(InitialProcessOptions params, File f) throws IOException {
        PrintStream out = new PrintStream(f);

        out.print("Input file: ");
        out.print(params.seqInfile.get(0).getName());
        for (int index = 1; index < params.seqInfile.size(); index++) {
            out.print(", ");
            out.print(params.seqInfile.get(index).getName());
        }

        if (params.qualInfile == null || !params.qualInfile.exists()) {
            out.println("Quality file: (None)");
        } else {
            out.println("Quality file: " + params.qualInfile.getName());
        }

        out.println("Gene name: " + params.genename.toString());
        out.println("Maximum number of Ns: " + params.noofns);
        out.println("Minimum sequence length: " + params.minSeqLength);
        out.println("Minimum exponential Q-score: " + params.minExpQualScore);
        out.println("Primers removed?: " + !params.keepPrimers);
        out.println("Process no tag?: " + params.processNoTag);
       //out.println("Forward primer(s): " + Arrays.asList(params.fPrimer).toString().replace("[", "").replace("]", ""));
        //out.println("Max forward primer distance: " + params.forwardMaxEditDist);

        if (params.fPrimer == null || params.fPrimer.isEmpty()) {
            out.println("Forward primer(s): (None)");
        } else {
            out.println("Forward primer(s): " + Arrays.asList(params.fPrimer).toString().replace("[", "").replace("]", ""));
            out.println("Max forward primer distance: " + params.forwardMaxEditDist);
        }

        if (params.rPrimer == null || params.rPrimer.isEmpty()) {
            out.println("Reverse primer(s): (None)");
        } else {
            out.println("Reverse primer(s): " + Arrays.asList(params.rPrimer).toString().replace("[", "").replace("]", ""));
            out.println("Max reverse primer distance: " + params.reverseMaxEditDist);
        }

        out.close();
    }

    public static void runSeqTrim(Map<String, List<File>> tagSortedFiles, File resultDir, File tagSortDir, InitialProcessOptions defaultOptions, Map<String, InitialProcessOptions> customOptsMap) throws IOException, SequenceParsingException {
        for (String tagName : tagSortedFiles.keySet()) {
            if (tagName.equals(BarcodeSorter.NoTag) && !defaultOptions.processNoTag) {
                continue;
            }

            File tagDir = new File(resultDir, tagName);
            if (!tagDir.mkdir()) {
                throw new IOException("Failed to make initial process directory " + tagDir.getAbsolutePath());
            }

            InitialProcessOptions options = defaultOptions;

            if (customOptsMap.containsKey(tagName)) {
                options = customOptsMap.get(tagName);
                printParams(options, new File(tagDir, "custom_options.txt"));
            }


            File trimSeqFastaOutfile = new File(tagDir, tagName + "_trimmed.fasta");
            File trimSeqFastqOutfile = new File(tagDir, tagName + "_trimmed.fastq");
            File bestScoreOutfile = new File(tagDir, tagName + "_bestscore.txt");
            File trimQualOutfile = new File(tagDir, tagName + "_trimmed.qual");
            File lengthHistogramFile = new File(tagDir, tagName + "_length_histo.png");
            File qualityChartFile = new File(tagDir, tagName + "_quality.png");
            File summaryFile = new File(tagDir, tagName + "_summary.txt");
            File lengthStatsFile = new File(tagDir, tagName + "_length_stats.txt");
            File qualStatsFile = new File(tagDir, tagName + "_qual_stats.txt");
            File fileredSeqsFile = new File(tagDir, tagName + "_dropped_seqs.txt");

            final PrintStream bestScoresOut = new PrintStream(bestScoreOutfile);
            PrimerFilterListener l = new PrimerFilterListener() {
                public void sequencePassed(String seqName, int seqStart, int seqStop, int forwardPrimer, int reversePrimer, int forwardScore, int reverseScore, int partialScore, int ns, int length, boolean reversed) {
                    bestScoresOut.println(seqName
                            + "|fPrimer|" + forwardPrimer
                            + "|fPrimer.bestScore|" + forwardScore
                            + "|rPrimer|" + reversePrimer
                            + "|rPrimer.bestScore|" + reverseScore
                            + "|dpscore|" + partialScore
                            + "|noofns|" + ns
                            + "|trimmed_len|" + length
                            + "|reverse|" + reversed);
                }
            };

            SeqFilterChain chain = new SeqFilterChain(
                    new CharReplaceFilter(CharReplaceFilter.dnaReplaceMap),
                    new ValidAlphabetFilter(ValidAlphabetFilterMode.DROP_SEQUENCE, SeqUtils.RNAAlphabet),
                    new PrimerFilter(options.fPrimer, options.rPrimer, options.forwardMaxEditDist, options.reverseMaxEditDist, options.genename != PipelineGene.OTHER, options.keepPrimers, l),
                    new NSeqFilter(options.noofns),
                    new SeqLengthFilter(options.minSeqLength, LengthFilterType.GreaterThan),
                    new SeqLengthFilter(options.maxSeqLength, LengthFilterType.LessThan),
                    new ExpQualityFilter(options.minExpQualScore) //new CaseTransformFilter()
                    );

            InitProcessOutput out = new InitProcessOutput(trimSeqFastqOutfile, trimSeqFastaOutfile);
            SeqFilteringResult results = chain.filterSeqs(new MultiFileSeqReader(tagSortedFiles.get(tagName)), out);
            out.writeStats(fileredSeqsFile, summaryFile, lengthStatsFile, lengthHistogramFile, qualStatsFile, qualityChartFile, tagName, results);
            out.close();
            bestScoresOut.close();

            if(trimSeqFastaOutfile.length() == 0) {
                trimSeqFastaOutfile.delete();
            }
            if(trimSeqFastqOutfile.length() == 0) {
                trimSeqFastqOutfile.delete();
            }
        }

        printParams(defaultOptions, new File(resultDir, "input_params.txt"));
    }

    public static void doInitialProcessing(File resultDir, File tagSortDir, File tagFile, InitialProcessOptions defaultOptions) throws IOException, SequenceParsingException, BarcodeInvalidException {

        if (!tagSortDir.exists() && !tagSortDir.mkdir()) {
            throw new IOException("Failed to make tagsort dir " + tagSortDir.getAbsolutePath());
        }

        if (!resultDir.exists() && !resultDir.mkdir()) {
            throw new IOException("Failed to make result dir " + resultDir.getAbsolutePath());
        }

        List<File> inputFiles = defaultOptions.seqInfile;
        Map<String, List<File>> tagSortedFiles;
        Map<String, InitialProcessOptions> customOptsMap;

        if (tagFile != null) {
            customOptsMap = loadCustomOpts(tagFile, defaultOptions);
            tagSortedFiles = BarcodeSorter.sort(inputFiles, tagFile, tagSortDir, (byte) defaultOptions.minExpQualScore);
        } else {
            customOptsMap = new HashMap();
            tagSortedFiles = new HashMap();
            tagSortedFiles.put(BarcodeSorter.NoTag, inputFiles);
        }

        // create
        runSeqTrim(tagSortedFiles, resultDir, tagSortDir, defaultOptions, customOptsMap);
    }
}
