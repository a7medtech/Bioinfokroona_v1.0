package edu.msu.cme.rdp.initprocess;

import edu.msu.cme.pyro.PipelineGene;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author siddiq15
 *
 */
public class InitialProcessOptions {

    public List<String> fPrimer = new ArrayList();
    public List<String> rPrimer = new ArrayList();  // assign default
    public int forwardMaxEditDist;
    public int reverseMaxEditDist;
    public int noofns;
    public int minSeqLength;
    public int maxSeqLength;
    public List<File> seqInfile;
    public File trimSeqOutfile;
    public File qualInfile;
    public File trimQualOutfile;
    public File bestScoreOutfile;
    public PipelineGene genename;
    public boolean keepPrimers = false;
    public int minExpQualScore;
    public boolean processNoTag = true;
}
