##RDP's Sequence Filtering pipeline and Sequence Initial Processor

### Intro
RDP Pipeline Initial processing steps include matching the raw reads to experimental samples, trimming off the tag and primer portions, 
and removing sequences of low quality. If the gene chosen is 16S RNA, the orientation of sequences will be checked and reverse complemented if needed.
Initial processing requires a sequence file and at least one forward primer. The required sequence file, 
obtained from the pyrosequencing center, and can be in FASTA format or SFF Format (which contains both the sequence and quality information).

<a name="Tutorial"></a>
The sample input and output files can be downloaded from RDP tutorial http://rdp.cme.msu.edu/tutorials/init_process/RDPtutorial_INITIAL-PROCESS.html.

### Setup
This project depends on ReadSeq, AlignmentTools, and ProbeMatch. See RDPTools (https://github.com/rdpstaff/RDPTools) to install.

### Usage

* Run Initial Processing

		java -jar /path/SeqFilters.jar
		USAGE: InitialProcessorMain <options>
 		-f,--forward-primers <arg>   Comma seperated list of forward primers
 		-F,--max-forward <arg>       Maximum forward edit distance (default=2)
 		-g,--gene-name <arg>         Gene name, possible values are RRNA16S and OTHER (default=RRNA16S)
 		-m,--min-length <arg>        Minimum sequence length after primer triming (default=0)
 		-n,--max-ns <arg>            Maxmimum number of Ns allowed in a sequence default=0)
 		-o,--outdir <arg>            Output directory (default=cwd)
 		-O,--result-dir-name <arg>   Result dir name (default=result_dir)
 		-p,--keep-primer <arg>       Don't trim primers (default=false)
		 -Q,--min-qual <arg>          Minimum sequence length after primer triming (default=20)
 		-q,--qual-file <arg>         Quality input file (default=null)
 		-r,--reverse-primers <arg>   Comma seperated list of reverse primers (default=null)
 		-R,--max-reverse <arg>       Maximum reverse edit distance (default=0)
 		-s,--seq-file <arg>          Sequence file to process
 		-S,--skip-notag              Don't process no tag sequences
 		-t,--tag-file <arg>          Tag file (default=null)

 		
 	An example command using the data from RDP tutorial, [see link above](#Tutorial). 
 	 		
 		java -jar /path/to/SeqFilters.jar --forward-primers AYTGGGYDTAAAGNG --max-forward 2 --reverse-primers CCGTCAATTCMTTTRAGT --max-reverse 1 --seq-file 1.TCA.454Reads.fna --qual-file 1.TCA.454Reads.qual --min-length 300 --tag-file region1_tag.txt --outdir initial_process


