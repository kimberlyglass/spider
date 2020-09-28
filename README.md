# SPIDER

Implementation of the SPIDER (Seeding PANDA Interactions to Derive Epigenetic Regulation) network reconstruction approach.
Authors: Abhijeet R. Sonawane; Kimberly Glass


SPIDER takes three primary inputs, including the location of:
(1) potential transcription factor binding sites, which can be defined by position weight matrices mapped onto the DNA
(2) open chromatin regions, based on epigenetic data
(3) gene regulatory regions, which can be defined based on proximity to transcriptional start sites. 


In the code, these are designated by the following variables:
(1) motifdir: A directory containing a series of BED files, one for each TF, listing the location where its corresponding sequence motif can be found in the genome. The names of these files should be of the format <TFname>.bed, the columns in these files are of the form: <chr>\t<start>\t<end>
(2) epifile: A BED file containing regions of open chromatin. File format is: <chr>\t<start>\t<end> (all other columns ignored)
(3) regfile: A BED file containing the regulatory regions of genes. File format is: <chr>\t<start>\t<end>\t<gene-name>

In addition to these input files, there is one tuning parameter associated with the PANDA message-passing framework:
alpha: by default is set equal to 0.1. Must be between 0 and 1, values between 0.05 and 0.25 are recommended.


The SPIDER code outputs two files containing the (1) epigenetically-informed motif information, and (2) SPIDER-predicted regulatory network. These are encoded by the variables:
(1) motifhitfile: A single BED file created from intersecting the locations of TF motifs in the genome with open chromatin data.
(2) outtag: A file containing the predicted SPIDER network.

IMPORTANT: bedtools must be installed in order to run SPIDER. Please see https://bedtools.readthedocs.io/ for more information.
If bedtools is not already on the system path, please update the following variable to indicate where it is installed: bedtoolspath



The SPIDER algorithm can be summarized by four main steps:

(1) STEP 1: Intersect open chromatin regions with motif locations:
Function: CreateEpigeneticMotif()
Input-parameter(s): motifdir, epifile
Summary: In this step, SPIDER uses bedtools to intersect a BED file containing regions of open chromatin with a series of BED files containing the locations of TF motifs (one BED file per TF).
Output: The output of this step is a single BED file that contains the locations of TF motifs that are in open chromatin regions. By default, each of these locations is given a default score of one. Note that, in practice, the file produced by this step could be produced in another manner and still used by SPIDER.

(2) STEP 2: Intersect motifs in open chromatin (from Step 1) with gene regulatory regions and create a seed regulatory network:
Function: BuildSPIDERprior()
Input-parameter(s): regfile, motifhitfile (created in STEP 1)
Summary: In this step, SPIDER uses bedtools to overlap a BED file containing the locations of TFs that are in open chromatin (created in Step 1) with a BED file containing the regulatory regions of genes. Note that a gene can have multiple associated regulatory regions in this second file. If a TFâ€™s motif falls within the regulatory region(s) of a gene, then an edge is created between that TF and gene. The maximum score across all TF motif instances associated with a gene is used to weight the edge; by default, this value is one.
Output: The output of this step is an epigenetically-informed seed regulatory network between all transcription factors and genes. 

(3) STEP 3: Degree normalize seed network
Function: DegreeAdjust()
Input-parameter(s): seed network created in STEP 2
Summary: The weights of edges in the seed network are degree-normalized to emphasize connections to high degree TFs and genes, based on the transformation A'=A*sqrt(k^2+l^2), where A is the seed network adjacency matrix, k is the degree of the TF in the seed network, and l is the degree of the gene in the seed network.
Output: The output of the step is an adjacency matrix where edges between high degree TFs and genes have been upweighted.

(4) STEP 4: Apply message-passing
Function: PANDA()
Input-parameter(s): alpha, degree-normalized network created in STEP 3
Summary: SPIDER applies the PANDA message-passing algorithm to the degree-normalized seed network A* calculated in Step 3. PANDAâ€™s message-passing framework integrates information from three networks, representing TF protein-protein interaction (P), TF/gene regulation (W) and gene co-expression (C). In SPIDER, P and C are set equal to the identity matrix; W is set equal to A*.
Output: This step returns a complete, bi-partite network with edge weights representing the likelihood that a TF regulates a gene; the distribution of these edge weights is similar to z-scores. 
