# SPIDER

Implementation of the SPIDER (Seeding PANDA Interactions to Derive Epigenetic Regulation) network reconstruction approach.\
Authors: Abhijeet R. Sonawane; Kimberly Glass

<space>\
<space>

SPIDER takes three primary inputs, including the location of:\
(1) potential transcription factor binding sites, which can be defined by position weight matrices mapped onto the DNA\
(2) open chromatin regions, based on epigenetic data\
(3) gene regulatory regions, which can be defined based on proximity to transcriptional start sites.

The output of SPIDER is a file containing an epigenetically informed gene regulatory network.

IMPORTANT: bedtools must be installed in order to run SPIDER. Please see https://bedtools.readthedocs.io/ for more information. 

For additional information on SPIDER, please visit: https://sites.google.com/a/channing.harvard.edu/kimberlyglass/tools/spider

<space>\
<space>
  
<space>\
<space>

The SPIDER algorithm can be summarized by four main steps:

(1) STEP 1: Intersect open chromatin regions with motif locations

In this step, SPIDER uses bedtools to intersect a BED file containing regions of open chromatin with a series of BED files containing the locations of TF motifs (one BED file per TF). The output of this step is a single BED file that contains the locations of TF motifs that are in open chromatin regions. By default, each of these locations is given a default score of one. Note that, in practice, the file produced by this step could be produced in another manner and still used by SPIDER.  

<space>\
<space>
  
(2) STEP 2: Intersect motifs in open chromatin (from Step 1) with gene regulatory regions and create a seed regulatory network

In this step, SPIDER uses bedtools to overlap a BED file containing the locations of TFs that are in open chromatin (created in Step 1) with a BED file containing the regulatory regions of genes. Note that a gene can have multiple associated regulatory regions in this second file. If a TF's motif falls within the regulatory region(s) of a gene, then an edge is created between that TF and gene. The maximum score across all TF motif instances associated with a gene is used to weight the edge; by default, this value is one. The output of this step is an epigenetically-informed seed regulatory network between all transcription factors and genes.  

<space>\
<space>


(3) STEP 3: Degree normalize seed network

In this step, the weights of edges in the seed network are degree-normalized to emphasize connections to high degree TFs and genes, based on the transformation A'=A\*sqrt(k^2+l^2), where A is the seed network adjacency matrix, k is the degree of the TF in the seed network, and l is the degree of the gene in the seed network. The output of the step is an adjacency matrix where edges between high degree TFs and genes have been upweighted.  

<space>\
<space>

(4) STEP 4: Apply message-passing

In this step, SPIDER applies the PANDA message-passing algorithm to the degree-normalized seed network A* calculated in Step 3. PANDA's message-passing framework integrates information from three networks, representing TF protein-protein interaction (P), TF/gene regulation (W) and gene co-expression (C). In SPIDER, P and C are set equal to the identity matrix; W is set equal to A*. This step returns a complete, bi-partite network with edge weights representing the likelihood that a TF regulates a gene; the distribution of these edge weights is similar to z-scores. 
