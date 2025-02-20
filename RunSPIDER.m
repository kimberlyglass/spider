%%%% PARAMETER REGION %%%%

% Input Data
motifdir='ToyData/MotifBedFiles_subset/'; % where the original motif scan files are stored (one bed file per motif)
epifile='ToyData/K562_DnasePeaks.bed'; % bed file with open chromatin regions
regfile='ToyData/RegulatoryRegions_0-1kb_subset.bed'; % bed file containing regulatory regions for genes

% message-passing parameter
alpha=0.1; % level of message-passing

% output files
motifhitfile='ToyData/dnasefiltered_motiflocations.bed'; % intermediate file storing epigenetically informed motif information; created with CreateEpigeneticMotif.m
outtag='ToyData_FinalSPIDERNetwork.pairs'; % name of file to print the final network info

% location of codes needed to run SPIDER
bedtoolspath=''; % where bedtools is installed, set equal to '' if bedtools is already on the system path
funcpath='./SPIDER/'; % path to where all the SPIDER functions are stored
addpath(funcpath);


%%%% Run SPIDER %%%%

% STEP 1: Create a file containing epigenetically-informated motif information
CreateEpigeneticMotif(epifile, motifdir, motifhitfile, bedtoolspath);

% STEP 2: Build SPIDER seed network
[PriorNet, TFNames, GeneNames]=BuildSPIDERprior(motifhitfile, regfile, bedtoolspath);

% STEP 3: Degree-normalize seed network
DANet=DegreeAdjust(PriorNet);

% STEP 4: Run message-passing
SpiderNet=PANDA(DANet, eye(length(GeneNames)), eye(length(TFNames)), alpha);


%%%% Print SPIDER-predcted Network to file %%%%

% reshape network information into vectors 
TF=repmat(TFNames, 1, length(GeneNames)); TF=TF(:);
gene=repmat(GeneNames', length(TFNames), 1); gene=gene(:);
PriorNet=PriorNet(:);
SpiderNet=SpiderNet(:);

% print to file
fid=fopen(outtag, 'wt');
fprintf(fid, 'TF\tgene\tSeed-network\tSPIDER-prediction\n');
for(cnt=1:length(TF))
	fprintf(fid, '%s\t%s\t%f\t%f\n', TF{cnt}, gene{cnt}, PriorNet(cnt), SpiderNet(cnt));
end
fclose(fid);
