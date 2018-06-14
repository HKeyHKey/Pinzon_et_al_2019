1. The list of metazoan "eukaryotic RdRP" proteins identified in Figure 1 was supplemented with sequences directly retrieved from the PFAM database for fungi, plants and other eukaryotic lineages. Sequence names are arbitrary but the capital letters correspond to species identifiers according to the Uniprot nomenclature (http://www.uniprot.org/help/taxonomy#organism-denomination). The complete sequence set is given in file 'RdRPMARSV1.fa'.

2. Sequences in 'RdRPMARSV1.fa' were aligned using HMMalign with the eukaryotic RdRP HMM profile obtained from PFAM (PFAM#05183). Sequences that were missing more than 40 amino acid at the N or C terminus of the domain (as defined by PFAM), or sequences that showed more than 20 adjacent amino acids missing inside the domain, were deleted from the alignment. Then, sites used for further phylogenetic reconstruction were selected using trimAl with the default parameters. The final alignment is given in file 'RDRPMARSV1alignHMMALIGNShortV3TRIM.fasta'.

3. The Bayesian inference tree was constructed using MrBayes 3.2.6 on that alignment file. MrBayes output files are copied here for convenience:
'done.txt', 'infile.nex', 'infile.nex.ckp', 'infile.nex.ckp%7E', 'infile.nex.con.tre', 'infile.nex.mcmc', 'infile.nex.parts', 'infile.nex.run1.p', 'infile.nex.run2.p', 'infile.nex.trprobs', 'infile.nex.tstat', 'infile.nex.vstat', '_JOBINFO.TXT', 'paramfile.txt', 'scheduler.conf', '_scheduler_stderr.txt', 'start.txt', 'STDERR', 'stderr.txt', 'STDOUT', 'stdout.txt' and 'term.txt'.
Note that MrBayes also generated two other files ('infile.nex.run1.t' and 'infile.nex.run2.t.'), which could not be copied on GitHub because of size limitations. These two files are available (in bz2-compressed form) at:

https://www.igh.cnrs.fr/images/microsite/herve-seitz/files/pinzon-et-al-2018/02_RdRP_phylogeny/infile.nex.run1.t.bz2

https://www.igh.cnrs.fr/images/microsite/herve-seitz/files/pinzon-et-al-2018/02_RdRP_phylogeny/infile.nex.run2.t.bz2
