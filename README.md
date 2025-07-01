# phytoSentinel

Principle - translate complete genome and get the position of conserved domains associated w/ NLRs and PRRs.

## To do list

### 1. Select domains 

Select hmm models with NLR/PRR domains (check NLR Express)

### 2. Training data 

1) get plant genomes and their annotation;
2) translate all material and characterize all domains via hmm;
3) get all NLRs from annotation and their position;
4) create the following datasets: 
- a set of NLR/PRR domains and their position along the plant genome (chr/domainID/start/stop)  - filter by max bitscore;
- the position of each NLR/PRR gene in genome and their domain architecture (e.g. CC-NBD-ARC1-ARC2-LRR);

### 3. How it will work (?)

I am thinking in the following steps:
1. Search for domains in genome;
2. Get the position of all domains and some flanking region;
3. Use ML to find the start/end postitions of the region with domains;
4. Execute BLASTn w/ custom database of NLR/PRRs genes;
5. Predict exonic and intronic regions (?)
