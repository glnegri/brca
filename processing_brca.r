rm(list=ls())
# load packages
library(data.table)
library(cleaver)
library(tidyverse)
library(stringr)
library(parallel)

setwd('workspace')

# load functions 
source("functions_brca.R")


# create brca_id
brca_id=data.frame(
  label=rep(c("126",  "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C"),38),
  set=gl(n = 38,k = 11),
  id=rep(1:11,38))
brca_id$id=unlist(apply(brca_id,1,function(X) gsub(' ','',paste(c('set_',X[2],'_',X[3]),collapse = '')) ))


# load search results
filenames=list.files(path = 'input/',include.dirs = T,recursive = T,full.names = T)
filenames_prot=filenames[grep('_Proteins.txt$',filenames)]
filenames_prot=filenames_prot[order(as.numeric(str_match(filenames_prot, "Set_([0-9]+)")[,2]))]
filenames=filenames[grep('PSM',filenames)]
filenames=filenames[order(as.numeric(str_match(filenames, "Set_([0-9]+)")[,2]))]


# isodoping library
isodoping=fread('input/isodoping_library_brca.txt',data.table = F)


## merge and normalize -----
Sys.time()
cl <- makeCluster(detectCores())
clusterEvalQ(cl, {library(data.table); library(ggplot2);library(tidyverse);library(plyr)})
#export variables in all nodes
clusterExport(cl,c("read_process","PSMaggregate","protein_aggregate",'read_process_psm',
                   "read_process_nonorm","PSMaggregate_nonorm",
                   "isodoping"))

master <- parLapply(cl, filenames, function(filename) read_process(filename,ctrl = 'pis'))
protein=merger(master)
master_psm <- parLapply(cl, filenames, function(filename) read_process_psm(filename,ctrl = 'pis'))
peptide=merger_psm(master_psm)
master_noN <- parLapply(cl, filenames, function(filename) read_process_nonorm(filename,ctrl = 'pis'))
psmNN=merger_psm_nonorm(master_noN)
stopCluster(cl)


# total S/N abundance before normalization
total_abundance=apply(psmNN[,1:(38*11)],2,function(X) sum(X,na.rm=T))
# samples with low abundance 
qc.exclude= names(total_abundance[total_abundance<2e06]);qc.exclude # exclude samples with total S/N < 2e06
 

# assign uniprot identifier to peptide sequences and count peptides X protein
peptidesXprotein=num_pep_x_protein(filenames = filenames )

 
# filter proteins
filenames_prot
protein_fdr=protein
for( filename in filenames_prot){
  Proteins=fread(filename,data.table=F)
  protein_fdr[!rownames(protein_fdr) %in%
                Proteins$Accession[grepl('High|Medium',Proteins$`Protein FDR Confidence Combined`)],
              grepl(tolower(str_match(filename,'BrCa_(Set_[0-9]+_)')[,2]),colnames(protein_fdr))]=NA
}
log.index=(unlist(apply(protein_fdr,1,function(X) length(X[!is.na(X)])>0 )));table(log.index)
protein_fdr=protein_fdr[log.index,]
dim(protein_fdr);dim(na.omit(protein_fdr))
protein=protein_fdr;dim(protein)

# keep only proteins with at least 2 peptides, remove PIS channels, log2 transform
protein = log2(protein[rownames(protein) %in% names(peptidesXprotein[peptidesXprotein > 1]), !grepl('_11$', colnames(protein))]);dim(protein)

peptide.total=peptide
#keep only peptides from filtered proteins, remove PIS channels, log2 transform
peptide=mutate_if(peptide[peptide$Accession %in% rownames(protein),!grepl('_11$', colnames(peptide))], is.numeric, log2);dim(peptide)

protein.final=protein_stats(filenames ,filenames_prot)


# write protein table
writexl::write_xlsx(protein.final,'protein_final.xlsx')
peptide$Sequence=rownames(peptide)
writexl::write_xlsx(peptide[,c(381,382,1:380)],'peptide_final.xlsx')






