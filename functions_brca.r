# protein aggregate ---------
protein_aggregate=function(PSM,ctrl){
  filter(PSM,
         `Confidence`=='High', 
         `Number of Proteins`==1,        # remove PSMs assigned to multiple proteins
         `Quan Info`!='NoQuanLabels',    # remove non unique quantification
         `Abundance 131C`>10,            # PIS filter
         Contaminant=='FALSE'            # remove contaminants
  ) -> PSM 
  colnames(PSM)=gsub('Abundance ','',colnames(PSM))
  PSM$Sequence=str_match(toupper(PSM$`Annotated Sequence`),'\\.([A-Z]+)')[,2]
  PSM[,grep('^1',colnames(PSM))][is.na(PSM[,grep('^1',colnames(PSM))])]=0
  PSM[,grep('^1',colnames(PSM))][PSM[,grep('^1',colnames(PSM))]==0]=0.01
  # calculate the isodoping normalization factor by using both isodoped and naive peptides from isodoped proteins  
  isodoping_factor= median(PSM[PSM$Sequence %in% isodoping$PeptideSequence,'131C'],na.rm=T)/
    median(PSM[PSM$`Master Protein Accessions` %in% isodoping$uniprot &
                 !(PSM$Sequence %in% isodoping$PeptideSequence),'131C'],na.rm=T)
  # correct isodoping peptides abundances  
  PSM[PSM$Sequence %in% isodoping$PeptideSequence,'131C']=
    PSM[PSM$Sequence %in% isodoping$PeptideSequence,'131C']/isodoping_factor
  # channel total intensity
  sum_int=apply(PSM[,grep('^1',colnames(PSM))],2,function(X) sum(X,na.rm=T))
  PSM[,grep('^1',colnames(PSM))]= t(t(PSM[,grep('^1',colnames(PSM))]) * (100000000/sum_int))
  # normalize on PIS 
  PSM[,grep('^1',colnames(PSM))]=t(apply(PSM[,grep('^1',colnames(PSM))],1,function(X) X/X[11])) 
  PSM=ddply(PSM,.(`Master Protein Accessions`),
            summarize,
            `126`=median(  `126`[order(`Average Reporter SN`,decreasing = T)[1:5]],na.rm=T),
            `127N`=median(`127N`[order(`Average Reporter SN`,decreasing = T)[1:5]],na.rm=T),
            `127C`=median(`127C`[order(`Average Reporter SN`,decreasing = T)[1:5]],na.rm=T),
            `128N`=median(`128N`[order(`Average Reporter SN`,decreasing = T)[1:5]],na.rm=T),
            `128C`=median(`128C`[order(`Average Reporter SN`,decreasing = T)[1:5]],na.rm=T),
            `129N`=median(`129N`[order(`Average Reporter SN`,decreasing = T)[1:5]],na.rm=T),
            `129C`=median(`129C`[order(`Average Reporter SN`,decreasing = T)[1:5]],na.rm=T),
            `130N`=median(`130N`[order(`Average Reporter SN`,decreasing = T)[1:5]],na.rm=T),
            `130C`=median(`130C`[order(`Average Reporter SN`,decreasing = T)[1:5]],na.rm=T),
            `131N`=median(`131N`[order(`Average Reporter SN`,decreasing = T)[1:5]],na.rm=T),
            `131C`=median(`131C`[order(`Average Reporter SN`,decreasing = T)[1:5]],na.rm=T))
  colnames(PSM)[1]='Accession'
  filter(PSM,Accession!='')-> PSM
  return(PSM)
}


# protein merger --------
merger=function(mat){
  master_prot=mat
  for( idx in 1:length(mat)){
    colnames(mat[[idx]])[2:12]=brca_id[brca_id$set==idx,'id']
    mat[[idx]]= mat[[idx]][,grep('set',colnames( mat[[idx]]))]
  }
  # get protein ids 
  prot.id=lapply(master_prot, function(X) X['Accession'])
  # take all IDS
  lab.id=as.character(unique(unlist(prot.id)))
  lab.id=lab.id[lab.id!='']
  # add them to normalized frames
  mat=Map(cbind, mat, Accession = prot.id)
  # create new frame with just Accession
  mat[[length(master_prot)+1]]=data.frame(Accession=lab.id,stringsAsFactors = F)
  # reorder
  mat=mat[c(length(master_prot)+1,1:length(master_prot))]
  # merge all batches with left_join
  mat= mat %>% purrr::reduce(left_join, by = "Accession")
  # assign uniprot to rownames and remove the column
  rownames(mat)=mat$Accession
  mat=mat[,- 1]
  return(mat)
}


# reads in PSMs files and apply PSMaggreate 
read_process= function(filename,ctrl)  {
  PSM=fread(filename,data.table = F)
  PSM=protein_aggregate(PSM,ctrl)
  return(PSM) 
}

# function to aggregate PSMs 
PSMaggregate=function(PSM,ctrl){
  # remove non unique PSMS
  filter(PSM,
         `Confidence`=='High', 
         `Number of Proteins`==1,        # remove PSMs assigned to multiple proteins
         `Quan Info`!='NoQuanLabels',    # remove non unique quantification
         `Abundance 131C`>10,            # PIS filter
         Contaminant=='FALSE'            # remove contaminants
  ) -> PSM 
  colnames(PSM)=gsub('Abundance ','',colnames(PSM))
  PSM$Sequence=str_match(toupper(PSM$`Annotated Sequence`),'\\.([A-Z]+)')[,2]
  PSM[,grep('^1',colnames(PSM))][is.na(PSM[,grep('^1',colnames(PSM))])]=0
  PSM[,grep('^1',colnames(PSM))][PSM[,grep('^1',colnames(PSM))]==0]=0.01
  # calculate the isodoping normalization factor by using both isodoped and naive peptides from isodoped proteins  
  isodoping_factor= median(PSM[PSM$Sequence %in% isodoping$PeptideSequence,'131C'],na.rm=T)/
    median(PSM[PSM$`Master Protein Accessions` %in% isodoping$uniprot &
                 !(PSM$Sequence %in% isodoping$PeptideSequence),'131C'],na.rm=T)
  # correct isodoping peptides abundances  
  PSM[PSM$Sequence %in% isodoping$PeptideSequence,'131C']=
    PSM[PSM$Sequence %in% isodoping$PeptideSequence,'131C']/isodoping_factor
  # channel total intensity
  sum_int=apply(PSM[,grep('^1',colnames(PSM))],2,function(X) sum(X,na.rm=T))
  PSM[,grep('^1',colnames(PSM))]= t(t(PSM[,grep('^1',colnames(PSM))]) * (100000000/sum_int))
  # normalize on PIS 
  PSM[,grep('^1',colnames(PSM))]=t(apply(PSM[,grep('^1',colnames(PSM))],1,function(X) X/X[11])) 
  PSM=ddply(PSM,.(`Master Protein Accessions`,`Sequence`), 
            summarize,
            `126`=median(  `126`[order(`Average Reporter SN`,decreasing = T)[1:5]],na.rm=T),
            `127N`=median(`127N`[order(`Average Reporter SN`,decreasing = T)[1:5]],na.rm=T),
            `127C`=median(`127C`[order(`Average Reporter SN`,decreasing = T)[1:5]],na.rm=T),
            `128N`=median(`128N`[order(`Average Reporter SN`,decreasing = T)[1:5]],na.rm=T),
            `128C`=median(`128C`[order(`Average Reporter SN`,decreasing = T)[1:5]],na.rm=T),
            `129N`=median(`129N`[order(`Average Reporter SN`,decreasing = T)[1:5]],na.rm=T),
            `129C`=median(`129C`[order(`Average Reporter SN`,decreasing = T)[1:5]],na.rm=T),
            `130N`=median(`130N`[order(`Average Reporter SN`,decreasing = T)[1:5]],na.rm=T),
            `130C`=median(`130C`[order(`Average Reporter SN`,decreasing = T)[1:5]],na.rm=T),
            `131N`=median(`131N`[order(`Average Reporter SN`,decreasing = T)[1:5]],na.rm=T),
            `131C`=median(`131C`[order(`Average Reporter SN`,decreasing = T)[1:5]],na.rm=T))
  colnames(PSM)[1]='Accession'
  filter(PSM,Accession!='')-> PSM
  return(PSM)
}


# reads in PSMs files and apply PSMaggreate 
read_process_psm= function(filename,ctrl)  {
  PSM=fread(filename,data.table = F)
  PSM=PSMaggregate(PSM,ctrl)
  return(PSM) 
}

# psm merger
merger_psm=function(mat){
  # get protein ids 
  pep.sequence=lapply(mat, function(X) X['Sequence'])
  accession=lapply(mat, function(X) X['Accession'])
  id=data.frame(Sequence=unlist(pep.sequence),Accession=unlist(accession))
  id=unique(id)
  for( idx in 1:length(mat)){
    colnames(mat[[idx]])[3:13]=brca_id[brca_id$set==idx,'id']
    mat[[idx]]= mat[[idx]][,grep('set',colnames( mat[[idx]]))]
  }
  # take all IDS
  lab.id=as.character(unique(unlist(pep.sequence)))
  lab.id=lab.id[lab.id!='']
  # # add them to normalized frames
  mat=Map(cbind, mat, Sequence = pep.sequence)
  # # create new frame with just Accession
  mat[[length(pep.sequence)+1]]=data.frame(Sequence=lab.id,stringsAsFactors = F)
  # # reorder so it's the first df in the list of df
  mat=mat[c(length(pep.sequence)+1,1:length(pep.sequence))]
  # # merge all batches with left_join
  mat= mat %>% purrr::reduce(left_join, by = "Sequence")
  mat$Accession=id$Accession[match(mat$Sequence,id$Sequence)]  
  # # # assign uniprot to rownames and remove the column
  rownames(mat)=mat$Sequence
  mat=mat[,- 1]
  return(mat)
}


# Aggregate PSMs without normalization 
PSMaggregate_nonorm=function(PSM,ctrl){
  # remove non unique PSMS
  filter(PSM,
         `Quan Info`!='NoQuanLabels', # remove non unique quantification
         Contaminant=='FALSE'  # remove contaminants
  ) -> PSM
  colnames(PSM)=gsub('Abundance ','',colnames(PSM))
  # aggregate psms by sum
  PSM=aggregate(cbind(`126`,`127N`,`127C`,`128N`,`128C`,`129N`,`129C`,`130N`,`130C`,`131N`,`131C`)~
                  `Annotated Sequence`,
                data = PSM,
                FUN= sum,
                na.action=na.pass,
                na.rm=TRUE)
  colnames(PSM)[1]='Sequence'
  return(PSM)
}




## merger functions ------
merger_psm_nonorm=function(mat){
  # get protein ids 
  pep.sequence=lapply(mat, function(X) X['Sequence'])
  id=data.frame(Sequence=unlist(pep.sequence))
  id=unique(id)
  for( idx in 1:length(mat)){
    colnames(mat[[idx]])[2:12]=brca_id[brca_id$set==idx,'id']
    mat[[idx]]= mat[[idx]][,grep('set',colnames( mat[[idx]]))]
  }
  # take all IDS
  lab.id=as.character(unique(unlist(pep.sequence)))
  lab.id=lab.id[lab.id!='']
  # add them to normalized frames
  mat=Map(cbind, mat, Sequence = pep.sequence)
  # create new frame with just Accession
  mat[[length(pep.sequence)+1]]=data.frame(Sequence=lab.id,stringsAsFactors = F)
  # reorder so it's the first df in the list of df
  mat=mat[c(length(pep.sequence)+1,1:length(pep.sequence))]
  # merge all batches with left_join
  mat= mat %>% purrr::reduce(left_join, by = "Sequence")
  # assign uniprot to rownames and remove the column
  rownames(mat)=mat$Sequence
  mat=mat[,- 1]
  return(mat)
}

## read functions -----
read_process_nonorm=function(filename,ctrl)  {
  PSM=fread(filename,data.table = F)
  PSM=PSMaggregate_nonorm(PSM,ctrl)
  return(PSM) 
}



### assign uniprot identifier to peptide sequences 
num_pep_x_protein=function(filenames){
  # digest
  proteome=readAAStringSet(file='/input/uniprot_homo_20180803_brca.fasta')
  ms1=data.frame(cleave(proteome, custom=c("K|R", "K(?=P)"), missedCleavages = 0, unique = TRUE))
  ms2=data.frame(cleave(proteome, custom=c("K|R", "K(?=P)"), missedCleavages = 1, unique = TRUE))
  ms3=data.frame(cleave(proteome, custom=c("K|R", "K(?=P)"), missedCleavages = 2, unique = TRUE))
  ms4=data.frame(cleave(proteome, custom=c("K|R", "K(?=P)"), missedCleavages = 3, unique = TRUE))
  df=rbind(ms1,ms2,ms3,ms4);dim(df)
  rm(ms1);rm(ms2);rm(ms3);rm(ms4);gc()
  df$size=unlist(lapply(df$value,function(X) nchar(X)))
  # filter size
  min_aa=5
  max_aa=60
  filter(df,size >= min_aa,size <=max_aa) -> df
  df$symbol=str_match(df$group_name,pattern = 'GN=(.+) PE')[,2]
  df$uniprot=str_match(df$group_name,pattern = 'sp\\|(.+)\\|')[,2]
  df=df[,c('value','uniprot')]
  # merge shared peptides by sequence 
  peptide_uniprot_map=data.frame(setDT(df)[,
                                           lapply(.SD, function(x) paste(x,collapse=';') ),
                                           by=list(`value`),
                                           .SDcols = c('uniprot')])
  # load psm files
  psm.merge = lapply(filenames, function(x) fread(x, data.table = F)) 
  psm.merge = do.call("rbind", psm.merge) 
  filter(psm.merge,  !is.na(`Abundance 131C`), `Quan Info`!='NoQuanLabels',Contaminant=='FALSE' ) -> psm.merge
  psm.merge$Sequence=str_match(toupper(psm.merge$`Annotated Sequence`),'\\.([A-Z]+)')[,2]
  peptides.identified=data.frame(Sequence=unique(psm.merge$Sequence));dim(peptides.identified)
  # match peptide sequence to uniprot ID
  peptides.identified$Accession=peptide_uniprot_map$uniprot[match(peptides.identified$Sequence,peptide_uniprot_map$value)]
  peptidesXprotein=lapply(peptides.identified$Accession,function(X)strsplit(X,';'))
  peptidesXprotein=table(unlist(peptidesXprotein))
  return(peptidesXprotein)
}



## calculate proteins stats and number of psms x plex ----
protein_stats=function(filenames,filenames_prot){
# load psm files
psm.merge = lapply(filenames, function(x) fread(x, data.table = F)) 
psm.merge = do.call("rbind", psm.merge) 
filter(psm.merge,
       `Number of Proteins`==1,        # remove PSMs assigned to multiple proteins
       `Quan Info`!='NoQuanLabels',    # remove non unique quantification
       `Abundance 131C`>10,            # PIS filter
       Contaminant=='FALSE'            # remove contaminants
) -> psm.merge.filtered
# add set name
psm.merge.filtered$set=tolower(str_match(psm.merge.filtered$`Spectrum File`, "(Set_[0-9]+_)")[,2])
# count psm x protein x plex
group_by(psm.merge.filtered,`Master Protein Accessions`,set) %>% summarize(n_psm=n()) -> psm.plex
pivot_wider(psm.plex,names_from = set, values_from = n_psm) %>% data.frame -> psm.plex
rownames(psm.plex)=psm.plex$Master.Protein.Accessions
psm.plex=psm.plex[,-1]
psm.plex=psm.plex[,order(as.numeric(str_match(colnames(psm.plex), "set_([0-9]+)_")[,2]))]
# remove psms counts for protein excluded by FDR filtering
for( filename in filenames_prot){
  Proteins=fread(filename,data.table=F)
  psm.plex[!rownames(psm.plex) %in%
                       Proteins$Accession[grepl('High|Medium',Proteins$`Protein FDR Confidence Combined`)],
                     grepl(tolower(str_match(filename,'BrCa_(Set_[0-9]+_)')[,2]),colnames(psm.plex))]=NA
}
colnames(psm.plex)=gsub('_$','_number_PSMs',colnames(psm.plex))
# merge protein df and psm count df
index=match(rownames(protein),rownames(psm.plex))
protein.final=cbind(protein,psm.plex[index,])
# add total number of peptides, unique peptides, and PSMS
protein.stat=data.table("Master Protein Accessions"=unique(psm.merge.filtered$`Master Protein Accessions`))
# assign  gene symbol
protein.merge = lapply(filenames_prot, function(x) fread(x, data.table = F,select =  c('Accession','Description'))) 
protein.merge = unique(do.call("rbind", protein.merge))
protein.merge$Name=str_match(protein.merge$Description, "GN=(.+) PE")[,2]
protein.stat$Symbol=protein.merge$Name[match(gsub('-[0-9]','',protein.stat$`Master Protein Accessions`),protein.merge$Accession)]

psm.merge.filtered$Sequence=str_match(toupper(psm.merge.filtered$`Annotated Sequence`),'\\.([A-Z]+)')[,2]
#number of peptides  
protein.stat$Number_peptides=as.numeric(peptidesXprotein[match(protein.stat$`Master Protein Accessions`,names(peptidesXprotein))])
#number of unique peptides  
Number_unique_peptides=setDT(filter(psm.merge.filtered,`Number of Proteins`==1))[, .(Number_unique_peptides = length(unique(`Sequence`))), by = "Master Protein Accessions"]
protein.stat=full_join(protein.stat,Number_unique_peptides,by="Master Protein Accessions")
#number of total PSMs  
Number_PSMs=setDT(psm.merge.filtered)[, .(Number_PSMs = .N) , by = `Master Protein Accessions`]
protein.stat=full_join(protein.stat,Number_PSMs,by="Master Protein Accessions")
colnames(protein.stat)[1]='Accession'
#protein.summary=protein_stats(psm.merge)
index=match(rownames(protein.final),protein.stat$Accession)
protein.final=cbind(protein.stat[index,],protein.final)
return(protein.final)
}



