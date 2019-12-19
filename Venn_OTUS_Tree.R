#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

setwd(paste("./",args[1],sep = ''))

#library("futile.logger", lib.loc="/ngs/mparras/lib_R")
#library("VennDiagram", lib.loc="/ngs/mparras/lib_R") #The script needed is VennDiagram
if (!require("VennDiagram")) install.packages("VennDiagram", repos='http://cran.us.r-project.org')

Cores_Tree<-read.table('./Tree/results.txt',sep="\t",header=T,na.strings=c("","NA"),fill=TRUE) #First, we get the information from Tree results
Reads_cores_Tree<-data.frame(matrix(nrow = dim(Cores_Tree)[1],ncol = 2))

files_read<-list.files(path="./Tree/",full.names = TRUE, recursive = TRUE)

for (i in 1:length(files_read))
{
  
  if(grepl("otus.txt.modified.txt",files_read[i]) & file.info(files_read[i])$size > 0) #The files with this termination have the info of the reads from each OTU
  {
    Tree_reads_core<-read.table(files_read[i],sep="\t")
  }
}


for(i in 1:(dim(Cores_Tree)[1]-1)) #With this loop we get all the sequence names from each core node that was anotated into the results.txt file
{
  Cores_list=strsplit(x = as.character(Cores_Tree[i,9]),split = ";")
  new_info_core<-NULL
  
  for(o in 1:length(Cores_list[[1]]))
  {
    new_info_core<-c(new_info_core,as.character(Tree_reads_core[Tree_reads_core$V1==Cores_list[[1]][o],2]))  
  }

  new_info_core<-paste(new_info_core,collapse = "")
  
  Reads_cores_Tree[i,1]<-as.character(Cores_Tree[i,1])
  Reads_cores_Tree[i,2]<-new_info_core
  
}

delete_row=length(Reads_cores_Tree[,1]) #We erase the last row as it is empty
Reads_cores_Tree=Reads_cores_Tree[-delete_row,]

OTU_cores<-read.table('./OTUs/Distr.txt',sep="\t") #The second thing wi will do is obtain the core information from denovo OTUs picking
Cores_list_OTUs<-colnames(OTU_cores)
OTUs_reads_core<-data.frame(matrix(nrow = length(Cores_list_OTUs),ncol = 2))

Cores_list_OTUs=strsplit(x = as.character(Cores_list_OTUs),split = "\\.")          

files_read<-list.files(path="./OTUs/",full.names = TRUE, recursive = TRUE)


for (o in 1:length(Cores_list_OTUs))
{
  for (i in 1:length(files_read))
  {
    if(grepl("otus.txt.modified.txt",files_read[i]) & file.info(files_read[i])$size > 0) #As we did with Tree, we need to find all the files with this termination into all the OTUs levels folder evaluated, that have the info off the reads contained into each OTU 
    {
      name2=strsplit(files_read[i],'\\/')[[1]][4]
      name2=strsplit(name2,'\\.')[[1]][2]
      
      if(name2==Cores_list_OTUs[[o]][3]) #Here we check that the file read match with any level that have core OTUs
      {
        OTUs_reads_cores<-read.table(files_read[i],sep="\t")
        new_info_core<-NULL
        new_info_core<-c(new_info_core,as.character(OTUs_reads_cores[OTUs_reads_cores$V1==Cores_list_OTUs[[o]][4],2]))
        new_info_core<-paste(new_info_core,collapse = "")
        OTUs_reads_core[o,1]<-as.character(colnames(OTU_cores)[o])
        OTUs_reads_core[o,2]<-new_info_core
      }
    }    
  }
}

#Venn global

Reads_cores_Tree_venn<-unlist(strsplit(x = as.character(Reads_cores_Tree[,2]),split = ";")) #All reads, from the two list, are concatenated with ; and should be splitted
OTUs_reads_core_venn<-unlist(strsplit(x = as.character(OTUs_reads_core[,2]),split = ";"))
venn.diagram(list(Arbol = Reads_cores_Tree_venn, OTUs = OTUs_reads_core_venn),paste("Reads_cores_Tree_ALL","OTUs_reads_core_ALL",".png",sep = "_"),print.mode="percent",imagetype = "png") #A Venn diagram is created from the sequences names of both methods 

#Venn pairwise

#for(i in 1:dim(Reads_cores_Tree)[1])
#{
#  for(o in 1:dim(OTUs_reads_core)[1])
#  {
#    Reads_cores_Tree_venn<-strsplit(x = as.character(Reads_cores_Tree[i,2]),split = ";")[[1]]
#    OTUs_reads_core_venn<-strsplit(x = as.character(OTUs_reads_core[o,2]),split = ";")[[1]]
#    venn.diagram(list(Arbol = Reads_cores_Tree_venn, OTUs = OTUs_reads_core_venn),paste(Reads_cores_Tree[i,1],OTUs_reads_core[o,1],".png",sep = "_"), print.mode="percent",imagetype = "png")
#    
#  }
#  unlink("*.log")
#}


