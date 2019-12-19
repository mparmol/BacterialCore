#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

setwd(paste("./",args[1],sep = ''))

if (!require("biomformat")){ #This package is necessary to read the .biom files
  source("http://bioconductor.org/biocLite.R")
  biocLite("biomformat")
}
if (!require("picante")) install.packages("picante", repos='http://cran.us.r-project.org')
if (!require("matrixStats")) install.packages("matrixStats", repos='http://cran.us.r-project.org')

All<-biom_data(read_biom(paste('./0.97/All_table_json_0.97.biom'))) #First of all we read the .biom table at 0.97 level to get the subjects names
All<-All[,order(colnames(All))] #We order the subjects names
Counts<-colSums(as.array(All)) #We get a sum of each column
Cores<-data.frame(matrix(ncol=length(All[1,]))) #Create a data frame to save de satitistical results
colnames(Cores)=colnames(All) #The colnames of the new data frame are the same as the original data frame

total_files<-list.files(path=".",full.names = TRUE, recursive = TRUE) #We get a list of all files containd into the OTUs folder

for (i in 1:length(total_files)) #With this loop we get the abundance information from each core for each subject
{
  
  if(grepl("changed_name",total_files[i]) & file.info(total_files[i])$size > 0)
  {
    Biom_table<-read_biom(total_files[i])
    #print(i)
    Row_saved<-rownames(Biom_table)
    Biom_table2<-t(as.matrix(biom_data(Biom_table)))
    Biom_table2<-t(as.matrix(Biom_table2[,order(colnames(Biom_table2))]))
    if(dim(Biom_table2)[1]==1)
    {
      rownames(Biom_table2)<-Row_saved
    }
    
    Cores<-rbind(data.frame(list(Biom_table2)),Cores) #We save the info in the Cores data frame
  }
}

erase_row=length(Cores[,1]) #We erase the last results row
Cores=Cores[-erase_row,]
Cores<-Cores[,order(colnames(Cores))]

#colnames(Biom94)==colnames(All)

#esto era para TWINS, variará un poco en función del dataset, camabiarlo par que se pueda estandarizar

#Cores<-do.call("rbind",list(Biom94,Biom92,as.matrix(Biom91),as.matrix(Biom90),Biom89,Biom88,as.matrix(Biom87),as.matrix(Biom85),Biom80)) #estandarizar mejor, creo que los as.matrix era para los Leveles que tenían más de uno, quizas se pueda aplicar a todos y ya está
Cores.freq<-t(t(Cores)/Counts) #Here we calculate the frequency for each subject in each core obtained
Cores.freq.means<-rowMeans(Cores.freq) #Mean
Cores.freq.Sds<-rowSds(Cores.freq) #Standard deviation
Cores.freq.Stats<-rbind(Cores.freq.means,Cores.freq.Sds) #An auxiliary data frame to join means and standar deviation
Cores.freq.quant<-rowQuantiles(Cores.freq) #Quantiles calculation for each subject
Cores.freq.CV<-Cores.freq.Sds/Cores.freq.means #Coefficient of variation
Cores.freq.Stats<-rbind(Cores.freq.means,Cores.freq.Sds,Cores.freq.CV)
#write.table(Cores.freq.Stats,'Core_twins_Stats2.txt',sep='\t',quote=F)

###############################################################################################

results_ecdf=data.frame(matrix(ncol=2,nrow = length(Cores.freq.means)))
colnames(results_ecdf)=c("ecdf-AV","ecdf-CV")
u=1

for (i in 1:length(total_files))
{
  
  if(grepl("json",total_files[i])) #Inside this loop we look for the tables that contains all the abundance OTUs information of all subjects 
  {
    Biom_table<-biom_data(read_biom(total_files[i])) #We calculate the same variables as we did for the biom cores tables
    Biom_table.freq<-t(t(as.matrix(Biom_table))/Counts)
    Biom_table.freq.means<-rowMeans(Biom_table.freq)
    Biom_table.freq.Sds<-rowSds(Biom_table.freq)
    Biom_table.freq.CV<-Biom_table.freq.Sds/Biom_table.freq.means
    Biom_table.f<-ecdf(Biom_table.freq.means) #The ecdf (empirical distribution function) calculation is performed here (the position of the mean)
    Biom_table.fcv<-ecdf(Biom_table.freq.CV) #The same but with the coefficient of variation
    
    name2=strsplit(total_files[i],'.biom')[[1]][1]
    name2=strsplit(name2,'_')[[1]][4]
    #name2=as.character(0.97)
    
    for(o in 1:length(Cores.freq.means)) #Inside this loop we merge the total info obtained from the original abundance table and the info obtained for the cores in each level
    {
      name1=strsplit(names(Cores.freq.means[o])," ")[[1]][2]
      #name1=strsplit(name1,'._')[[1]][1]
      #print(total_files[i])
      #name2=as.character(0.97)
      #print(name1[[1]][1])
      
      if(name1==name2)
      {
        #print(names(Cores.freq.means[o]))
        denovos<-c(Biom_table.f(Cores.freq.means[o]),Biom_table.fcv(Cores.freq.CV[o])) #What we need for each core is are the mean and the frequency
        #print(denovos)  
        
        results_ecdf[u,1:2]=denovos
        rownames(results_ecdf)[u]=names(Cores.freq.means[o])
        #print(rownames(results_ecdf)[u])
        u=u+1
        
      }
      
    }
  }
}

results_ecdf2<-results_ecdf #Auxiliary variable to avoid the overwrite of the original results_ecdf

#rownames(results_ecdf2)<-results_ecdf[,1]



if(dim(results_ecdf2)[1]==1)
{
  results_ecdf2<-data.frame(t(results_ecdf2))
  results_ecdf3<-rbind(data.frame(Cores.freq.Stats),data.frame(results_ecdf2))
  
}else{
  results_ecdf2<-t(results_ecdf2)
  results_ecdf2<-results_ecdf2[,order(colnames(results_ecdf2))] 
  Cores.freq.Stats<-Cores.freq.Stats[,order(colnames(Cores.freq.Stats))]
  results_ecdf3<-rbind(Cores.freq.Stats,results_ecdf2) #We order the results of the mean, coefficient of variation and ecdf for each core at each level an bind them to a new data frame, printed as results.
}

write.table(results_ecdf3,'Distr.txt',sep='\t',quote=F) ###Statistical data for each core

#Cores.freq.final=rbind(Cores.freq,colSums(Cores.freq))
#rownames(Cores.freq.final)[dim(Cores.freq.final)[1]]="Sum"
Cores.freq.final=Cores.freq
write.table(Cores.freq.final,'Core_Otus_freq.txt',sep='\t',quote=F) ###Cores and abundance per subject

if(dim(results_ecdf2)[2]==1)
{
  core_list=unlist(strsplit(colnames(results_ecdf3), "\\."))
  core_list[2]=paste(core_list[2],core_list[3],sep = ".")
  core_list=list(core_list)
}else
{
  core_list=strsplit(colnames(results_ecdf3)," ")
}
  
OTUs_total <- data.frame(matrix(ncol=2))

for (i in 1:length(core_list))
{
  OTUs_total[i]=core_list[[i]][2]
}

OTUs_total=data.frame(table(as.numeric(OTUs_total)))
OTUs_total$Var1=as.character(OTUs_total$Var1)
OTUs_total=rbind(OTUs_total,c("Total",sum(OTUs_total$Freq)))
colnames(OTUs_total)[1]="Level"
colnames(OTUs_total)[2]="Core groups"
write.table(OTUs_total,'General.txt',sep='\t',quote=F,row.names = FALSE) ###Total number of cores per level

OTUs_total=OTUs_total[dim(OTUs_total)[1],]
OTUs_total[dim(OTUs_total)[1],3]=mean(colSums(Cores.freq))
OTUs_total[dim(OTUs_total)[1],4]=sd(colSums(Cores.freq))
OTUs_total[dim(OTUs_total)[1],5]=min(colSums(Cores.freq))
OTUs_total[dim(OTUs_total)[1],6]=max(colSums(Cores.freq))
colnames(OTUs_total)[3:6]=c("Mean","SD","Min","Max")
write.table(OTUs_total,'Statistics.txt',sep='\t',quote=F,row.names = FALSE) ###Phylogenetic total core statistics, for all cores

######################### Control permutations. In this module we perform 1000 random permutations for each abundance biom table and calculate the number of cores we obtain, at each level, to get a negative control

Cores<-data.frame(matrix(ncol=100,nrow = 1))

for (i in 1:length(total_files))
{
  if(grepl("Filt_",total_files[i]))
  {
    #print(total_files[i])
    All<-as.matrix(biom_data(read_biom(total_files[i])))
    
    Vec<-data.frame(matrix(ncol=100,nrow = 1))
    
    for (o in 1:100) 
    {
      Simu<-randomizeMatrix(All,null.model='frequency',iterations=1000)
      Simu[Simu>0]<-1
      Simu.rowsums<-rowSums(as.matrix(Simu))
      Vec[,o]<-length(Simu.rowsums[which(as.numeric(Simu.rowsums)==dim(All)[2])])
    }
    #a<-sum(Vec)
    
    name2=strsplit(total_files[i],'.Filt')[[1]][1]
    name2=strsplit(name2,'/')[[1]][2]
    
    rownames(Vec)<-name2
    
    Cores<-rbind(Vec,Cores)
    
    #print (name2)
  }
}  

erase_row=length(Cores[,1])
Cores=Cores[-erase_row,]
Cores_sumado=Cores
Cores_sumado=data.frame(Cores_sumado[,1])
rownames(Cores_sumado)=rownames(Cores)
colnames(Cores_sumado)="Number of positives (1000 perm)"
Cores_sumado[,1]=rowSums(Cores)

write.table(Cores_sumado,paste("control_random_table_cores.txt",sep = "_"),sep = "\t")

#warnings()
