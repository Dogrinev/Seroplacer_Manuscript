#Load in GenomeTrakr Data file
GenomeTrakrData<-read.csv(file = "isolates.tsv",sep = "\t",header=FALSE,fill=TRUE,stringsAsFactors = FALSE)
#Remove top row with labels
GenomeTrakrData<-GenomeTrakrData[-1,]

setwd("~/Desktop/EC_Database_Rebuild/16s_Files")
Needs_Download<-GenomeTrakrData[,12]

#File Download Prep ===============================

library(rentrez)
for(y in Needs_Download[1:3545])
{
  assembly_search<-entrez_search(db = "assembly",term=y)
  refseq_data<-entrez_link(dbfrom = "assembly",id = assembly_search$ids,db = "nucleotide")
  RefSeq_DataFile<-entrez_fetch(db = "nuccore",id = refseq_data[[1]][1][[1]][1],rettype = "gb")
  write(RefSeq_DataFile,file=paste(y,".txt",sep = ""))
  Sys.sleep(5)
}

#Initial 16s filtering criteria ============
library(tidyverse)

setwd("~/Desktop/EC_Database_Rebuild/16s_Files")
Files_To_Parse<-list.files(path = "~/Desktop/EC_Database_Rebuild/16s_Files")
source('~/Desktop/EC_Database_Rebuild/Scripts/Parse_16s_Numbers.R')
Raw_DataTable16s<-map_dfr(Files_To_Parse, Parse_16s_Numbers)
DataTable_16s<-as.data.frame(Raw_DataTable16s[,2:4])
Numbers_16s_GoodAssemblies<-as.numeric(as.vector(DataTable_16s[,2]))
hist(Numbers_16s_GoodAssemblies)
#Table of resulting data
table(Numbers_16s_GoodAssemblies)
#3347 out of 3545 assemblies have 7 16s sequences
#=================

#First filter out the files which do not have 7 6s copy numbers
Files_With_7Copies<-Files_To_Parse[which(Numbers_16s_GoodAssemblies==7)]
New_Assemblies2<-gsub(x = Files_With_7Copies,pattern = ".txt",replacement = "")
#New_Assemblies2 is now the correct vector to use with assemblies moving forward, remaining
#assemblies should be cut down to this 3347 value.

#Need to fetch all of the FASTAs, this loop will pull the largest complete genome FASTA files based 
#on sequence length to make sure the biggest assembly is downloaded
source('~/Desktop/EC_Database_Rebuild/Scripts/Sequence_Length_Check.R')

setwd("~/Desktop/EC_Database_Rebuild/Assembly_FASTAs")

for(y in New_Assemblies2[2345:3347])
{
  assembly_search<-entrez_search(db = "assembly",term=y)
  refseq_data<-entrez_link(dbfrom = "assembly",id = assembly_search$ids,db = "nucleotide")
  Sequences_To_Scan<-refseq_data[[1]][1][[1]]
  Sequence_Length_Results<-map_int(Sequences_To_Scan, Sequence_Length_Check)
  names(Sequence_Length_Results)<-Sequences_To_Scan
  Sequence_Lengths_Sorted<-sort(Sequence_Length_Results,decreasing = TRUE)
  CDS_FASTA<-entrez_fetch(db = "nuccore",id = names(Sequence_Lengths_Sorted[1]),rettype = "fasta")
  write(CDS_FASTA,file=paste(y,".txt",sep = ""))
}

setwd("~/Desktop/EC_Database_Rebuild/CDS_FASTAs")

LastFile<-list.files()[length(list.files())]
LastFileIndex<-which(New_Assemblies2==gsub(pattern = ".txt",replacement = "",x = LastFile))
for(y in New_Assemblies2[LastFileIndex:3347])
{
  assembly_search<-entrez_search(db = "assembly",term=y)
  refseq_data<-entrez_link(dbfrom = "assembly",id = assembly_search$ids,db = "nucleotide")
  Sequences_To_Scan<-refseq_data[[1]][1][[1]]
  Sequence_Length_Results<-map_int(Sequences_To_Scan, Sequence_Length_Check)
  names(Sequence_Length_Results)<-Sequences_To_Scan
  Sequence_Lengths_Sorted<-sort(Sequence_Length_Results,decreasing = TRUE)
  CDS_FASTA<-entrez_fetch(db = "nuccore",id = names(Sequence_Lengths_Sorted[1]),rettype = "fasta_cds_na")
  write(CDS_FASTA,file=paste(y,".txt",sep = ""))
}

#7-12 night loop
setwd("~/Desktop/EC_Database_Rebuild")

#redownload GCA_013805885.1.gff3

LastFile<-list.files("~/Desktop/EC_Database_Rebuild/GFF3s")[length(list.files("~/Desktop/EC_Database_Rebuild/GFF3s"))]
LastFileIndex<-which(New_Assemblies2==gsub(pattern = ".gff3",replacement = "",x = LastFile))
for(y in New_Assemblies2[2797])
{
  assembly_search<-entrez_search(db = "assembly",term=y)
  refseq_data<-entrez_link(dbfrom = "assembly",id = assembly_search$ids,db = "nucleotide")
  Sequences_To_Scan<-refseq_data[[1]][1][[1]]
  Sequence_Length_Results<-map_int(Sequences_To_Scan, Sequence_Length_Check)
  names(Sequence_Length_Results)<-Sequences_To_Scan
  Sequence_Lengths_Sorted<-sort(Sequence_Length_Results,decreasing = TRUE)
  Nuc_Input<-names(Sequence_Lengths_Sorted[1])
  system(command = paste('wget -O GFF3s/',y,'.gff3 "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=',Nuc_Input,'"',sep = ""))
}

#Database Construction ===============

# Core Gene Determination and Filtering ============

setwd("~/Desktop/EC_Database_Rebuild/CDS_FASTAs")
library(Biostrings)
Potential_CG_List<-as.list(NULL)
for(i in 1:3347)
{
  fastaFile <- readDNAStringSet(list.files("~/Desktop/EC_Database_Rebuild/CDS_FASTAs")[i],format = "fasta")
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  SplitFASTA <- data.frame(seq_name, sequence,stringsAsFactors = FALSE)
  PotentialCGs<-str_extract(string = SplitFASTA[,1],pattern = "(?<=\\[gene=)[:alpha:]*")
  Potential_CG_List[[i]]<-PotentialCGs
}
#Make collapsed version for string searching in a further loop
Potential_CG_List_Collapsed<-map(.x = Potential_CG_List,.f = ~paste(.,collapse=""))
#Pick any large assembly as the list of tentative core genes, does not matter which is picked
#since in the end we will filter down to core genes which are in most.
Potential_Core_Genes<-Potential_CG_List[[1]]
Potential_Core_Genes2<-na.omit(Potential_Core_Genes)

CoreGeneList<-as.list(NULL)

for(i in Potential_Core_Genes2)
{
  True_Value<-table(map_chr(.x = Potential_CG_List_Collapsed,.f = ~str_detect(string = .,pattern = i)))[2]
  CoreGeneList[[i]]<-True_Value
}

CoreGeneList2<-unlist(CoreGeneList)

#The following table shows how many assemblies does each core gene appear in, we are assuming that
#a core gene appears in most if not all assemblies. We want to use as many core genes as 
#possible even if it requires cutting some assemblies. The first filtering step is to select a
#primary core gene pool. 
table(CoreGeneList2)
#For this version of database construction, we cut at core genes represented in >3346 assemblies
CoreGeneList3<-CoreGeneList2[which(CoreGeneList2>3345)]
EC_Core_Genes_Pre<-names(CoreGeneList3)
EC_Core_Genes<-gsub(pattern = ".TRUE",replacement = "",x = EC_Core_Genes_Pre)
#core genes list is set now
###########################

source('~/Desktop/EC_Database_Rebuild/Scripts/Core_Gene_Extractor.R')
EC_Core_Genes_Adj<-paste("gene=",EC_Core_Genes,sep="")
Core_Genes_List<-as.list(NULL)

for(i in 1:3347)
{
  fastaFile <- readDNAStringSet(list.files("~/Desktop/EC_Database_Rebuild/CDS_FASTAs")[i],format = "fasta")
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  SplitFASTA <- data.frame(seq_name, sequence,stringsAsFactors = FALSE)
  Core_Genes_Table<-map_dfr(EC_Core_Genes_Adj,Core_Gene_Extractor)
  Core_Genes_Table_Final<-Core_Genes_Table[!duplicated(Core_Genes_Table$seq_name),]
  FixedGeneNames<-str_extract(string = Core_Genes_Table_Final[,1],pattern = "(?<=\\[gene=)[:alpha:]*")
  AssemblyName<-paste(list.files("~/Desktop/EC_Database_Rebuild/CDS_FASTAs")[i])
  AssemblyNameFixed<-gsub(pattern = ".txt",replacement = "",x = AssemblyName)
  Core_Genes_Table_Final[,1]<-paste(AssemblyNameFixed,FixedGeneNames,sep = " ")
  Core_Genes_Table_Final2<-Core_Genes_Table_Final[!(duplicated(Core_Genes_Table_Final[,1]) | duplicated(Core_Genes_Table_Final[,1], fromLast = TRUE)), ]
  Core_Genes_List[[i]]<-Core_Genes_Table_Final2
}

setwd("~/Desktop/EC_Database_Rebuild")
save(Core_Genes_List,file = "Core_Genes_List_EC.Rdata")
Core_Genes_List_Cut<-Core_Genes_List

#Core Gene Assemblies is a list containing all potential core genes with which assemblies have those genes,
#next we have to filter appropriately
source('~/Desktop/EC_Database_Rebuild/Scripts/CoreGeneBuilder.R')
Core_Gene_Assemblies<-map(.x = EC_Core_Genes, .f = ~ Core_Gene_Builder(.))
save(Core_Gene_Assemblies,file = "Core_Gene_Assemblies.Rdata")
Core_Gene_Assemblies_NoZero<-map_int(.x = Core_Gene_Assemblies,.f = ~ length(rownames(.)))

#2 options here:
#The following calculation is reached by doing a manual adjustment of the cut range. A minimum is set
#at what level Core_Gene_Assemblies_NoZero should be trimmed, this trimming is completely flexible
#and up to the user during database construction. It is reccomended to test a variety of ranges
#and determine what is best. Higher numbers will result in fewer core genes but more assemblies 
#preserved, while lowering this number will expand the pool of core genes but cut out more assemblies.

#EDIT THESE FOR EC
#The bottom line will result in 92 core genes and 1703 assemblies remaining in the dataset.
#The top line will result in 349 core genes and 1681 assemblies remaining in the dataset.
#For this build, we will focus on 92 core genes and retain as many assemblies as possible, so we will
#choose to execute the bottom line:

#The remaining assembly number is calculated in the following lines by evaluating table(Assemblies_To_Keep)
Core_Gene_Assemblies_Cut<-Core_Gene_Assemblies[which(Core_Gene_Assemblies_NoZero>3345 & Core_Gene_Assemblies_NoZero<3347)]

source('~/Desktop/EC_Database_Rebuild/Scripts/Assemblies_Missing_Genes.R')
Assemblies_To_Keep<-map_int(.x = New_Assemblies2,.f = Assemblies_Missing_Genes)
table(Assemblies_To_Keep)
names(Assemblies_To_Keep)<-New_Assemblies2
Assemblies_With_CGs_Index<-which(Assemblies_To_Keep==101)
#3336 total assemblies remain

Core_Genes_List_Cut<-Core_Genes_List[Assemblies_With_CGs_Index]

Core_Gene_Assemblies2<-map(.x = EC_Core_Genes, .f = ~ Core_Gene_Builder(.))
Core_Gene_Assemblies_NoZero2<-map_int(.x = Core_Gene_Assemblies2,.f = ~ length(rownames(.)))

#second version to focus on final stuff 
Core_Gene_Assemblies_Cut2<-Core_Gene_Assemblies2[which(Core_Gene_Assemblies_NoZero2==3336)]
source('~/Desktop/EC_Database_Rebuild/Scripts/FixColNames.R')
Core_Gene_Assemblies_Cut3<-map(.x = Core_Gene_Assemblies_Cut2,.f = FixColNames)

Core_Gene_Assemblies_Cut4<-Core_Gene_Assemblies_Cut3


##### NEXT SECTION

library(phylotools)
setwd("~/Desktop/EC_Database_Rebuild/Alignment_FASTAs")
for(i in 1:101)
{
  Data_Table<-Core_Gene_Assemblies_Cut4[[i]]
  dat2fasta(dat = Data_Table,outfile = paste(i,".fasta",sep = ""))
}

#generating mafft alignments using command line
for(i in 1:101)
{
  system(command = paste("/usr/local/bin/mafft --retree 2 --inputorder /Users/dogrinev/Desktop/EC_Database_Rebuild/Alignment_FASTAs/",i,".fasta > /Users/dogrinev/Desktop/EC_Database_Rebuild/MAFFT_Outputs/",i,"_aligned.fasta",sep = ""))
}

MAFFT_List<-as.list(NULL)

for (i in 1:101)
{
  fastaFile <- readDNAStringSet(paste("~/Desktop/EC_Database_Rebuild/MAFFT_Outputs/",i,"_aligned.fasta",sep=""),format = "fasta")
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  Aligned_FASTA <- data.frame(seq_name, sequence,stringsAsFactors = FALSE)
  MAFFT_List[[i]]<-Aligned_FASTA
}

StringMaker<-function(xx)
{
  FullString<-map_chr(.x = MAFFT_List,.f = ~paste(.[xx,2]))
  FullString2<-paste(FullString,sep="",collapse="")
  return(data.frame(FullString2))
}

Test<-1:3336
Maybe<-map_dfr(.x = Test,.f = StringMaker)
Test2<-str_extract(string = seq_name,pattern = "GCA_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].[0-9]")
Test3<-gsub(pattern = ".",replacement = "",x = Test2,fixed = TRUE)

MaybeFinal<-cbind(Test3,Maybe)
colnames(MaybeFinal)=c("seq.name","seq.text")
library(phylotools)
dat2fasta(dat = MaybeFinal,outfile = "RAXML_Input.fasta")

#Getting into constructing the 16s parsing table
#Use Assembly_Names_Subset fo downloads
#Use Limited_FASTAs to test

#first need to rename headers of raw fasta files
library(Biostrings)
library(seqinr)
Tree_Assemblies<-Test2

setwd("~/Desktop/EC_Database_Rebuild")
for(y in Tree_Assemblies[2:3336])
{
  FASTA_To_Rename<-readDNAStringSet(filepath = paste("Assembly_FASTAs/",y,".txt",sep=""),format="fasta")
  names(FASTA_To_Rename)<-"chr"
  write.fasta(sequences = paste(FASTA_To_Rename),names = names(FASTA_To_Rename),file.out = paste("Raw_FASTAs_CHR/",y,".txt",sep=""))
}

#now looping to generate BED files
library(stringr)
library(tidyverse)

setwd("~/Desktop/EC_Database_Rebuild")

for(y in Tree_Assemblies[2787:2796])
{
  Test<-read.csv(file = paste("GFF3s/",y,".gff3",sep=""),sep = "\r",header = FALSE)
  Test2<-as.data.frame(Test[which(str_detect(string = Test$V1,pattern = "16S ribosomal RNA")),])
  Test3<-as.data.frame(Test2[c(1,3,5,7,9,11,13),])
  Position_Table<-map_dfc(.x = Test3[,1],.f = ~str_split(string = .,pattern = "\t"))
  Position_Table_Cut<-t(Position_Table[4:5,1:7])
  First_Column<-rep("chr",7)
  Last_Column<-c("16s_1","16s_2","16s_3","16s_4","16s_5","16s_6","16s_7")
  Final_Position_Table<-as.data.frame(cbind(First_Column,Position_Table_Cut,Last_Column))
  write.table(x = Final_Position_Table,file = paste("GFF3s_Positions/",y,".bed",sep=""),sep="\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
}

#bedtools section to make 16s data table

setwd("~/Desktop/EC_Database_Rebuild")

for(y in Tree_Assemblies)
{
  system(command = paste("bedtools getfasta -fi Raw_FASTAs_CHR/",y,".txt -bed GFF3s_Positions/",y,".bed -fo 16s_Sequences_BT/",y,".txt",sep=""))
  system(command = paste("rm Raw_FASTAs_CHR/",y,".txt.fai",sep=""))
}

#script to merge 16s files into one data frame
setwd("~/Desktop/EC_Database_Rebuild/16s_Sequences_BT")
source('~/Desktop/EC_Database_Rebuild/Scripts/Merger_16s.R')
Full_16s_Data<-map(Tree_Assemblies, Merger_16s)

#need to fix bad nucleotides
source('~/Desktop/EC_Database_Rebuild/Scripts/GSubber.R')
Full_16s_Data_Fixed<-map(.x = Full_16s_Data,.f = ~map_chr(.x = .,.f = GSubber))

#Finding index values for strands with - that need to be RCed
source('~/Desktop/EC_Database_Rebuild/Scripts/Strand_Searcher.R')
setwd("~/Desktop/EC_Database_Rebuild")
Strand_Indices<-map(.x = Tree_Assemblies,.f = Strand_Searcher)
#actual functions that do the reverse complement fixing - input range equal to number of files
source('~/Desktop/EC_Database_Rebuild/Scripts/ReverseComplementFixer.R')
source('~/Desktop/EC_Database_Rebuild/Scripts/ReverseComplementer.R')

Full_16s_Data_RC<-map(.x = c(1:3336),.f = ReverseComplementFixer)
#List number 1097 had to be removed because the FASTA disagrees with the GFF file, the final
#allele positions are outside of the range of the total file causing the 7th 16S allele to 
#not get pulled. This entry will be removed because it is missing a 16S allele.
#We need to check if any others are missing 16S alleles:
Bad_Assembly_Check<-map(.x = 1:3336,.f = ~length(Full_16s_Data_RC[[.]]))
Bad_Assembly_Check2<-unlist(Bad_Assembly_Check)
table(Bad_Assembly_Check2)
Remove<-which(Bad_Assembly_Check2<8)
#Entry 500 and 732 need to be removed because they are missing an allele
Full_16s_Data_RC2<-Full_16s_Data_RC[-Remove]

save(Full_16s_Data_RC2,file = "Full_16s_Data_RC.Rdata")

#aligned maker
Prep_For_MSA<-as.data.frame(unlist(Full_16s_Data_RC2))

Set1<-seq.int(from = 2,to = 26666,by = 8)
Set2<-seq.int(from = 3,to = 26667,by = 8)
Set3<-seq.int(from = 4,to = 26668,by = 8)
Set4<-seq.int(from = 5,to = 26669,by = 8)
Set5<-seq.int(from = 6,to = 26670,by = 8)
Set6<-seq.int(from = 7,to = 26671,by = 8)
Set7<-seq.int(from = 8,to = 26672,by = 8)

Alleles1<-as.data.frame(Prep_For_MSA[Set1,])
Alleles2<-as.data.frame(Prep_For_MSA[Set2,])
Alleles3<-as.data.frame(Prep_For_MSA[Set3,])
Alleles4<-as.data.frame(Prep_For_MSA[Set4,])
Alleles5<-as.data.frame(Prep_For_MSA[Set5,])
Alleles6<-as.data.frame(Prep_For_MSA[Set6,])
Alleles7<-as.data.frame(Prep_For_MSA[Set7,])

Tree_Assemblies_Prep<-map(.x = 1:3334,.f = ~Full_16s_Data_RC2[[.]][1])
Tree_Assemblies2<-unlist(Tree_Assemblies_Prep)

Alleles1_MSA<-cbind(Tree_Assemblies2,Alleles1)
colnames(Alleles1_MSA)<-c("seq.name","seq.text")
Alleles2_MSA<-cbind(Tree_Assemblies2,Alleles2)
colnames(Alleles2_MSA)<-c("seq.name","seq.text")
Alleles3_MSA<-cbind(Tree_Assemblies2,Alleles3)
colnames(Alleles3_MSA)<-c("seq.name","seq.text")
Alleles4_MSA<-cbind(Tree_Assemblies2,Alleles4)
colnames(Alleles4_MSA)<-c("seq.name","seq.text")
Alleles5_MSA<-cbind(Tree_Assemblies2,Alleles5)
colnames(Alleles5_MSA)<-c("seq.name","seq.text")
Alleles6_MSA<-cbind(Tree_Assemblies2,Alleles6)
colnames(Alleles6_MSA)<-c("seq.name","seq.text")
Alleles7_MSA<-cbind(Tree_Assemblies2,Alleles7)
colnames(Alleles7_MSA)<-c("seq.name","seq.text")

Alleles_MSA<-rbind(Alleles1_MSA,Alleles2_MSA,Alleles3_MSA,Alleles4_MSA,Alleles5_MSA,Alleles6_MSA,Alleles7_MSA)

setwd("~/Desktop/EC_Database_Rebuild/MSA")

dat2fasta(dat = Alleles_MSA,outfile = "Alleles_MSA.fasta")

#Alignment Command:
#mafft --auto --inputorder Full_Alignment.txt > Full_Alignment_Finished.txt

#first reload the aligned 16s files
fastaFile <- readDNAStringSet("~/Desktop/EC_Database_Rebuild/MSA/Full_Alignment_EC.txt",format = "fasta")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
Alleles_Aligned <- data.frame(seq_name, sequence,stringsAsFactors = FALSE)

Unset1<-seq.int(from = 1,to = 3334,by = 1)
Unset2<-seq.int(from = 3335,to = 6668,by = 1)
Unset3<-seq.int(from = 6669,to = 10002,by = 1)
Unset4<-seq.int(from = 10003,to = 13336,by = 1)
Unset5<-seq.int(from = 13337,to = 16670,by = 1)
Unset6<-seq.int(from = 16671,to = 20004,by = 1)
Unset7<-seq.int(from = 20005,to = 23338,by = 1)

Full_16s_Data_Aligned<-map(.x = 1:3334,.f = ~c(Tree_Assemblies2[.],Alleles_Aligned[Unset1[.],2],Alleles_Aligned[Unset2[.],2],Alleles_Aligned[Unset3[.],2],Alleles_Aligned[Unset4[.],2],Alleles_Aligned[Unset5[.],2],Alleles_Aligned[Unset6[.],2],Alleles_Aligned[Unset7[.],2]))

save(Core_Gene_Assemblies,Core_Gene_Assemblies_Cut2,Core_Gene_Assemblies_Cut3,Core_Gene_Assemblies2,Core_Genes_List,Core_Genes_List_Cut,Potential_CG_List,file = "Large_Files_Backup.Rdata")

RAXML_Output<-readChar("RAXML_Input.fasta.raxml.bestTree",nchars = 1e7)
vert.tree<-read.tree(text=RAXML_Output)

Tree_Assemblies_Cut<-gsub(pattern = ".",replacement = "",x = Tree_Assemblies,fixed=TRUE)

TipLabel<-vert.tree$tip.label
TipLabelsFixed<-Tree_Assemblies[match(TipLabel,Tree_Assemblies_Cut)]
vert.tree.2<-vert.tree
vert.tree.2$tip.label<-TipLabelsFixed

#Cut the 2 tips which were removed during 16S database construction
vert.tree.cut<-drop.tip(phy = vert.tree.2,tip = which(vert.tree.2$tip.label %in% Tree_Assemblies2 == FALSE))
vert.tree.EC<-vert.tree.cut

save(vert.tree.EC,Tree_Assemblies2,Full_16s_Data_Aligned,Full_16s_Data_RC2,file = "AnalysisFiles_1.Rdata")
