##read files from gtf

setwd('/Users/rui/Google Drive/Rui/1 data analysis/Annotation_WT26 from Micheal/')


gff<-read.table('4c_anno_gencodeOverlap.gtf',sep='\t',header=F)
#colnames(gff)<-c('Seqname',	'Source',	'Feature',	'Start',	'End',	'Score',	'Strand',	'Frame',	'Group')


#detect duplicates combining 'Seqname',	'Start',	'End',	'Strand' 
duplicates<-duplicated(gff[,c(1,4,5,7)])

#remove duplicates  
de_gff<-gff[!duplicates,]

#Use 'rtracklayer' to save into gtf format
#source("http://bioconductor.org/biocLite.R")
#biocLite('rtracklayer')
rtracklayer::export(de_gff,'de_duplicates_4c_anno.gtf')
