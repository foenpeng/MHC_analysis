library(data.table)

setwd("/Users/pengfoen/OneDrive - University of Connecticut/MHC")

# randomly choose 60 nt sequence to run paml, >1000 would be too slow. Actually 60 nt sequence takes ~2 days.
# nt_seq<-fread("./6. R scripts to genotype fish and produce metadata/6.2 Bayesian analysis/6.2.1. Input/sequence_to_align.csv")
# nt_seq_sample<-nt_seq[sample(1:1254,60,replace=F),]
# fwrite(nt_seq_sample,"./6. R scripts to genotype fish and produce metadata/6.2 Bayesian analysis/6.2.1. Input/60_sample_nt_sequence_to_align.csv")

aas_seq <- fread("./6. R scripts to genotype fish and produce metadata/6.2 Bayesian analysis/6.2.1. Input/aas_to_align.csv")
aas_abbv<- fread("./Supertype/aas_abbreviation.csv")
aas_prop<- fread("./Supertype/aas_property.csv")
aas_prop_short<-aas_prop[aas_abbv,on=c("abbrev"="Abbreviation (3 Letter)")][!is.na(z1)]
aas_prop_short[,(4:8):=lapply(.SD, function(x){as.numeric(gsub("_","-",x))}),.SDcols=4:8]
setnames(aas_prop_short,"Abbreviation (1 Letter)","Abbrev_1")
# pos<-c(7,9,14,16,26,32,38,39,44,46,49,56) # 1st run
# pos<-c(7,14,26,32,39,44,46,49,56,65) # 2nd run
# pos<-c(7,9,14,16,26,32,38,39,44,46,49,56,65) # either significant
pos<-c(7,14,26,32,39,44,46,49,56) #both significant
aas_seq[,paste0("pos",pos):=lapply(pos,function(x) substr(aa.sequence,x,x))]

aas_matrix<-as.matrix(aas_seq[,sapply(.SD, function(x){aas_prop_short[match(x, aas_prop_short[,Abbrev_1]),.(z1,z2,z3,z4,z5)]}),.SDcols=3:11])
dim(aas_matrix)<-c(1204, 5, 9)
dimnames(aas_matrix)<-list(aas_seq[,aa.ID],colnames(aas_prop_short)[4:8],colnames(aas_seq)[3:11])

library(adegenet)
grp<-find.clusters.matrix(aas_matrix,choose.n.clust=F, criterion="diffNgroup")
grp<-find.clusters.matrix(aas_matrix,max.n.clust=120)
dapc1 <- dapc.matrix(aas_matrix, grp$grp)
