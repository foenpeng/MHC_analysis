library(data.table)
library(glmnet)
library(Hmisc)
library(corrplot)

setwd("/Users/pengfoen/OneDrive - University of Connecticut/MHC")
MHC_protein <- fread("./7. Metadata and MHC data combined/7.0.BayesProteinsMergedMetadata.csv")
MHC_gene <- fread("./7. Metadata and MHC data combined/7.0.BayesGenotypesMergedMetadata.csv")
head(MHC_protein)

MHC_name_all<-colnames(MHC_protein)[7:1210]
Parasite_name_all<-colnames(MHC_protein)[1271:1311]
Lake_name_all<-unlist(unique(MHC_protein[,"site_name"]))

### clean parasite data
# remove rows with no parasite info and lakes with too few samples
MHC_protein_cleaned<-na.omit(MHC_protein, cols=c("Parasite.richness"))
anyNA(MHC_protein_cleaned[,..Parasite_name_all])
MHC_protein_cleaned<-MHC_protein_cleaned[,if(.N>1) .SD, by="site_name"]
Lake_name<-unlist(unique(MHC_protein_cleaned[,"site_name"]))

{ 
  # old clean
# count how many times each parasite appear in the table, there are many unfrequent ones
Parasite_count_bylake<-MHC_protein[,lapply(.SD, function(x) sum(x>0,na.rm=T)),.SDcols=Parasite_name_all, by="site_name"][,rowSums(.SD, na.rm = TRUE),.SDcols=Parasite_name_all, by="site_name"]
Parasite_appearance_bylake<-MHC_protein[,c(.N,lapply(.SD, function(x) sum(x>0,na.rm=T))),.SDcols=Parasite_name_all, by="site_name"][,lapply(.SD, function(x) sum(x>0,na.rm=T)),.SDcols=Parasite_name_all]
Parasite_appearance_bylake<-data.table(parasite = names(Parasite_appearance_bylake), transpose(Parasite_appearance_bylake))
#Parasite_name<-Parasite_appearance_bylake[V1>10,parasite] # keep parasite that appears in at least 11 lakes
Parasite_name<-Parasite_name_all

Parasite_count<-MHC_protein[,lapply(.SD, function(x) sum(x>0,na.rm=T)),.SDcols=Parasite_name]
Parasite_count<-data.table(cn = names(Parasite_count), transpose(Parasite_count))



### clean MHC data
MHC_appearance_bylake<-MHC_protein_cleaned[,c(.N,lapply(.SD, function(x) sum(x>0,na.rm=T))),.SDcols=MHC_name_all, by="site_name"][,lapply(.SD, function(x) sum(x>0,na.rm=T)),.SDcols=MHC_name_all]
MHC_appearance_bylake<-data.table(MHC = names(MHC_appearance_bylake), transpose(MHC_appearance_bylake))
#MHC_name<-MHC_appearance_bylake[V1>1,MHC] # require MHC appears in at least two lakes
MHC_name<-MHC_name_all

# remove MHC with 0 appearance
MHC_count<-MHC_protein_cleaned[,lapply(.SD, function(x) sum(x>0,na.rm=T)),.SDcols=MHC_name]
MHC_count<-data.table(cn = names(MHC_count), transpose(MHC_count))

# look for highly correlated MHC and delete one of them
temp_MHC<-MHC_protein_cleaned[,..MHC_name]
temp_cor<-cor(temp_MHC)
temp_cor[!lower.tri(temp_cor)] <- 0
length(colnames(temp_MHC)[apply(temp_cor,2,function(x) any(abs(x) > 0.99))])
MHC_name<-colnames(temp_MHC)[!apply(temp_cor,2,function(x) any(abs(x) > 0.99))]}

{
  # new clean 8/8/2019
  
  # choose parasite in each lake that has frequency ranging from 0.1-0.9
  Parasite_prevalence_bylake<-MHC_protein_cleaned[,lapply(.SD, function(x) sum(x>0,na.rm=T)/(.N)), by='site_name',.SDcols=Parasite_name_all]
  Parasite_prevalence_bylake<-data.table(Parasite_name = names(Parasite_prevalence_bylake), transpose(Parasite_prevalence_bylake))
  colnames(Parasite_prevalence_bylake)<-as.character(unlist(Parasite_prevalence_bylake[1,]))
  Parasite_prevalence_bylake<-Parasite_prevalence_bylake[-1,]
  Parasite_prevalence_bylake[,lapply(.SD, function(x) sum(x>0.05 & x <0.95)),.SDcols=2:27]
  setnames(Parasite_prevalence_bylake,"site_name","Parasite_name")
  
  # choose MHC in each lake that has frequency ranging from 0.1-0.9
  MHC_prevalence_bylake<-MHC_protein_cleaned[,lapply(.SD, function(x) sum(x>0,na.rm=T)/(.N)), by='site_name',.SDcols=MHC_name_all]
  MHC_prevalence_bylake<-data.table(MHC_name = names(MHC_prevalence_bylake), transpose(MHC_prevalence_bylake))
  colnames(MHC_prevalence_bylake)<-as.character(unlist(MHC_prevalence_bylake[1,]))
  MHC_prevalence_bylake<-MHC_prevalence_bylake[-1,]
  MHC_prevalence_bylake[,lapply(.SD, function(x) sum(x>0.05 & x <0.95)),.SDcols=2:27]
  setnames(MHC_prevalence_bylake,"site_name","MHC_name")
}

############ Regression ############

# Simple glm
MHC_protein_subset<-MHC_protein_cleaned[,c(MHC_name,Parasite_name),with=F]
# fit glm for all lakes
fits <- lapply(Parasite_name,
               function(y)glm(paste(y , paste("~", (paste(MHC_name, collapse = " + ")))),data=MHC_protein_subset,family = gaussian()))


# GLM with regularization and sparse matrix

Parasite_name_bi = paste("bi", Parasite_name, sep = ".")
MHC_protein_cleaned[,c(Parasite_name_bi) := lapply(.SD, function(x) as.numeric(x>0,na.rm=T)),.SDcols=Parasite_name]

y_train <- as.matrix(MHC_protein_cleaned[,bi.Unionidae_internal])
x_train <- sparse.model.matrix(~ . -1, data=MHC_protein_cleaned[,..MHC_name])

# For example for logistic regression using L1 norm (lasso) 
cv.fit <- cv.glmnet(x=x_train, y=y_train, family='binomial', alpha=1, nfolds=5, parallel=TRUE)

plot(cv.fit)

# lm with single MHC as dependent variable
res_cor<-list()
MHC_name_lake<-list()
Parasite_name_lake<-list()

for(lake in Lake_name){
  MHC_name_lake[[lake]]<-unlist(MHC_prevalence_bylake[get(lake)>0.05 & get(lake)<0.95,MHC_name])
  Parasite_name_lake[[lake]]<-unlist(Parasite_prevalence_bylake[get(lake)>0.05 & get(lake)<0.95,Parasite_name])
  MHC_cor_temp<-MHC_protein_cleaned[site_name==lake,c(MHC_name_lake$lake,Parasite_name_lake$lake),with=F]
  cor_temp<-rcorr(as.matrix(MHC_cor_temp))
  res_cor[[lake]]<-cor_temp
  corrplot(cor_temp$r[MHC_name_temp,Parasite_name_temp], type="full", order="original", 
           p.mat = cor_temp$P[MHC_name_temp,Parasite_name_temp], sig.level = 0.05, insig = "blank",
           cl.lim=c(-1,1), col=colorRampPalette(c("blue","white","red"))(200))
}

names(sort(table(unlist(Parasite_name_lake)),decreasing=TRUE)[1:10])
sort(table(unlist(MHC_name_lake)),decreasing=TRUE)[1:10]


# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
  
