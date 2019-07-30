library(data.table)
library(glmnet)


setwd("/Users/pengfoen/OneDrive - University of Connecticut/MHC")
MHC_protein <- fread("./7. Metadata and MHC data combined/7.0.BayesProteinsMergedMetadata.csv")
MHC_gene <- fread("./7. Metadata and MHC data combined/7.0.BayesGenotypesMergedMetadata.csv")
head(MHC_protein)

MHC_name_all<-colnames(MHC_protein)[7:1210]
Parasite_name_all<-colnames(MHC_protein)[1271:1311]

### clean parasite data
# count how many times each parasite appear in the table, there are many unfrequent ones
Parasite_count<-MHC_protein[,lapply(.SD, function(x) sum(x>0,na.rm=T)),.SDcols=Parasite_name_all]
Parasite_count<-data.table(cn = names(Parasite_count), transpose(Parasite_count))
Parasite_name<-Parasite_count[V1>27,cn] # require parasite appear more than 27 times, the number of lakes

# remove rows with no parasite info and lakes with too few samples
MHC_protein_cleaned<-na.omit(MHC_protein, cols=c("Parasite.richness"))
anyNA(MHC_protein_cleaned[,..Parasite_name])
MHC_protein_cleaned<-MHC_protein_cleaned[,if(.N>1) .SD, by="site_name"]

### clean MHC data
# remove MHC with 0 appearance
MHC_count<-MHC_protein_cleaned[,lapply(.SD, function(x) sum(x>0,na.rm=T)),.SDcols=MHC_name_all]
MHC_count<-data.table(cn = names(MHC_count), transpose(MHC_count))
MHC_name<-MHC_count[V1>0,cn]

# look for highly correlated MHC and delete one of them
temp_MHC<-MHC_protein_cleaned[,..MHC_name]
temp_cor<-cor(temp_MHC)
temp_cor[!lower.tri(temp_cor)] <- 0
length(colnames(temp_MHC)[apply(temp_cor,2,function(x) any(abs(x) > 0.99))])
MHC_name<-colnames(temp_MHC)[!apply(temp_cor,2,function(x) any(abs(x) > 0.99))]

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


