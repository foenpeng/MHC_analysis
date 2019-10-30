library(data.table)
library(glmnet)
library(Hmisc)
library(corrplot)

setwd("/Users/pengfoen/OneDrive - University of Connecticut/MHC")
MHC_protein <- fread("./7. Metadata and MHC data combined/7.0.BayesProteinsMergedMetadata.csv")
MHC_gene <- fread("./7. Metadata and MHC data combined/7.0.BayesGenotypesMergedMetadata.csv")

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
MHC_name<-colnames(temp_MHC)[!apply(temp_cor,2,function(x) any(abs(x) > 0.99))]
}

{
  # new clean 8/8/2019
  
  # choose parasite in each lake that has frequency ranging from 0.05-0.95
  Parasite_prevalence_bylake<-MHC_protein_cleaned[,lapply(.SD, function(x) sum(x>0,na.rm=T)/(.N)), by='site_name',.SDcols=Parasite_name_all]
  Parasite_prevalence_bylake<-data.table(Parasite_name = names(Parasite_prevalence_bylake), transpose(Parasite_prevalence_bylake))
  colnames(Parasite_prevalence_bylake)<-as.character(unlist(Parasite_prevalence_bylake[1,]))
  Parasite_prevalence_bylake<-Parasite_prevalence_bylake[-1,]
  Parasite_prevalence_bylake[,lapply(.SD, function(x) sum(x>0.1 & x <0.9)),.SDcols=2:27]
  setnames(Parasite_prevalence_bylake,"site_name","Parasite_name")
  
  # choose MHC in each lake that has frequency ranging from 0.05-0.95
  MHC_prevalence_bylake<-MHC_protein_cleaned[,lapply(.SD, function(x) sum(x>0,na.rm=T)/(.N)), by='site_name',.SDcols=MHC_name_all]
  MHC_prevalence_bylake<-data.table(MHC_name = names(MHC_prevalence_bylake), transpose(MHC_prevalence_bylake))
  colnames(MHC_prevalence_bylake)<-as.character(unlist(MHC_prevalence_bylake[1,]))
  MHC_prevalence_bylake<-MHC_prevalence_bylake[-1,]
  MHC_prevalence_bylake[,lapply(.SD, function(x) sum(x>0.1 & x <0.9)),.SDcols=2:27]
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
MHC_name_lake<-list()
Parasite_name_lake<-list()
lake_fits<-list()
fit_res<-list()
Parasite_MHC_lake<-list()

for(lake in Lake_name){
  MHC_name_lake[[lake]]<-unlist(MHC_prevalence_bylake[get(lake)>0.1 & get(lake)<0.9,MHC_name])
  Parasite_name_lake[[lake]]<-unlist(Parasite_prevalence_bylake[get(lake)>0.1 & get(lake)<0.9,Parasite_name])
  MHC_temp<-MHC_protein_cleaned[site_name==lake,c(MHC_name_lake[[lake]],Parasite_name_lake[[lake]],"log_std_length"),with=F]
  lake_fits[[lake]]<-lapply(Parasite_name_lake[[lake]],function(y) glm(paste0(y,"~ log_std_length +",paste(MHC_name_lake[[lake]], collapse="+")),family="poisson",data=MHC_protein_cleaned[site_name==lake,]))
  fit_res[[lake]]<-sapply(lake_fits[[lake]], function(x) c(p=summary(x)$coef[,4],coef=summary(x)$coef[,1], z=summary(x)$coef[,3]))
  colnames(fit_res[[lake]])<-Parasite_name_lake[[lake]]
  Parasite_MHC_lake[[lake]]<-expand.grid(Parasite_name_lake[[lake]],MHC_name_lake[[lake]])
  
  
  #cor_temp<-rcorr(as.matrix(MHC_temp))
  #res_cor[[lake]]<-cor_temp
  #corrplot(cor_temp$r[MHC_name_lake[[lake]],Parasite_name_lake[[lake]]], type="full", order="original", 
           #p.mat = cor_temp$P[MHC_name_lake[[lake]],Parasite_name_lake[[lake]]], sig.level = 0.05, insig = "blank",
           #cl.lim=c(-1,1), col=colorRampPalette(c("blue","white","red"))(200))
}
p_res<-lapply(fit_res, function(x) t(x[grep("^p.prot",rownames(x)),]))
z_res<-lapply(fit_res, function(x) t(x[grep("^z.prot",rownames(x)),]))
coef_res<-lapply(fit_res, function(x) t(x[grep("^coef.prot",rownames(x)),]))

z_res_updated<-do.call(rbind,lapply(z_res,as.data.frame.table))
z_res_updated[,"lake"]<-sub("\\..*","", rownames(z_res_updated))
z_res_updated[,"Var2"]<-sub(".*z.","", z_res_updated[,"Var2"])
z_res_updated<-setDT(z_res_updated)
setnames(z_res_updated,c("Var1","Var2","Freq"),c("parasite","MHC","beta"))

MHC_prevalence_bylake_df<-as.data.frame(MHC_prevalence_bylake[,2:27])
rownames(MHC_prevalence_bylake_df)<-MHC_prevalence_bylake[,MHC_name]
MHC_prevalence_bylake_df<-setDT(as.data.frame.table(as.matrix(MHC_prevalence_bylake_df)))
setnames(MHC_prevalence_bylake_df,c("Var1","Var2","Freq"),c("MHC","lake","allele_freq"))

z_res_updated<-MHC_prevalence_bylake_df[z_res_updated,on=c("MHC","lake")]
z_res_updated[,allele_freq:=as.numeric(as.character(allele_freq))]

summary(glm("beta~allele_freq+lake+lake*allele_freq", data=z_res_updated, family="gaussian"))

for(i in unique(z_res_updated[,lake])){
  plot_data<-z_res_updated[lake==i,]
  plot(plot_data$allele_freq, plot_data$beta, data=plot_data, main=i)
  plot_model<-lm(plot_data$beta~plot_data$allele_freq)
  abline(plot_model)
  print(c(i,summary(plot_model)$coef[2,4]))
}

#### calculate parasite-MHC combinations
Parasite_MHC_combined<-rbindlist(lapply(seq_along(Parasite_MHC_lake),function(x) {
  df_tem<-Parasite_MHC_lake[[x]]
  df_tem["lake"]<-names(Parasite_MHC_lake)[[x]]
  return(df_tem)}))
Parasite_MHC_combined[,"combination":=as.factor(paste(Var1,Var2,sep="_"))]
freq_combination<-table(Parasite_MHC_combined[,combination])
freq_combination<-sort(freq_combination[freq_combination>1],decreasing=T)
freq_combination_df<-as.data.frame(t(sapply(names(freq_combination), function(x) unlist(strsplit(x,split="_prot_")))))
freq_combination_df["No.lakes"]<-freq_combination
freq_combination_df<-setDT(freq_combination_df,keep.rownames = T)
freq_combination_by_parasite<-freq_combination_df[,c("combinations"=list(list(rn))),by=V1]
setnames(freq_combination_by_parasite,"V1","parasite")
  
comb_fits<-list()
comb_fits_res<-list()
anova_res<-list()
for(comb in freq_combination_df[,rn]){
  Par_temp<-as.character(Parasite_MHC_combined[combination==comb,Var1[1]])
  MHC_allele_temp<-as.character(Parasite_MHC_combined[combination==comb,Var2[1]])
  lakes_temp<-Parasite_MHC_combined[combination==comb,lake]
  MHC_temp<-MHC_protein_cleaned[site_name %in% lakes_temp,c(Par_temp,MHC_allele_temp,"log_std_length","site_name"),with=F]
  comb_fits[[comb]]<-glm(paste0(Par_temp,"~ log_std_length + site_name +",MHC_allele_temp,"+ site_name*",MHC_allele_temp),
                         family="poisson",data=MHC_temp)
  anova_temp<-anova(comb_fits[[comb]], test = "Chisq")
  anova_res[[comb]]<- c(anova_temp[4:5,"Deviance"]/anova_temp[1,"Resid. Dev"],anova_temp[4:5,"Pr(>Chi)"])
  comb_fits_res[[comb]]<-as.matrix(c(p=summary(comb_fits[[comb]])$coef[,4],coef=summary(comb_fits[[comb]])$coef[,1], z=summary(comb_fits[[comb]])$coef[,3]))
}

div_per_res<-as.data.frame(do.call(rbind, anova_res))
colnames(div_per_res)<-c("D_a","D_as","p_a","p_as")
div_per_res$col<-"black"
div_per_res$log_D_a<-div_per_res$D_a
div_per_res$log_D_as<-div_per_res$D_as
div_per_res[which(div_per_res$p_a<0.05 & div_per_res$p_as<0.05),"col"]<-"purple"
div_per_res[which(div_per_res$p_a<0.05 & div_per_res$p_as>=0.05),"col"]<-"red"
div_per_res[which(div_per_res$p_a>=0.05 & div_per_res$p_as<0.05),"col"]<-"blue"

plot(div_per_res[,"log_D_as"]~div_per_res[,"log_D_a"],ylim=c(0,0.5),xlim=c(0,0.5), col=div_per_res$col, ylab = "allele*site deviance %", xlab="allele deviance %")
#abline(lm(div_per_res[,"log_D_as"]~div_per_res[,"log_D_a"]))
abline(0,1)

# select the significant interaction effect of each parasite
lapply(freq_combination_by_parasite[1:2,combinations], function(x){
  lapply(comb_fits_res[unlist(x)], function(y) {
    temp<-as.data.frame(y)
    temp[grepl("^p.prot",rownames(temp)) & temp["V1"]<0.05,,drop=F]})
})














Parasite_reg<-names(sort(table(unlist(Parasite_name_lake)),decreasing=TRUE)[1:10])
sort(table(unlist(MHC_name_lake)),decreasing=TRUE)[1:10]

for(par in Parasite_reg){
  
}


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
  
