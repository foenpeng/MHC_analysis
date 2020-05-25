library(data.table)
library(glmnet)
library(Hmisc)
library(corrplot)
library(reshape2)
library(seriation)
library(lme4)
library(lmerTest)
library(ape)

setwd("/Users/pengfoen/OneDrive - University of Connecticut/MHC")
MHC_protein <- fread("./7. Metadata and MHC data combined/7.0.BayesProteinsMergedMetadata.csv")
#MHC_gene <- fread("./7. Metadata and MHC data combined/7.0.BayesGenotypesMergedMetadata.csv")

MHC_name_all<-colnames(MHC_protein)[7:1210]
Parasite_name_all<-colnames(MHC_protein)[1271:1311]

### clean parasite data
# remove rows with no parasite info and lakes with too few samples
MHC_protein_cleaned<-na.omit(MHC_protein, cols=c("Parasite.richness"))
anyNA(MHC_protein_cleaned[,..Parasite_name_all])

# create site_ID from sample ID, becasue site name has many NAs (It turned out these NAs are the ones without metadata)
MHC_protein_cleaned[,site_ID:=unlist(lapply(SAMPLE_ID, function(x) strsplit(x, "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", perl=TRUE)[[1]][1]))]
MHC_protein_cleaned[,N_aas:=rowSums(.SD>0),.SDcols=MHC_name_all]

### explore the cut off threshold for low depth samples
# plot fig S1
# png("./Figures/figS1a.png",res=300, width=1000, height = 800)
# ggplot(MHC_protein_cleaned,aes(x=DEPTH_ALLELES)) +
#   geom_histogram(bins=80) +
#   labs(y="Sample Counts",x="Sequencing depth") +
#   geom_rect(aes(xmin=-Inf, xmax=450, ymin=-Inf, ymax=Inf), fill="grey",alpha=0.01) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
#   scale_x_continuous(breaks = c(0,450,1000,2000,3000,4000,5000))
# dev.off()

# png("./Figures/figS1b.png",res=300, width=1000, height = 800)
# ggplot(MHC_protein_cleaned,aes(x=DEPTH_ALLELES, y=N_aas)) +
#   geom_point(shape=1, size=1) +
#   geom_rect(aes(xmin=-Inf, xmax=450, ymin=-Inf, ymax=Inf), fill="grey",alpha=0.01) +
#   labs(x="Sequencing depth", y="Number of MHC alleles retrieved") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   scale_x_continuous(breaks = c(0,450,1000,2000,3000,4000,5000))
# dev.off()

# remove samples with low depth and lakes with too few data for regression
MHC_protein_cleaned<-MHC_protein_cleaned[,if(.N>1) .SD, by="site_ID"][DEPTH_ALLELES>450]
Lake_name<-unlist(unique(MHC_protein_cleaned[,"site_name"]))

# remove 89 MHC alleles that only present in the individuals removed from the previous step
MHC_remove<-MHC_name_all[t(MHC_protein_cleaned[,lapply(.SD, function(x) sum(x>0,na.rm=T)), .SDcols=MHC_name_all])==0]
MHC_protein_cleaned[,(MHC_remove):=NULL]
MHC_name<-colnames(MHC_protein_cleaned)[8:1122]

# convert MHC copy data into presence/absence
MHC_protein_cleaned[,(MHC_name):=lapply(.SD, function(x) as.numeric(x>0,na.rm=T)),.SDcols=MHC_name]

# { 
#   # old clean
# # count how many times each parasite appear in the table, there are many unfrequent ones
# Parasite_count_bylake<-MHC_protein[,lapply(.SD, function(x) sum(x>0,na.rm=T)),.SDcols=Parasite_name_all, by="site_name"][,rowSums(.SD, na.rm = TRUE),.SDcols=Parasite_name_all, by="site_name"]
# Parasite_appearance_bylake<-MHC_protein[,c(.N,lapply(.SD, function(x) sum(x>0,na.rm=T))),.SDcols=Parasite_name_all, by="site_name"][,lapply(.SD, function(x) sum(x>0,na.rm=T)),.SDcols=Parasite_name_all]
# Parasite_appearance_bylake<-data.table(parasite = names(Parasite_appearance_bylake), transpose(Parasite_appearance_bylake))
# #Parasite_name<-Parasite_appearance_bylake[V1>10,parasite] # keep parasite that appears in at least 11 lakes
# Parasite_name<-Parasite_name_all
# 
# Parasite_count<-MHC_protein[,lapply(.SD, function(x) sum(x>0,na.rm=T)),.SDcols=Parasite_name]
# Parasite_count<-data.table(cn = names(Parasite_count), transpose(Parasite_count))
# 
# 
# 
# ### clean MHC data
# MHC_appearance_bylake<-MHC_protein_cleaned[,c(.N,lapply(.SD, function(x) sum(x>0,na.rm=T))),.SDcols=MHC_name_all, by="site_name"][,lapply(.SD, function(x) sum(x>0,na.rm=T)),.SDcols=MHC_name]
# MHC_appearance_bylake<-data.table(MHC = names(MHC_appearance_bylake), transpose(MHC_appearance_bylake))
# #MHC_name<-MHC_appearance_bylake[V1>1,MHC] # require MHC appears in at least two lakes
# MHC_name<-MHC_name_all
# 
# # remove MHC with 0 appearance
# MHC_count<-MHC_protein_cleaned[,lapply(.SD, function(x) sum(x>0,na.rm=T)),.SDcols=MHC_name]
# MHC_count<-data.table(cn = names(MHC_count), transpose(MHC_count))
# 
# # look for highly correlated MHC and delete one of them
# temp_MHC<-MHC_protein_cleaned[,..MHC_name]
# temp_cor<-cor(temp_MHC)
# temp_cor[!lower.tri(temp_cor)] <- 0
# length(colnames(temp_MHC)[apply(temp_cor,2,function(x) any(abs(x) > 0.99))])
# MHC_name<-colnames(temp_MHC)[!apply(temp_cor,2,function(x) any(abs(x) > 0.99))]
#}

{
  # new clean 8/8/2019
  
  # choose parasite in each lake that has frequency ranging from 0.05-0.95
  Parasite_prevalence_bylake<-MHC_protein_cleaned[,lapply(.SD, function(x) sum(x>0,na.rm=T)/(.N)), by='site_name',.SDcols=Parasite_name_all]
  Parasite_prevalence_bylake<-data.table(Parasite_name = names(Parasite_prevalence_bylake), transpose(Parasite_prevalence_bylake))
  colnames(Parasite_prevalence_bylake)<-as.character(unlist(Parasite_prevalence_bylake[1,]))
  Parasite_prevalence_bylake<-Parasite_prevalence_bylake[-1,]
  Parasite_prevalence_bylake[,lapply(.SD, function(x) sum(x>0.05 & x <0.95)),.SDcols=2:27]
  setnames(Parasite_prevalence_bylake,"site_name","Parasite_name")
  
  # choose MHC in each lake that has frequency ranging from 0.05-0.95
  MHC_prevalence_bylake<-MHC_protein_cleaned[,lapply(.SD, function(x) sum(x>0,na.rm=T)/(.N)), by='site_name',.SDcols=MHC_name]
  MHC_prevalence_bylake<-data.table(MHC_name = names(MHC_prevalence_bylake), transpose(MHC_prevalence_bylake))
  colnames(MHC_prevalence_bylake)<-as.character(unlist(MHC_prevalence_bylake[1,]))
  MHC_prevalence_bylake<-MHC_prevalence_bylake[-1,]
  MHC_prevalence_bylake[,lapply(.SD, function(x) sum(x>0.05 & x <0.95)),.SDcols=2:27]
  setnames(MHC_prevalence_bylake,"site_name","MHC_name")
}

##### Basic statistics of MHC diversity
summary(MHC_protein_cleaned$DEPTH_ALLELES)
summary(MHC_protein_cleaned$N_aas)

N_samples<-t(MHC_protein_cleaned[,.N,by=site_ID][,N])
mean(N_samples)
sd(N_samples)

# calculate how many MHC alleles are present in only one site
MHC_prevalance_byallele<-t(MHC_protein_cleaned[,lapply(.SD, function(x) sum(x>0,na.rm=T)), by=site_ID, .SDcols=MHC_name][,lapply(.SD, function(x) sum(x>0,na.rm=T)), .SDcols=MHC_name])
sum(MHC_prevalance_byallele==1)/1115

changeCols <- c("site_ID","watershed")
MHC_protein_cleaned[,(changeCols):= lapply(.SD, as.factor), .SDcols = changeCols]
# test if data are normally distributed
shapiro.test(MHC_protein_cleaned$N_aas)
# aov does not work because shapiro test is significant
kruskal.test(N_aas ~ watershed,data=MHC_protein_cleaned)
kruskal.test(N_aas ~ site_ID,data=MHC_protein_cleaned)

# for parasite summary statistics
summary(MHC_protein_cleaned$Parasite.richness)
sd(MHC_protein_cleaned$Parasite.richness)
shapiro.test(MHC_protein_cleaned$Parasite.richness)
kruskal.test(Parasite.richness ~ watershed,data=MHC_protein_cleaned)
kruskal.test(Parasite.richness ~ site_ID,data=MHC_protein_cleaned)

############ Mantel test of matrix distance between parasite and MHC ################

P_dist<-dist(t(as.matrix(Parasite_prevalence_bylake[,-1],label=TRUE)))
M_dist<-dist(t(as.matrix(MHC_prevalence_bylake[,-1],label=TRUE)))
library(ecodist)
mantel(P_dist~M_dist)
dist_table<-data.table(Parasite=as.vector(P_dist), MHC=as.vector(M_dist))

## plot matrix
set.seed(2)
o <- seriate(as.matrix(P_dist), method="BEA_TSP")

P_longData<-melt(as.matrix(P_dist))
P_longData<-P_longData[P_longData$value!=0,]

P_longData$Var1 <- factor(P_longData$Var1, levels=names(unlist(o[[1]][]))) 
P_longData$Var2 <- factor(P_longData$Var2, levels=names(unlist(o[[2]][])))

M_longData<-melt(as.matrix(M_dist))
M_longData<-M_longData[M_longData$value!=0,]

M_longData$Var1 <- factor(M_longData$Var1, levels=names(unlist(o[[1]][]))) 
M_longData$Var2 <- factor(M_longData$Var2, levels=names(unlist(o[[2]][])))

png("./Figures/fig3a_P.png",res=300, height = 800, width = 1000)
ggplot(P_longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="grey20") +
  labs(x="",y="",title="Distance matrix of parasite") +
  theme_bw() + theme(axis.text.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks = element_blank(),
                     legend.title=element_blank(),
                     plot.title=element_text(size=11,hjust = 0.5))
dev.off()

png("./Figures/fig3a_M.png",res=300, height = 800, width = 1000)
ggplot(M_longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="grey20") +
  labs(x="",y="",title="Distance matrix of MHC alleles") +
  theme_bw() + theme(axis.text.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks = element_blank(),
                     legend.title=element_blank(),
                     plot.title=element_text(size=11,hjust = 0.5))
dev.off()


png("./Figures/fig3b.png",res=300, width=1000,height = 800)
ggplot(dist_table, aes(x=MHC, y =Parasite)) + 
  geom_point(size=0.3, color="grey30") +
  stat_ellipse() +  
  labs(x = "Value in MHC distance matrix", y = "Value in parasite distance matrix") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#scale_x_continuous(limits = c(0,3)) + 
#scale_y_continuous(limits = c(0,3))
dev.off()

###### use Principal coordinates of neighbour matrices (PCNM) to analyze multiple matrices

### prepare genomic distance matrix
Genomic_dist<-as.matrix(fread("PairwiseFst_LS45All.csv"))
rownames(Genomic_dist)<-Genomic_dist[,1]
genomic_dist<-Genomic_dist[,-1]
class(genomic_dist)<-"numeric"
# correct some names in the original file
setdiff(Lake_name,rownames(genomic_dist))
names_to_correct<-rownames(genomic_dist)
names_to_correct<-replace(names_to_correct,names_to_correct=="Higgins Lake", "Higgens Lake")
names_to_correct<-replace(names_to_correct,names_to_correct=="Pye Stream", "Pye Creek")
names_to_correct<-replace(names_to_correct,names_to_correct=="Pye Estuary", "Pye Outlet")
intersect(Lake_name,names_to_correct)
colnames(genomic_dist)<-names_to_correct
rownames(genomic_dist)<-names_to_correct
flip_value_to_triangular<-function(M){
    M2<-M
    for(i in 1:dim(M)[1]){
      for(j in 1:dim(M)[2]){
        if(i<j & is.na(M[j,i]) ){
          M[j,i]<-M2[i,j]
          M[i,j]<-NA
        }
      }
    }
  return(M)
}
genomic_dist_corrected<-as.dist(flip_value_to_triangular(genomic_dist[Lake_name,Lake_name]))

### prepare swimming distance matrix
Swim_dist<-as.matrix(fread("FishSwimDist.csv"))
rownames(Swim_dist)<-Swim_dist[,1]
swim_dist<-Swim_dist[,-1]
class(swim_dist)<-"numeric"
swim_dist_corrected<-as.dist(flip_value_to_triangular(swim_dist[Lake_name,Lake_name]))

# getting pcnm vectors from genomic, swimming and parasite distance matrices
library(vegan)
pcnm_parasite<-pcnm(P_dist)$vectors
pcnm_genome<-pcnm(genomic_dist_corrected)$vectors
pcnm_swim<-pcnm(swim_dist_corrected)$vectors
pcnm_mhc<-pcnm(M_dist)





############ Regression ############

# regression between lake averages of MHC and parasite
MHC_parasite_avg<-MHC_protein_cleaned[,.("avg_MHC"=mean(N_aas),"avg_parasite"=mean(Parasite.richness),"site_name"=site_name[1],"lake_area"=surface_area_ha[1],"benthicdiet"=mean(benthicdiet.score,na.rm=T)),by="site_ID"]
summary(lm(avg_MHC~avg_parasite,data=MHC_parasite_avg))

png("./Figures/fig1b.png",res=300, width=1000,height = 800)
ggplot(MHC_parasite_avg, aes(x=avg_parasite, y =avg_MHC)) + 
  geom_point(size=0.5) +
  stat_smooth(method="lm",formula= y~x, color = "black", size= 1) + 
  labs(x = "Mean value of parasite richness", y = "Mean number of MHC allelels") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()




# regression between MHC diversity and genomic heterozygosity
het<-fread("Heterozygosity_by_site_LS45.csv")
het[Pop=="Higgins Lake", Pop:="Higgens Lake"]
het[Pop=="Pye Stream", Pop:="Pye Creek"]
het[Pop=="Pye Estuary", Pop:="Pye Outlet"]
MHC_het<-het[MHC_parasite_avg,on="Pop==site_name"]
summary(lm(MHC_het$avg_MHC~MHC_het$MeanHet))
summary(lm(avg_MHC~MeanHet + avg_parasite + lake_area + benthicdiet, data=MHC_het))

png("./Figures/fig_MHC_het.png",res=300, width=1000,height = 800)
ggplot(MHC_het, aes(x=MeanHet, y =avg_MHC)) + 
  geom_point(size=0.5) +
  stat_smooth(method="lm",formula= y~x, color = "black", size= 1) + 
  labs(x = "Mean heterozygosity", y = "Mean number of MHC allelels") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()


# phylogeny
site_Fst<-fread("PairwiseFst_LS45All.csv")
NJ_tree<-njs(as.matrix(as.numeric(site_Fst)))


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
P_res<-list()
Z_res<-list()
Parasite_MHC_lake<-list()

for(lake in Lake_name){
  MHC_name_lake[[lake]]<-unlist(MHC_prevalence_bylake[get(lake)>0.05 & get(lake)<0.95,MHC_name])
  Parasite_name_lake[[lake]]<-unlist(Parasite_prevalence_bylake[get(lake)>0.05 & get(lake)<0.95,Parasite_name])
  MHC_temp<-MHC_protein_cleaned[site_name==lake,c(MHC_name_lake[[lake]],Parasite_name_lake[[lake]],"log_std_length"),with=F]
  #MHC_temp[,(Parasite_name_lake[[lake]]):=lapply(.SD, function(x) as.numeric(x>0,na.rm=T)),.SDcols=Parasite_name_lake[[lake]]]
  
  temp_mat<-matrix(nrow = length(Parasite_name_lake[[lake]]),
                   ncol = length(MHC_name_lake[[lake]]),
                   dimnames = list(Parasite_name_lake[[lake]],MHC_name_lake[[lake]]))
  z_temp<-temp_mat
  p_temp<-temp_mat
  for(par in Parasite_name_lake[[lake]]){
    temp_fit<-list()
    tryCatch({
      temp_fit <-lapply(MHC_name_lake[[lake]],function(x) glm.nb(paste0(par,"~ log_std_length +",x),link=log,data=MHC_temp))
      z_temp[par,]<-unlist(lapply(temp_fit, function(x) z=summary(x)$coef[3,3]))
      p_temp[par,]<-unlist(lapply(temp_fit, function(x) p=summary(x)$coef[3,4]))
    }, error=function(e){cat("ERROR :",conditionMessage(e),lake,par, "\n")})
  }
  P_res[[lake]]<-p_temp
  Z_res[[lake]]<-z_temp
  Parasite_MHC_lake[[lake]]<-expand.grid(Parasite_name_lake[[lake]],MHC_name_lake[[lake]])
 
  #cor_temp<-rcorr(as.matrix(MHC_temp))
  #res_cor[[lake]]<-cor_temp
  #corrplot(cor_temp$r[MHC_name_lake[[lake]],Parasite_name_lake[[lake]]], type="full", order="original", 
           #p.mat = cor_temp$P[MHC_name_lake[[lake]],Parasite_name_lake[[lake]]], sig.level = 0.05, insig = "blank",
           #cl.lim=c(-1,1), col=colorRampPalette(c("blue","white","red"))(200))
}

p_res_updated<-do.call(rbind,lapply(P_res,as.data.frame.table))
p_res_updated[,"lake"]<-sub("\\..*","", rownames(p_res_updated))
p_res_updated<-setDT(p_res_updated)
setnames(p_res_updated,c("Var1","Var2","Freq"),c("parasite","MHC","p"))

z_res_updated<-do.call(rbind,lapply(Z_res,as.data.frame.table))
z_res_updated[,"lake"]<-sub("\\..*","", rownames(z_res_updated))
z_res_updated<-setDT(z_res_updated)
setnames(z_res_updated,c("Var1","Var2","Freq"),c("parasite","MHC","z"))

z_p_combined<-z_res_updated[p_res_updated,on=c("lake","parasite","MHC")]

# test some models
pois<-summary(glm(formula=Diplostomum_spp2~log_std_length+prot_825,family="poisson",data=MHC_protein_cleaned[site_name=="Lawson Lake"]))
nb<-summary(glm.nb(formula=Diplostomum_spp2~log_std_length+prot_825,link=log,data=MHC_protein_cleaned[site_name=="Lawson Lake"]))
ggplot(aes(x=prot_825,y=Diplostomum_spp2),data=MHC_protein_cleaned[site_name=="Lawson Lake"]) +
  geom_jitter(width=0.1, height=0.05) +
  labs(title = "Lawson Lake, z=-3.11")

MHC_prevalence_bylake_df<-as.data.frame(MHC_prevalence_bylake[,2:27])
rownames(MHC_prevalence_bylake_df)<-MHC_prevalence_bylake[,MHC_name]
MHC_prevalence_bylake_df<-setDT(as.data.frame.table(as.matrix(MHC_prevalence_bylake_df)))
setnames(MHC_prevalence_bylake_df,c("Var1","Var2","Freq"),c("MHC","lake","allele_freq"))

z_MHC<-MHC_prevalence_bylake_df[z_p_combined,on=c("MHC","lake")]
z_MHC[,allele_freq:=as.numeric(as.character(allele_freq))]
z_MHC[,abs_z:=abs(z)]
  
allele_freq_reg<-summary(lm("p~allele_freq", data=z_MHC))

png("./Figures/fig2a.png",res=300, width=1000,height = 800)
ggplot(z_MHC, aes(x=allele_freq, y =p, x2 = lake)) + 
  geom_point(size=0.3, color="grey30") +
  stat_smooth(method="lm",formula= y~x, color = "black", size= 1, aes(group=1)) + 
  labs(x = "Allele frequency", y = "p value", title=paste("p=",round(allele_freq_reg$coef[2,4],5))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous( breaks = seq(0.1,0.9,by=0.1))
dev.off()


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
div_per_res$p_a.adjust<-p.adjust(div_per_res$p_a,method = "fdr", n=nrow(div_per_res))
div_per_res$p_as.adjust<-p.adjust(div_per_res$p_as,method = "fdr", n=nrow(div_per_res))
div_per_res$col<-"grey75"
div_per_res$log_D_a<-log10(div_per_res$D_a*100)
div_per_res$log_D_as<-log10(div_per_res$D_as*100)
div_per_res[which(div_per_res$p_a.adjust<0.05 & div_per_res$p_as.adjust<0.05),"col"]<-"purple"
div_per_res[which(div_per_res$p_a.adjust<0.05 & div_per_res$p_as.adjust>=0.05),"col"]<-"red"
div_per_res[which(div_per_res$p_a.adjust>=0.05 & div_per_res$p_as.adjust<0.05),"col"]<-"blue"

plot(div_per_res[,"D_as"]~div_per_res[,"D_a"], col=div_per_res$col, ylab = "allele*site deviance %", xlab="allele deviance %")
#abline(lm(div_per_res[,"log_D_as"]~div_per_res[,"log_D_a"]))
abline(0,1)

png("./Figures/fig2b.png", res = 300, width = 1000, height =800)
ggplot(div_per_res, aes(x = D_a*100, y = D_as*100)) + 
  geom_point(size = 0.3, color = div_per_res$col) + 
  labs(x = "Percentage of Deviation\n explained by MHC", y = "Percentage of Deviation\n explained by MHC * site") +
  geom_abline(intercept = 0, slope = 1, color="black", linetype="dashed", size=0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(limits = c(0,15)) + 
  scale_y_continuous(limits = c(0,15))
dev.off()

png("./Figures/fig2b_log_0.01.png", res = 300, width = 1000, height =800)
ggplot(div_per_res, aes(x = log_D_a, y = log_D_as)) + 
  geom_point(size = 0.7, color = div_per_res$col) + 
  labs(x = "Percentage of Deviation\n explained by MHC (log scale)", y = "Pergentage of Deviation\n explained by MHC * site (log scale)") +
  geom_abline(intercept = 0, slope = 1, color="black", linetype="dashed", size=0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(limits = c(-2,2)) + 
  scale_y_continuous(limits = c(-2,2))
dev.off()

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
  
