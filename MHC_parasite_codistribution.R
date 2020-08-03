library(data.table)
library(glmnet)
library(Hmisc)
library(corrplot)
library(reshape2)
library(seriation)
library(lme4)
library(lmerTest)


setwd("/Users/pengfoen/OneDrive - University of Connecticut/MHC")
MHC_protein <- fread("./7. Metadata and MHC data combined/7.0.BayesProteinsMergedMetadata.csv")

### clean parasite name and remove parasite with no appearance 
{
  setnames(MHC_protein, old = c("Dip_spath","Unionidae_internal","Ergasilus_internal","Proteocephalus.spp2"), 
         new = c("Diplostomum_spathaceum", "Unionidae","Ergasilus","Proteocephalus_spp2"))

  # Diplostomum spp1 and spp2 are actually a single species and should be merged into Diplostomum spp
  MHC_protein[,Diplostomum_spp:=Diplostomum_spp1+Diplostomum_spp2]
  MHC_protein[,c("Diplostomum_spp1","Diplostomum_spp2"):=NULL]
  
  Parasite_name_all<-colnames(MHC_protein)[1271:1311]
  Parasite_remove<-Parasite_name_all[unlist(MHC_protein[,lapply(.SD, function(x) sum(x>0,na.rm=T)==0), .SDcols=c(1271:1311)])]
  MHC_protein[, (Parasite_remove):= NULL]
  Parasite_name<-colnames(MHC_protein)[c(1271:1302,1372)]
  Parasite_name_plot<- unlist(lapply(Parasite_name, gsub, pattern = "_", replacement = " ", fixed = TRUE))
}  
  
# remove rows with no parasite info and lakes with too few samples
MHC_protein_cleaned<-na.omit(MHC_protein, cols=c("Parasite.richness"))
anyNA(MHC_protein_cleaned[,..Parasite_name])
MHC_protein_cleaned<-MHC_protein_cleaned[,if(.N>1) .SD, by="site_name"]


### Clean MHC data ###

### explore the cut off threshold for low depth samples
#plot fig S1
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

{
  # calculate the total number of alleles for each sample
  MHC_name_all<-colnames(MHC_protein)[7:1210]
  MHC_protein_cleaned[,N_aas:=rowSums(.SD>0),.SDcols=MHC_name_all]
  
  # remove lakes with too few data for regression and samples with low depth 
  MHC_protein_cleaned<-MHC_protein_cleaned[DEPTH_ALLELES>450,]
  Lake_name<-unlist(unique(MHC_protein_cleaned[,"site_name"]))
  MHC_protein_cleaned[,site_ID:=unlist(lapply(SAMPLE_ID, function(x) strsplit(x, "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", perl=TRUE)[[1]][1]))]
  
  # remove 89 MHC alleles that only present in the individuals removed from the previous step
  MHC_remove<-MHC_name_all[t(MHC_protein_cleaned[,lapply(.SD, function(x) sum(x>0,na.rm=T)), .SDcols=MHC_name_all])==0]
  MHC_protein_cleaned[,(MHC_remove):=NULL]
  MHC_name<-colnames(MHC_protein_cleaned)[8:1122]
  
  # convert MHC copy data into presence/absence
  MHC_protein_cleaned[,(MHC_name):=lapply(.SD, function(x) as.numeric(x>0,na.rm=T)),.SDcols=MHC_name]
}



### calculate paraste and MHC prevalance for each site ###
{
  # choose parasite in each lake that has frequency ranging from 0.05-0.95
  Parasite_prevalence_bylake<-MHC_protein_cleaned[,lapply(.SD, function(x) sum(x>0,na.rm=T)/(.N)), by='site_name',.SDcols=Parasite_name]
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
summary(MHC_protein_cleaned$Parasite.richness)

N_samples<-MHC_protein_cleaned[,.N,by=site_ID][,N]
sd(N_samples)
summary(N_samples)
#tableS1<-MHC_protein_cleaned[,.("site_ID"=site_ID[1],"habitat_type"=habitat_type[1],"watershed"=watershed[1],"sample_size"=.N),by=site_name]
#write.csv(tableS1,"./Tables/Table S1 Sample sites and sample size information.csv")

# calculate how many MHC alleles are present in only one site
MHC_prevalance_byallele<-t(MHC_protein_cleaned[,lapply(.SD, function(x) sum(x>0,na.rm=T)), by=site_ID, .SDcols=MHC_name][,lapply(.SD, function(x) sum(x>0,na.rm=T)), .SDcols=MHC_name])
sum(MHC_prevalance_byallele==1)/1115

# calculate how many parasite are present in only one site
Parasite_prevalance_byallele<-t(MHC_protein_cleaned[,lapply(.SD, function(x) sum(x>0,na.rm=T)), by=site_ID, .SDcols=Parasite_name][,lapply(.SD, function(x) sum(x>0,na.rm=T)), .SDcols=Parasite_name])

# ANOVA on MHC and parasite diversity
changeCols <- c("site_ID","watershed")
MHC_protein_cleaned[,(changeCols):= lapply(.SD, as.factor), .SDcols = changeCols]
summary(aov(N_aas ~ habitat_type + watershed + site_ID %in% watershed,data=MHC_protein_cleaned))
summary(aov(Parasite.richness ~ habitat_type + watershed + site_ID %in% watershed,data=MHC_protein_cleaned))

############ Regression Analysis of the number of MHC alleles ################
# quadratic regression for all data
full_model<-glmer( Parasite.richness ~ log_std_length + N_aas+I(N_aas^2) +  (1|site_name) + (N_aas|site_name) + (I(N_aas^2)|site_name), family = "poisson", data = MHC_protein_cleaned)
model1<-glmer( Parasite.richness ~ log_std_length + N_aas+I(N_aas^2) +  (1|site_name) + (N_aas|site_name), family = "poisson", data = MHC_protein_cleaned)
model2<-glmer( Parasite.richness ~ log_std_length + N_aas +  (1|site_name) + (N_aas|site_name) + (I(N_aas^2)|site_name) , family = "poisson", data = MHC_protein_cleaned)
model3<-glmer( Parasite.richness ~ log_std_length + I(N_aas^2) +  (1|site_name) + (I(N_aas^2)|site_name), family = "poisson", data = MHC_protein_cleaned)
model4<-glmer( Parasite.richness ~ log_std_length + I(N_aas^2) +  (1|site_name) , family = "poisson", data = MHC_protein_cleaned)
model5<-glmer( Parasite.richness ~ log_std_length + N_aas +  (1|site_name) + (N_aas|site_name) , family = "poisson", data = MHC_protein_cleaned)
model6<-glmer( Parasite.richness ~ log_std_length + N_aas + (1|site_name) , family = "poisson", data = MHC_protein_cleaned)
base_model<-glmer( Parasite.richness ~ log_std_length + (1|site_name) , family = "poisson", data = MHC_protein_cleaned)

no_int_model<-glm( Parasite.richness ~ log_std_length + x , family = "poisson", data = MHC_protein_cleaned)

residualrichness <- resid( glmer( Parasite.richness ~ log_std_length + (1|site_name) , "poisson", data= MHC_protein_cleaned,na.action=na.exclude))
MHC_protein_cleaned[,residual_rich:=residualrichness]

summary(lm(residual_rich~N_aas+I(N_aas^2),data=MHC_protein_cleaned[habitat_type=="Lake",]))

# plot full data between number of MHC alleles and parasite load
png("./Figures/fig_hd_qua.png",res = 300, width = 1000, height = 800)
ggplot(MHC_protein_cleaned, aes(group=N_aas,x=N_aas, y=residual_rich)) + 
  geom_jitter(shape=16, size = 0.2, position=position_jitter(0.2)) +
  stat_smooth(method = "lm", formula = y~x + I(x^2) , size = 1, aes(group=1), color='black') +
  stat_smooth(method = "lm", formula = y~x + I(x^2) , size = 0.1, aes(group=site_name), color='grey20',fill=NA) +
  labs(title="",  y = "Residual parasite richness") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous("Number of MHC alleles per individual", breaks = unique(MHC_protein_cleaned$N_aas))
dev.off()

# calculate the p value of the quadratic term for each individual fit line
p.vals = sapply(Lake_name, function(i) {
  coef(summary(lm(residual_rich~N_aas + I(N_aas^2), data=MHC_protein_cleaned[site_name==i, ])))[3,4]})

# regression between lake averages of MHC and parasite
MHC_parasite_avg<-MHC_protein_cleaned[,.("avg_MHC"=mean(N_aas),"avg_parasite"=mean(Parasite.richness),"site_name"=site_name[1],"lake_area"=surface_area_ha[1],"benthicdiet"=mean(benthicdiet.score,na.rm=T), "habitat_type"=habitat_type[1]),by="site_ID"]
popMHC_popParasite<-summary(lm(avg_MHC~avg_parasite  ,data=MHC_parasite_avg))

png("./Figures/fig_hd_pop.png",res=300, width=1000,height = 800)
ggplot(MHC_parasite_avg, aes(x=avg_parasite, y =avg_MHC)) + 
  geom_point(size=0.5) +
  stat_smooth(method="lm",formula= y~x, color = "black", size= 1) + 
  labs(x = "Average parasite richness per population", y = "Average number of MHC allelels \n per population") +
  #annotate(geom="text", x=4, y=9.5, label=paste("p =",round(popMHC_popParasite$coef[2,4],2))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size=10))
dev.off()


# regression between MHC diversity and genomic heterozygosity
het<-fread("Heterozygosity_by_site_LS45.csv")
het[Pop=="Higgins Lake", Pop:="Higgens Lake"]
het[Pop=="Pye Stream", Pop:="Pye Creek"]
het[Pop=="Pye Estuary", Pop:="Pye Outlet"]
MHC_het<-het[MHC_parasite_avg,on="Pop==site_name"]
popMHC_popHet<-summary(lm(MHC_het$avg_MHC~MHC_het$MeanHet))
summary(lm(avg_MHC~MeanHet , data=MHC_het))
summary(lm(avg_MHC~MeanHet + avg_parasite + lake_area + benthicdiet, data=MHC_het))

png("./Figures/fig_MHC_het.png",res=300, width=1000,height = 800)
ggplot(MHC_het, aes(x=MeanHet, y =avg_MHC)) + 
  geom_point(size=0.5) +
  stat_smooth(method="lm",formula= y~x, color = "black", size= 1) + 
  labs(x = "Population mean heterozygosity", y = "Average number of MHC allelels \n per population") +
  #annotate(geom="text", x=0.08, y=9.5, label=paste("p =",round(popMHC_popHet$coef[2,4],2))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size=10))
dev.off()


############################################################
# lm with single MHC as dependent variable
MHC_name_lake<-list()
Parasite_name_lake<-list()
lake_fits<-list()
P_res<-list()
Z_res<-list()
Parasite_MHC_lake<-list()
library(MASS)
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
    temp_fit2<-list()
    tryCatch({
      temp_fit <-lapply(MHC_name_lake[[lake]],function(x) glm.nb(paste0(par,"~ log_std_length +",x),link=log,data=MHC_temp))
      #temp_fit <-lapply(MHC_name_lake[[lake]],function(x) glm(paste0(par,"~ log_std_length +",x),family="binomial",data=MHC_temp))
      z_temp[par,]<-unlist(lapply(temp_fit, function(x) z=summary(x)$coef[3,3]))
      p_temp[par,]<-unlist(lapply(temp_fit, function(x) p=summary(x)$coef[3,4]))
      # test if the two regression are similar enough, if not, make it NA
      temp_fit2<-lapply(MHC_name_lake[[lake]],function(x) glm.nb(paste0(par,"~ log_std_length +",x),link=log,data=MHC_temp[!which.max(get(par))])) # only remove ONE data point
      z_temp2<-unlist(lapply(temp_fit2, function(x) tryCatch(summary(x)$coef[3,3],error=function(e) return(NA))))
      z_temp[par,]<-ifelse(abs(z_temp2-z_temp[par,])<0.5,z_temp[par,],NA)
      p_temp[par,]<-ifelse(is.na(z_temp[par,]),NA,p_temp[par,])
    }, error=function(e){cat("ERROR :",conditionMessage(e),lake,par, "\n")})
  }
  P_res[[lake]]<-p_temp
  Z_res[[lake]]<-z_temp
  Parasite_MHC_lake[[lake]]<-expand.grid(Parasite_name_lake[[lake]],MHC_name_lake[[lake]])
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
z_p_combined[,p_adjust:=p.adjust(p,method = "BH")]

################## plot a heat map for a particular lake
z_p_combined[,parasite:=Parasite_name_plot[match(z_p_combined[,parasite],Parasite_name)]]
png("./Figures/fig_heatmap.png",res=300, width=1200,height = 2000)
ggplot(z_p_combined[lake=="Lawson Lake"], aes(x=parasite, y=MHC, fill=z)) +
  geom_tile() + theme_bw() + coord_equal() +
  scale_fill_distiller(palette="RdBu", direction=1,na.value = "grey80") +
  labs(x="Parasite species", y="MHC alleles",title = "Lawson Lake") +
  guides(fill=guide_colorbar("Z value")) +
  theme(axis.text.x = element_text(angle = 90), axis.ticks = element_blank())
dev.off()

# g_color<-ggplot_build(g)
# g_data<-setDT(g_color$data[[1]])
# g_data[x==12&y==10,]
# g_data[x==12&y==17,]

# plot examples for the heatmap
png("./Figures/fig_eg1_heatmap.png",res=300, width=500,height = 800)
ggplot(aes(as.factor(x=prot_577),y=Unionidae),data=MHC_protein_cleaned[site_name=="Lawson Lake"]) +
  geom_jitter(width=0.1, height=0.05, size=0.5, col="#B2182B") +  ## for eg2: #3372B3 eg1:#B2182B
  geom_smooth(method="lm",formula= y~x, color = "#B2182B", size= 1,aes(group=1), se=F) + 
  labs(title = "", x = "prot_577", y = "No. of Unionidae parasite\n per fish") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_discrete( breaks = c("0","1"),labels=c("absent","present")) + 
  scale_y_continuous(trans = 'log2')
dev.off()


MHC_prevalence_bylake_df<-as.data.frame(MHC_prevalence_bylake[,2:27])
rownames(MHC_prevalence_bylake_df)<-MHC_prevalence_bylake[,MHC_name]
MHC_prevalence_bylake_df<-setDT(as.data.frame.table(as.matrix(MHC_prevalence_bylake_df)))
setnames(MHC_prevalence_bylake_df,c("Var1","Var2","Freq"),c("MHC","lake","allele_freq"))

z_MHC<-MHC_prevalence_bylake_df[z_p_combined,on=c("MHC","lake")]
z_MHC[,allele_freq:=as.numeric(as.character(allele_freq))]
z_MHC[,abs_z:=abs(z)]
z_MHC[,log_P:=-log10(p)]
  
summary(lmer("z~allele_freq + (allele_freq|lake) + (allele_freq|parasite)", data=z_MHC))

png("./Figures/fig_fd_ra.png",res=300, width=1000,height = 800)
ggplot(data=subset(z_MHC, !is.na(z)), aes(x=allele_freq, y =z, x2 = lake))+
  geom_point(size=0.2, color="grey50") +
  stat_smooth(method="lm",formula= y~x, color = "black", size= 0.7, aes(group=1)) + 
  labs(x = "Allele prevalence", y = "Effect size (Z value)") +
  #annotate(geom="text", x=0.8, y=2.7, label=paste("p=",round(allele_freq_reg$coef[2,4],2)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous( breaks = seq(0.1,0.9,by=0.1))
dev.off()

#### calculate parasite-MHC combinations
Parasite_MHC_combined<-rbindlist(lapply(seq_along(Parasite_MHC_lake),function(x) {
  df_tem<-Parasite_MHC_lake[[x]]
  df_tem["lake"]<-names(Parasite_MHC_lake)[[x]]
  return(df_tem)}))
Parasite_MHC_combined[,"combination":=as.factor(paste(Var1,Var2,sep="-"))]
freq_combination<-table(Parasite_MHC_combined[,combination])
freq_combination<-sort(freq_combination[freq_combination>1],decreasing=T)
freq_combination_df<-as.data.frame(t(sapply(names(freq_combination), function(x) unlist(strsplit(x,split="-prot_")))))
freq_combination_df["No.lakes"]<-freq_combination
freq_combination_df<-setDT(freq_combination_df,keep.rownames = T)
freq_combination_by_parasite<-freq_combination_df[,c("combinations"=list(list(rn))),by=V1]
setnames(freq_combination_by_parasite,"V1","parasite")

comb_fits<-list()
comb_fits_res<-list()
anova_res<-list()
for(comb in freq_combination_df[,rn]){
  tryCatch({
  Par_temp<-as.character(Parasite_MHC_combined[combination==comb,Var1[1]])
  MHC_allele_temp<-as.character(Parasite_MHC_combined[combination==comb,Var2[1]])
  lakes_temp<-Parasite_MHC_combined[combination==comb,lake]
  MHC_temp<-MHC_protein_cleaned[site_name %in% lakes_temp,c(Par_temp,MHC_allele_temp,"log_std_length","site_name"),with=F]
  comb_fits[[comb]]<-glm.nb(paste0(Par_temp,"~ log_std_length + site_name +",MHC_allele_temp,"+ site_name*",MHC_allele_temp),
                         link=log,data=MHC_temp)
  anova_temp<-anova(comb_fits[[comb]], test = "Chisq")
  anova_res[[comb]]<- c(anova_temp[4:5,"Deviance"]/anova_temp[1,"Resid. Dev"],anova_temp[4:5,"Pr(>Chi)"])
  comb_fits_res[[comb]]<-as.matrix(c(p=summary(comb_fits[[comb]])$coef[,4],coef=summary(comb_fits[[comb]])$coef[,1], z=summary(comb_fits[[comb]])$coef[,3]))
  },error=function(e){cat("ERROR :",conditionMessage(e),comb, "\n")})
}

div_per_res<-as.data.frame(do.call(rbind, anova_res))
colnames(div_per_res)<-c("D_a","D_as","p_a","p_as")
div_per_res$p_a.adjust<-p.adjust(div_per_res$p_a,method = "fdr", n=nrow(div_per_res))
div_per_res$p_as.adjust<-p.adjust(div_per_res$p_as,method = "fdr", n=nrow(div_per_res))
div_per_res$col<-"grey90"
div_per_res$log_D_a<-log10(div_per_res$D_a*100)
div_per_res$log_D_as<-log10(div_per_res$D_as*100)
div_per_res[which(div_per_res$p_a<0.05 & div_per_res$p_as<0.05),"col"]<-"black"
div_per_res[which(div_per_res$p_a<0.05 & div_per_res$p_as>=0.05),"col"]<-"black"
div_per_res[which(div_per_res$p_a>=0.05 & div_per_res$p_as<0.05),"col"]<-"black"

plot(div_per_res[,"D_as"]~div_per_res[,"D_a"], col=div_per_res$col, ylab = "allele*site deviance %", xlab="allele deviance %")
#abline(lm(div_per_res[,"log_D_as"]~div_per_res[,"log_D_a"]))
abline(0,1)

png("./Figures/fig_fd_ia.png", res = 300, width = 1000, height =800)
ggplot(div_per_res, aes(x = D_a*100, y = D_as*100)) + 
  geom_point(size = 0.3, color = div_per_res$col) + 
  geom_point(data=div_per_res[ "Nematode_spp4-prot_110",], aes(x = D_a*100, y = D_as*100), colour="black", size=1.2) +
  geom_text(data=div_per_res[ "Nematode_spp4-prot_110",], aes(x = D_a*100, y = D_as*100, label="C"),hjust=-0.2,vjust=-0.3,size=3) +
  geom_point(data=div_per_res[ "Bunoderina-prot_1215",], aes(x = D_a*100, y = D_as*100), colour="black", size=1.2) +
  geom_text(data=div_per_res[ "Bunoderina-prot_1215",], aes(x = D_a*100, y = D_as*100, label="B"),hjust=-0.2,vjust=-0.2,size=3) +
  labs(x = "Percentage of deviation\n explained by MHC", y = "Percentage of deviation\n explained by MHC x site") +
  geom_abline(intercept = 0, slope = 1, color="black", linetype="dashed", size=0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(limits = c(0,15)) + 
  scale_y_continuous(limits = c(0,15))
dev.off()

comb_plot<-"Bunoderina-prot_1215"
lake_plot<-as.factor(Parasite_MHC_combined[combination==comb_plot,lake])
par_plot<-as.character(Parasite_MHC_combined[combination==comb_plot,Var1[1]])
MHC_plot<-as.character(Parasite_MHC_combined[combination==comb_plot,Var2[1]])

png("./Figures/fig_fd_ia_eg1.png",res=300, width=1050,height = 800)
ggplot(aes(x=as.factor(get(MHC_plot)),y=get(par_plot),color=site_name),data=MHC_protein_cleaned[site_name %in%lake_plot,]) +
  geom_jitter(width=0.1, height=0.05, size=0.5) +  
  geom_smooth(method="lm",formula= y~x, size= 1,aes(group=site_name),se=F) + 
  labs(title = "", x = MHC_plot, y = "No. of Bunoderina\n per fish", color="site") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="top") +
  scale_x_discrete( breaks = c("0","1"),labels=c("absent","present")) +
  scale_y_continuous(trans = "log2") +
  scale_color_manual(values=c("#3372B3", "#B2182B"))
dev.off()

############ Mantel test of matrix distance between parasite and MHC ################
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

P_dist<-dist(t(as.matrix(Parasite_prevalence_bylake[,-1],label=TRUE)))
M_dist<-dist(t(as.matrix(MHC_prevalence_bylake[,-1],label=TRUE)))

library(ecodist)
mantel(P_dist~M_dist)
mantel(genomic_dist_corrected~M_dist)
dist_table<-data.table(Parasite=as.vector(P_dist), MHC=as.vector(M_dist), genome=as.vector(genomic_dist_corrected))

png("./Figures/fig3b_parasite.png",res=300, width=1000,height = 800)
ggplot(dist_table, aes(x=Parasite, y =MHC)) + 
  geom_point(size=0.3, color="grey30") +
  stat_ellipse() +  
  labs(x = "Parasite distance between populations", y = "MHC distance\n between populations") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#scale_x_continuous(limits = c(0,3)) + 
#scale_y_continuous(limits = c(0,3))
dev.off()

png("./Figures/fig3b_genome.png",res=300, width=1000,height = 800)
ggplot(dist_table, aes(x=genome, y =MHC)) + 
  geom_point(size=0.3, color="grey30") +
  stat_ellipse() +  
  labs(x = "Fst", y = "MHC distance\n between populations") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
#scale_x_continuous(limits = c(0,3)) + 
#scale_y_continuous(limits = c(0,3))
dev.off()


####### plot map and phylogeny ###############
###plot phylogeny
library(ape)
tree<-nj(genomic_dist_corrected)
nj_tree<-root(tree,"Sayward Estuary")
colors2 <- rep("black", Nedge(nj_tree))
colors2[which(nj_tree$edge[,2] %in% 1:26)] <- topo.colors(26)
is_tip <- nj_tree$edge[,2] <= length(nj_tree$tip.label)
ordered_tips <- nj_tree$edge[is_tip, 2]

png("./Figures/fig_phylogeny.png",res=300, width=1500,height = 2400)
plot(nj_tree, edge.color=colors2,edge.width=1.5, font=0.8, label.offset = 0.02)
tiplabels(pch = 21, bg=topo.colors(26)[match(c(1:26),ordered_tips)], cex = 1,adj=0.5)
dev.off()

#####plot the map
site_pos<-unique(MHC_protein_cleaned[,.(site_name,site_ID,utm_N,utm_E)])

library(proj4)
proj4string <- "+proj=utm +zone=10 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs "
pj <- project(site_pos[,.(utm_E,utm_N)], proj4string, inverse=TRUE)
site_pos[,c("long","lat"):=pj]

png("./Figures/fig_mab_text.png",res=300, width=600,height = 1200)
layout(matrix(c(1,1,1,
                2,2,2,
                3,3,3),nrow=3,byrow=T),heights=c(6,0.5,5))
long <- c(-126,-125)
lat <- c(48.5,50.7)
center <- c(mean(lat), mean(long))
zoom <- 7
library(RgoogleMaps)
terrmap <- GetMap(center=center, zoom=zoom, maptype= "terrain", destfile = "terrain.png")
PlotOnStaticMap(terrmap, lat = c(49.9,49.9,50.5,50.5,49.9), lon = c(-126.2,-125,-125,-126.2,-126.2), lwd=1, FUN = lines)
plot.new()
long <- c(-126.2,-125)
lat <- c(49.9,50.5)
center <- c(mean(lat), mean(long))
zoom <- 10
terrmap2 <- GetMap(center=center, zoom=zoom, maptype= "terrain", destfile = "terrain2.png")
PlotOnStaticMap(terrmap2, lat = site_pos$lat, lon = site_pos$long, pch = 21, cex = 1, col="black",bg=topo.colors(26)[match(c(1:26),ordered_tips)])
PlotOnStaticMap(terrmap2, lat = site_pos$lat-0.015, lon = site_pos$long, FUN= text,labels = site_pos$site_ID,cex=0.5, add=TRUE) 
box(lwd = 1.5,col="black")
par(mfrow=c(1,1))
dev.off()
