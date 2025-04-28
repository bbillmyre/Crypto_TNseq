setwd("E:/rb2091/tnseq/Crypto_TN_seq_paper_final_analysis/predicting_essentials")
library(randomForest)
library("ggplot2")
library(readxl)
library("cowplot")

genes=read.table("all_unique_genes.txt")
GFF=read.table("FungiDB-44_CneoformansH99.gff",sep="\t")
GFF=subset(GFF,V3=="CDS")

MSEI= read.table("T0_MseI_out.labeled.bed")
CVIAII=read.table("T0_CviAII_out.labeled.bed")

combined_T0 = merge(MSEI,CVIAII,by=c("V1","V2","V3","V6","V8"),all=TRUE)

combined_T0$V4.x[ is.na(combined_T0$V4.x)] = 0
combined_T0$V5.x[ is.na(combined_T0$V5.x)] = 0
combined_T0$V4.y[ is.na(combined_T0$V4.y)] = 0
combined_T0$V5.y[ is.na(combined_T0$V5.y)] = 0

combined_T0=combined_T0[c("V1","V2","V6","V4.x","V5.x","V4.y","V5.y")]
combined_T0$total_inserts=combined_T0$V5.x+combined_T0$V5.y
combined_T0$ave_freq=(combined_T0$V4.x+combined_T0$V4.y)/2

combined_T0=combined_T0[c("V1","V2","V6","total_inserts","ave_freq")]
colnames(combined_T0)=c("chr","coord","strand","total_inserts_T0","ave_freq_T0")

MSEI_DMSO= read.table("../../flucon_test_redo/strict_matching/T2_DMSO_MseI_out.labeled.bed")
CVIAII_DMSO=read.table("../../flucon_test_redo/strict_matching/T2_DMSO_CviAII_out.labeled.bed")

combined_DMSO = merge(MSEI_DMSO,CVIAII_DMSO,by=c("V1","V2","V3","V6","V8"),all=TRUE)

combined_DMSO$V4.x[ is.na(combined_DMSO$V4.x)] = 0
combined_DMSO$V5.x[ is.na(combined_DMSO$V5.x)] = 0
combined_DMSO$V4.y[ is.na(combined_DMSO$V4.y)] = 0
combined_DMSO$V5.y[ is.na(combined_DMSO$V5.y)] = 0

combined_DMSO=combined_DMSO[c("V1","V2","V6","V4.x","V5.x","V4.y","V5.y")]
combined_DMSO$total_inserts=combined_DMSO$V5.x+combined_DMSO$V5.y
combined_DMSO$ave_freq=(combined_DMSO$V4.x+combined_DMSO$V4.y)/2

combined_DMSO=combined_DMSO[c("V1","V2","V6","total_inserts","ave_freq")]
colnames(combined_DMSO)=c("chr","coord","strand","total_inserts_DMSO","ave_freq_DMSO")


MSEI_IC50= read.table("../../flucon_test_redo/strict_matching/T2_Flu_IC50_MseI_out.labeled.bed")
CVIAII_IC50=read.table("../../flucon_test_redo/strict_matching/T2_Flu_IC50_CviAII_out.labeled.bed")

combined_IC50 = merge(MSEI_IC50,CVIAII_IC50,by=c("V1","V2","V3","V6","V8"),all=TRUE)

combined_IC50$V4.x[ is.na(combined_IC50$V4.x)] = 0
combined_IC50$V5.x[ is.na(combined_IC50$V5.x)] = 0
combined_IC50$V4.y[ is.na(combined_IC50$V4.y)] = 0
combined_IC50$V5.y[ is.na(combined_IC50$V5.y)] = 0

combined_IC50=combined_IC50[c("V1","V2","V6","V4.x","V5.x","V4.y","V5.y")]
combined_IC50$total_inserts=combined_IC50$V5.x+combined_IC50$V5.y
combined_IC50$ave_freq=(combined_IC50$V4.x+combined_IC50$V4.y)/2

combined_IC50=combined_IC50[c("V1","V2","V6","total_inserts","ave_freq")]
colnames(combined_IC50)=c("chr","coord","strand","total_inserts_IC50","ave_freq_IC50")


##have 3 combo files, need to start combining them, 


df_list <- list(combined_T0, combined_DMSO, combined_IC50)

all_combined_frame=Reduce(function(x, y) merge(x, y, all=TRUE), df_list)


all_combined_frame$total_inserts_T0[is.na(all_combined_frame$total_inserts_T0)] = 0
all_combined_frame$ave_freq_T0[is.na(all_combined_frame$ave_freq_T0)] = 0

all_combined_frame$total_inserts_DMSO[is.na(all_combined_frame$total_inserts_DMSO)] = 0
all_combined_frame$ave_freq_DMSO[is.na(all_combined_frame$ave_freq_DMSO)] = 0

all_combined_frame$total_inserts_IC50[is.na(all_combined_frame$total_inserts_IC50)] = 0
all_combined_frame$ave_freq_IC50[is.na(all_combined_frame$ave_freq_IC50)] = 0



##do some depth filtering and set 0 values to min

all_combined_frame=subset(all_combined_frame, total_inserts_T0>5)

##set 0 values to min
all_combined_frame$total_inserts_DMSO[all_combined_frame$total_inserts_DMSO == 0] = min(subset(all_combined_frame, total_inserts_DMSO!=0)$total_inserts_DMSO)
all_combined_frame$total_inserts_IC50[all_combined_frame$total_inserts_IC50 == 0] = min(subset(all_combined_frame, total_inserts_IC50!=0)$total_inserts_IC50)

all_combined_frame$ave_freq_DMSO[all_combined_frame$ave_freq_DMSO == 0] = min(subset(all_combined_frame, ave_freq_DMSO!=0)$ave_freq_DMSO)
all_combined_frame$ave_freq_IC50[all_combined_frame$ave_freq_IC50 == 0] = min(subset(all_combined_frame, ave_freq_IC50!=0)$ave_freq_IC50)



##label with genes

all_combined_frame$Gene="Intergenic"

for (gene in genes$V1){
  
  cur_gene=GFF[grep(gene,GFF$V9), ] 
  
  cur_min=min(cur_gene$V4)
  cur_max=max(cur_gene$V5)
  cur_chr=head(cur_gene$V1,1)
  cur_name=substr(head(cur_gene$V9,1),4,13)
  
  all_combined_frame$Gene[all_combined_frame$chr==cur_chr & all_combined_frame$coord<cur_max & all_combined_frame$coord>cur_min]=cur_name
  
}


stats_frame= data.frame(GeneID=NULL, Mean50_DMSO=NULL, P50_DMSO=NULL,
                        Mean50_DMSO_5=NULL, P50_DMSO_5=NULL, 
                        LenSite=NULL, Len_5=NULL)

intergenic= subset(all_combined_frame, Gene=="Intergenic")

int_gly_freq   = intergenic$ave_freq_T0
int_dmso_freq  = intergenic$ave_freq_DMSO
int_ic50_freq  = intergenic$ave_freq_IC50
int_ic50_dmso= log10(int_ic50_freq/int_dmso_freq)

for (gene in genes$V1){
  
  cur_gene=GFF[grep(gene,GFF$V9), ] 
  
  cur_min=min(cur_gene$V4)
  cur_max=max(cur_gene$V5)
  cur_chr=head(cur_gene$V1,1)
  
  
  if(head(cur_gene,1)$V7=="+"){
    cur_gene_5_min= min(cur_gene$V4)-301
    cur_gene_5_max= min(cur_gene$V4)-1 
    cur_gene_3_min= max(cur_gene$V5)+1
    cur_gene_3_max= max(cur_gene$V5)+301
  }else{
    cur_gene_5_min= max(cur_gene$V5)+1
    cur_gene_5_max= max(cur_gene$V5)+301  
    cur_gene_3_min= min(cur_gene$V4)-301
    cur_gene_3_max= min(cur_gene$V4)-1
    
  }
  
  cur_gene = subset(all_combined_frame, chr==cur_chr & coord < cur_max & coord > cur_min)
  cur_gene_5 = subset(all_combined_frame, chr==cur_chr & coord < cur_gene_5_max & coord > cur_gene_5_min)
  
  gly_freq   = cur_gene$ave_freq_T0
  dmso_freq  = cur_gene$ave_freq_DMSO
  ic50_freq  = cur_gene$ave_freq_IC50
  
  gly_freq_5   = cur_gene_5$ave_freq_T0
  dmso_freq_5  = cur_gene_5$ave_freq_DMSO
  ic50_freq_5  = cur_gene_5$ave_freq_IC50
  
  ic50_to_dmso= log10(ic50_freq/dmso_freq) 
  ic50_to_dmso_5= log10(ic50_freq_5/dmso_freq_5) 
  
  
  ave_50_to_dmso= mean(ic50_to_dmso)
  ave_50_to_dmso_5= mean(ic50_to_dmso_5)
  
  

  
 
  
  if(ave_50_to_dmso!="NaN" & ave_50_to_dmso!="Inf" & ave_50_to_dmso!="-Inf" & length(ic50_freq)>1){
    
    stat_result= wilcox.test(ic50_to_dmso,int_ic50_dmso,alternative="two.sided")
    pval_50_DMSO=stat_result$p.value
  } else { pval_50_DMSO=NA }
  

  
  if(ave_50_to_dmso_5!="NaN" & ave_50_to_dmso_5!="Inf" & ave_50_to_dmso_5!="-Inf" & length(ic50_freq_5)>1){
    
    stat_result= wilcox.test(ic50_to_dmso_5,int_ic50_dmso,alternative="two.sided")
    pval_50_DMSO_5=stat_result$p.value
  } else { pval_50_DMSO_5=NA }
  
  
  
  
  stats_frame=rbind(stats_frame,data.frame(GeneID=gene,  Mean50_DMSO=ave_50_to_dmso, P50_DMSO=pval_50_DMSO, 
                                            Mean50_DMSO_5=ave_50_to_dmso_5, P50_DMSO_5=pval_50_DMSO_5,
                                           LenSite=nrow(cur_gene),Len_5=length(ic50_to_dmso_5) ))
  print(tail(stats_frame,1))                                                             
}


gene_stats_frame=data.frame(Gene=NULL,Chromosome=NULL,Length=NULL,Num_Inserts=NULL,Num_Inserts_per_Length=NULL,
                            Frequency= NULL, Frequency_per_Length=NULL, Largest_gap= NULL, Largest_gap_perc=NULL,
                            Middle_Frequency=NULL,Middle_Frequency_per_Length=NULL,Freq_Local_100k=NULL)


for (gene in genes$V1){
  
  cur_gene=GFF[grep(gene,GFF$V9), ] 
  
  cur_min=min(cur_gene$V4)
  cur_max=max(cur_gene$V5)
  cur_chr=head(cur_gene$V1,1)
  
  cur_MSEI= subset(MSEI,V1==cur_chr & V2 < cur_max & V2 > cur_min)
  cur_CVIAII= subset(CVIAII,V1==cur_chr & V2 < cur_max & V2 > cur_min)
  
  
  ##combine frequencies from above files
  
  cur_merge= merge(cur_MSEI,cur_CVIAII,by=c("V1","V2","V3","V6","V8"),all=TRUE)
  
  cur_merge$V4.x[ is.na(cur_merge$V4.x)] = 0
  cur_merge$V5.x[ is.na(cur_merge$V5.x)] = 0
  cur_merge$V4.y[ is.na(cur_merge$V4.y)] = 0
  cur_merge$V5.y[ is.na(cur_merge$V5.y)] = 0
  
  cur_merge=cur_merge[c("V1","V2","V6","V4.x","V5.x","V4.y","V5.y")]
  cur_merge$total_inserts=cur_merge$V5.x+cur_merge$V5.y
  cur_merge$ave_freq=(cur_merge$V4.x+cur_merge$V4.y)/2
  
  ##local stats
  local_MSEI = subset(MSEI,V1==cur_chr & V2 < cur_max+50000 & V2 > cur_min-50000)
  local_CVIAII= subset(CVIAII,V1==cur_chr & V2 < cur_max+50000 & V2 > cur_min-50000)
  local_merge= merge(local_MSEI,local_CVIAII,by=c("V1","V2","V3","V6","V8"),all=TRUE)
  
  local_merge$V4.x[ is.na(local_merge$V4.x)] = 0
  local_merge$V5.x[ is.na(local_merge$V5.x)] = 0
  local_merge$V4.y[ is.na(local_merge$V4.y)] = 0
  local_merge$V5.y[ is.na(local_merge$V5.y)] = 0
  
  local_merge=local_merge[c("V1","V2","V6","V4.x","V5.x","V4.y","V5.y")]
  local_merge$total_inserts=local_merge$V5.x+local_merge$V5.y
  local_merge$ave_freq=(local_merge$V4.x+local_merge$V4.y)/2
  
  
  
  ##build stats
  
  length = cur_max-cur_min
  
  ##build largest gap
  gap=min(cur_merge$V2)-cur_min
  if (cur_max-max(cur_merge$V2)>gap){
    gap=cur_max-max(cur_merge$V2)
  }
  
  cur_merge_gap=cur_merge
  cur_merge_gap$gap <- cur_merge$V2 - c(1L, head(cur_merge$V2, -1))
  cur_merge_gap = cur_merge_gap[-1,]
  
  if (max(cur_merge_gap$gap)>gap){
    gap=max(cur_merge_gap$gap)
  }
  
  if (gap==Inf){
    gap=length
  }
  
  ##build middle 80% set
  mid_start= cur_min + (0.1*length)
  mid_end = cur_max - (0.1*length)
  
  mid_merge=subset(cur_merge,V1==cur_chr & V2 < mid_end & V2 > mid_start)
  
  cur_stats=data.frame(Gene=gene,Chromosome=cur_chr,Length=length,Num_Inserts=sum(cur_merge$total_inserts),Num_Inserts_per_Length=(sum(cur_merge$total_inserts)/length),
                       Frequency= sum(cur_merge$ave_freq), Frequency_per_Length= (sum(cur_merge$ave_freq)/length),Largest_gap= gap , Largest_gap_perc=gap/length,
                       Middle_Frequency=sum(mid_merge$ave_freq),Middle_Frequency_per_Length=(sum(mid_merge$ave_freq)/length) , Freq_Local_100k= sum(local_merge$ave_freq))
  
  gene_stats_frame=rbind(gene_stats_frame,cur_stats)
  
}


##stat building done, now build machine learning model





##for conserved genes

ortho_table=read.table(file = "../../machine_learning/GenesByOrthologs_candida.txt",header=TRUE,sep="\t")
other_essentials = read.table(file="../../machine_learning/candida_essentials.csv",header=TRUE,sep="\t")

colnames(ortho_table)=c("Crypto_name","tran","paralog","orthocount","name","empty")

crypto_ess_table=merge(ortho_table,other_essentials,by="name")
crypto_ess_table=subset(crypto_ess_table,select=-c(tran,empty))


conserved_ne=subset(crypto_ess_table, CaTn =="NE" & ScTn=="NE" & SpTn=="NE")
conserved_ne$essential=FALSE
conserved_ess=subset(crypto_ess_table, CaTn =="Ess" & ScTn=="Ess" & SpTn=="Ess")
conserved_ess$essential=TRUE
conserved=rbind(conserved_ne, conserved_ess)

training_all_data=subset(gene_stats_frame, Gene %in%conserved$Crypto_name)
training_all_data=training_all_data[order(training_all_data$Gene),]
conserved=conserved[!duplicated(conserved), ]
conserved=conserved[order(conserved$Crypto_name),]

training_all_data$conserved_essentials = conserved$essential  

##end conserved genes




##precision recall and full prediction/features
library(PRROC)

pre_rec_curve=data.frame(run_no =NULL, Threshold= NULL, Prec_Cons = NULL, Rec_Cons=NULL, AUROC=NULL, AUPRC=NULL)
essentials_scores=data.frame(Gene=NULL,essential=NULL)
features_total=data.frame(Parameter=NULL,MeanDecreaseAccuracy=NULL)
set.seed(100)

for(rep in 1:100){
  
  train <- sample(nrow(training_all_data),0.8*nrow(training_all_data),replace=FALSE)
  TrainSet <- training_all_data[train,]
  ValidSet <- training_all_data[-train,]
  
  
  
  
  
  TrainSet_Cons=TrainSet
  ValidSet_Cons=ValidSet
  TrainSet_Cons$conserved_essentials=as.factor(TrainSet$conserved_essentials)
  TrainSet_Cons = subset(TrainSet_Cons, select = -c(Gene) )
  rf_model_Cons <- randomForest(conserved_essentials ~ ., data =TrainSet_Cons,mtry=5, ntree =10000, importance =TRUE)
  rf_model_Cons
  predAllprob_Cons <- predict(rf_model_Cons, ValidSet_Cons, type="prob")
  ValidSet_Cons$pred=predAllprob_Cons[,2]
  
  
  
  for(cutoff in seq(0.01, 0.99, by=0.01) ){
    
    
    True_Pos_Cons=nrow(subset(ValidSet_Cons,pred>cutoff & conserved_essentials ==TRUE))
    False_Pos_Cons=nrow(subset(ValidSet_Cons,pred>cutoff & conserved_essentials ==FALSE))
    False_Neg_Cons=nrow(subset(ValidSet_Cons,pred<cutoff & conserved_essentials ==TRUE))
    
    Precision_Cons= (True_Pos_Cons)/(True_Pos_Cons+False_Pos_Cons)
    Recall_Cons= (True_Pos_Cons)/(True_Pos_Cons+False_Neg_Cons)
    
    Pos_cons=subset(ValidSet_Cons,conserved_essentials==TRUE)$pred
    Neg_cons=subset(ValidSet_Cons,conserved_essentials==FALSE)$pred
    auroc=roc.curve(scores.class0 = Pos_cons, scores.class1 = Neg_cons)$auc
    auprc=pr.curve(scores.class0 = Pos_cons, scores.class1 = Neg_cons)$auc.integral
    
    pre_rec_curve=rbind(pre_rec_curve, data.frame(run_no =rep, Threshold= cutoff, Prec_Cons = Precision_Cons, Rec_Cons=Recall_Cons, AUROC=auroc, AUPRC=auprc))
    
    
    
    
    
  }
  
  #run full predictions
  predAllprob <- predict(rf_model_Cons, gene_stats_frame, type="prob")
  gene_stats_frame$essential=predAllprob[,2]
  temp_essentials=subset(gene_stats_frame, select= c("Gene","essential"))
  essentials_scores=rbind(essentials_scores,temp_essentials)
  
  ##features
  features=data.frame(rf_model_Cons$importance)
  features=cbind(rownames(features),data.frame(features,row.names=NULL))
  colnames(features)=c("Parameter","FALSE","TRUE","MeanDecreaseAccuracy","MeanDecreaseGini")
  temp_features= subset(features, select= c("Parameter","MeanDecreaseAccuracy"))
  features_total=rbind(features_total,temp_features)
  
  print(rep)
}


mean_AUROC=mean(pre_rec_curve$AUROC)
mean_AUPRC=mean(pre_rec_curve$AUPRC)


##find means for PRC curve building
prec_recall_means=data.frame(Threshold=NULL, Mean_Precision_Cons=NULL, Mean_Recall_Cons=NULL)
for(thresh in  levels(as.factor(pre_rec_curve$Threshold))){
  
  mean_prec_cons= mean(subset(pre_rec_curve,Threshold==thresh)$Prec_Cons)
  mean_rec_cons =mean(subset(pre_rec_curve,Threshold==thresh)$Rec_Cons)
  
  prec_recall_means=rbind(prec_recall_means, data.frame(Threshold=thresh,Mean_Precision_Cons=mean_prec_cons, Mean_Recall_Cons=mean_rec_cons))
  
}



baseline=nrow(subset(training_all_data,conserved_essentials==TRUE))/(nrow(training_all_data))

ggplot(prec_recall_means)  + geom_point(aes(x=Mean_Recall_Cons,y=Mean_Precision_Cons))+geom_hline(aes(yintercept=baseline),linetype="dashed",color="blue")  + theme_bw()


##find means for gene predictions
essential_score_table=data.frame(Gene=NULL, Essential = NULL, SD = NULL , pval= NULL)
for(genes in levels(as.factor(essentials_scores$Gene))){
  
  print (genes)
  
  mean_essential=mean(subset(essentials_scores,Gene==genes)$essential)
  sd_essential=sd(subset(essentials_scores,Gene==genes)$essential)
  stat_result=t.test(subset(essentials_scores,Gene==genes)$essential,mu=0.5,alternative="two.sided")
  
  essential_score_table=rbind(essential_score_table,data.frame(Gene=genes,Essential=mean_essential,SD=sd_essential,pval=stat_result$p.value))
}

num_genes=nrow(essential_score_table)

essential_score_table$ess= ifelse(essential_score_table$Essential<0.5 & (essential_score_table$pval *num_genes) <0.05, "NESS", ifelse(essential_score_table$Essential>0.5 & (essential_score_table$pval * num_genes) <0.05, "ESS","UNK"))





##find means for features

parameters_score_table=data.frame(Parameter=NULL, Accuracy = NULL, SD = NULL )
for(parameters in levels(as.factor(features_total$Parameter))){
  
  print (parameters)
  
  mean_accuracy=mean(subset(features_total,Parameter==parameters)$MeanDecreaseAccuracy)
  sd_accuracy=sd(subset(features_total,Parameter==parameters)$MeanDecreaseAccuracy)
  
  parameters_score_table=rbind(parameters_score_table,data.frame(Parameter=parameters,Accuracy=mean_accuracy,SD=sd_accuracy))
}







##build master table

##find ortho information
orthos=read.csv(file= "../../Crypto_TN_seq_paper_final_analysis/predicting_essentials/GenesByText_Summary.txt",sep="\t")
colnames(orthos)=c("Gene","Location","Description","Orthologs","Paralogs","OG_Group")

orthos=subset(orthos, Gene %in%essential_score_table$Gene)


master_table=merge(orthos,essential_score_table,by="Gene")

names(stats_frame)[names(stats_frame) == "GeneID"] = "Gene"
master_table =merge(master_table, stats_frame, by ="Gene")

#write.csv(master_table,file="master_table.tsv", quote=TRUE,sep="\t")



##pull in madhani deletions
hiten_valid=read.delim("../../Crypto_TN_seq_paper_final_analysis/predicting_essentials/hiten_validated_kos.csv",sep=",")


master_table_deleted=subset(master_table, Gene %in%hiten_valid$Gene)
master_table_not_deleted=subset(master_table, !Gene %in%hiten_valid$Gene)

master_table_deleted$Deleted="TRUE"
master_table_not_deleted$Deleted="FALSE"
master_table=rbind(master_table_deleted,master_table_not_deleted)



##compare to other species
ortho_table=read.table(file = "../../machine_learning/GenesByOrthologs_candida.txt",header=TRUE,sep="\t")
other_essentials = read.table(file="../../machine_learning/candida_essentials.csv",header=TRUE,sep="\t")

colnames(ortho_table)=c("Gene","tran","paralog","orthocount","name","empty")

crypto_ess_table=merge(ortho_table,other_essentials,by="name")

crypto_ess_table=subset(crypto_ess_table,select=-c(tran,empty))




crypto_full_table = merge(crypto_ess_table,master_table,by="Gene")


##start round V2 model
##for conserved genes

ortho_table=read.table(file = "../../machine_learning/GenesByOrthologs_candida.txt",header=TRUE,sep="\t")
other_essentials = read.table(file="../../machine_learning/candida_essentials.csv",header=TRUE,sep="\t")

colnames(ortho_table)=c("Crypto_name","tran","paralog","orthocount","name","empty")

crypto_ess_table=merge(ortho_table,other_essentials,by="name")
crypto_ess_table=subset(crypto_ess_table,select=-c(tran,empty))


conserved_ne=subset(crypto_ess_table, CaTn =="NE" & ScTn=="NE" & SpTn=="NE")
conserved_ne$essential=FALSE
conserved_ess=subset(crypto_ess_table, CaTn =="Ess" & ScTn=="Ess" & SpTn=="Ess")
conserved_ess$essential=TRUE
conserved=rbind(conserved_ne, conserved_ess)

gene_stats_frame$essential=NULL
gene_stats_frame_2=gene_stats_frame

training_all_data=subset(gene_stats_frame, Gene %in%conserved$Crypto_name)
training_all_data=training_all_data[order(training_all_data$Gene),]
conserved=conserved[!duplicated(conserved), ]
conserved=conserved[order(conserved$Crypto_name),]

training_all_data$conserved_essentials = conserved$essential  

##end conserved genes
temp1=subset(master_table,Gene %in% unique(subset(crypto_full_table,(CaTn=="NE" & ScTn=="NE" & SpTn=="NE" & ess=="UNK")))$Gene)
temp2=subset(master_table,Gene %in% unique(subset(crypto_full_table,(CaTn=="Ess" & ScTn=="Ess" & SpTn=="Ess" & ess=="UNK")))$Gene)
discordant=rbind(temp1,temp2)
training_sub_data=training_all_data[! training_all_data$Gene %in% discordant$Gene,]


pre_rec_curve_2=data.frame(run_no =NULL, Threshold= NULL, Prec_Cons = NULL, Rec_Cons=NULL, AUROC=NULL, AUPRC=NULL)
essentials_scores_2=data.frame(Gene=NULL,essential=NULL)
features_total_2=data.frame(Parameter=NULL,MeanDecreaseAccuracy=NULL)
set.seed(100)

for(rep in 1:100){
  
  train <- sample(nrow(training_sub_data),0.8*nrow(training_sub_data),replace=FALSE)
  TrainSet <- training_sub_data[train,]
  ValidSet <- training_sub_data[-train,]
  
  
  
  
  
  TrainSet_Cons=TrainSet
  ValidSet_Cons=ValidSet
  TrainSet_Cons$conserved_essentials=as.factor(TrainSet$conserved_essentials)
  TrainSet_Cons = subset(TrainSet_Cons, select = -c(Gene) )
  rf_model_Cons_2 <- randomForest(conserved_essentials ~ ., data =TrainSet_Cons,mtry=5, ntree =10000, importance =TRUE)
  rf_model_Cons_2
  predAllprob_Cons_2 <- predict(rf_model_Cons_2, ValidSet_Cons, type="prob")
  ValidSet_Cons$pred=predAllprob_Cons_2[,2]
  
  
  
  for(cutoff in seq(0.01, 0.99, by=0.01) ){
    
    
    True_Pos_Cons=nrow(subset(ValidSet_Cons,pred>cutoff & conserved_essentials ==TRUE))
    False_Pos_Cons=nrow(subset(ValidSet_Cons,pred>cutoff & conserved_essentials ==FALSE))
    False_Neg_Cons=nrow(subset(ValidSet_Cons,pred<cutoff & conserved_essentials ==TRUE))
    
    Precision_Cons= (True_Pos_Cons)/(True_Pos_Cons+False_Pos_Cons)
    Recall_Cons= (True_Pos_Cons)/(True_Pos_Cons+False_Neg_Cons)
    
    Pos_cons=subset(ValidSet_Cons,conserved_essentials==TRUE)$pred
    Neg_cons=subset(ValidSet_Cons,conserved_essentials==FALSE)$pred
    auroc=roc.curve(scores.class0 = Pos_cons, scores.class1 = Neg_cons)$auc
    auprc=pr.curve(scores.class0 = Pos_cons, scores.class1 = Neg_cons)$auc.integral
    
    pre_rec_curve_2=rbind(pre_rec_curve_2, data.frame(run_no =rep, Threshold= cutoff, Prec_Cons = Precision_Cons, Rec_Cons=Recall_Cons, AUROC=auroc, AUPRC=auprc))
    
    
    
    
    
  }
  
  #run full predictions
  predAllprob_2 <- predict(rf_model_Cons_2, gene_stats_frame_2, type="prob")
  gene_stats_frame_2$essential=predAllprob_2[,2]
  temp_essentials=subset(gene_stats_frame_2, select= c("Gene","essential"))
  essentials_scores_2=rbind(essentials_scores_2,temp_essentials)
  
  ##features
  features=data.frame(rf_model_Cons_2$importance)
  features=cbind(rownames(features),data.frame(features,row.names=NULL))
  colnames(features)=c("Parameter","FALSE","TRUE","MeanDecreaseAccuracy","MeanDecreaseGini")
  temp_features= subset(features, select= c("Parameter","MeanDecreaseAccuracy"))
  features_total_2=rbind(features_total_2,temp_features)
  
  print(rep)
}

mean_AUROC_2=mean(pre_rec_curve_2$AUROC)
mean_AUPRC_2=mean(pre_rec_curve_2$AUPRC)

##find means for PRC curve building
prec_recall_means_2=data.frame(Threshold=NULL, Mean_Precision_Cons=NULL, Mean_Recall_Cons=NULL)
for(thresh in  levels(as.factor(pre_rec_curve_2$Threshold))){
  
  mean_prec_cons= mean(subset(pre_rec_curve_2,Threshold==thresh)$Prec_Cons)
  mean_rec_cons =mean(subset(pre_rec_curve_2,Threshold==thresh)$Rec_Cons)
  
  prec_recall_means_2=rbind(prec_recall_means_2, data.frame(Threshold=thresh,Mean_Precision_Cons=mean_prec_cons, Mean_Recall_Cons=mean_rec_cons))
  
}



baseline_2=nrow(subset(training_sub_data,conserved_essentials==TRUE))/(nrow(training_sub_data))

ggplot(prec_recall_means_2)  + geom_point(aes(x=Mean_Recall_Cons,y=Mean_Precision_Cons))+geom_hline(aes(yintercept=baseline_2),linetype="dashed",color="blue")  + theme_bw()


##find means for gene predictions
essential_score_table_2=data.frame(Gene=NULL, Essential = NULL, SD = NULL , pval= NULL)
for(genes in levels(as.factor(essentials_scores_2$Gene))){
  
  print (genes)
  
  mean_essential=mean(subset(essentials_scores_2,Gene==genes)$essential)
  sd_essential=sd(subset(essentials_scores_2,Gene==genes)$essential)
  stat_result=t.test(subset(essentials_scores_2,Gene==genes)$essential,mu=0.5,alternative="two.sided")
  
  essential_score_table_2=rbind(essential_score_table_2,data.frame(Gene=genes,Essential=mean_essential,SD=sd_essential,pval=stat_result$p.value))
}

num_genes_2=nrow(essential_score_table_2)

essential_score_table_2$ess= ifelse(essential_score_table_2$Essential<0.5 & (essential_score_table_2$pval *num_genes_2) <0.05, "NESS", ifelse(essential_score_table_2$Essential>0.5 & (essential_score_table_2$pval * num_genes_2) <0.05, "ESS","UNK"))





##find means for features

parameters_score_table_2=data.frame(Parameter=NULL, Accuracy = NULL, SD = NULL )
for(parameters in levels(as.factor(features_total_2$Parameter))){
  
  print (parameters)
  
  mean_accuracy=mean(subset(features_total_2,Parameter==parameters)$MeanDecreaseAccuracy)
  sd_accuracy=sd(subset(features_total_2,Parameter==parameters)$MeanDecreaseAccuracy)
  
  parameters_score_table_2=rbind(parameters_score_table_2,data.frame(Parameter=parameters,Accuracy=mean_accuracy,SD=sd_accuracy))
}




orthos=read.csv(file= "../../Crypto_TN_seq_paper_final_analysis/predicting_essentials/GenesByText_Summary.txt",sep="\t")
colnames(orthos)=c("Gene","Location","Description","Orthologs","Paralogs","OG_Group")

orthos=subset(orthos, Gene %in%essential_score_table$Gene)


master_table_2=merge(orthos,essential_score_table_2,by="Gene")

names(stats_frame)[names(stats_frame) == "GeneID"] = "Gene"
master_table_2 =merge(master_table_2, stats_frame, by ="Gene")




##pull in Madhani deletions
hiten_valid=read.delim("../../Crypto_TN_seq_paper_final_analysis/predicting_essentials/hiten_validated_kos.csv",sep=",")

master_table_deleted_2=subset(master_table_2, Gene %in%hiten_valid$Gene)
master_table_not_deleted_2=subset(master_table_2, !Gene %in%hiten_valid$Gene)

master_table_deleted_2$Deleted="TRUE"
master_table_not_deleted_2$Deleted="FALSE"
master_table_2=rbind(master_table_deleted_2,master_table_not_deleted_2)



##compare to other species
ortho_table=read.table(file = "../../machine_learning/GenesByOrthologs_candida.txt",header=TRUE,sep="\t")
other_essentials = read.table(file="../../machine_learning/candida_essentials.csv",header=TRUE,sep="\t")

colnames(ortho_table)=c("Gene","tran","paralog","orthocount","name","empty")

crypto_ess_table=merge(ortho_table,other_essentials,by="name")

crypto_ess_table=subset(crypto_ess_table,select=-c(tran,empty))




crypto_full_table_2 = merge(crypto_ess_table,master_table_2,by="Gene")

crypto_full_table_2 = unique(crypto_full_table_2)

##essentials_human
no_human=unique(read.table(file="crypto_WITHOUT_human_ortho.csv",header=FALSE,sep="\t"))


master_no_human_2=subset(master_table_2, Gene %in%no_human$V1)
master_in_human_2=subset(master_table_2, !Gene %in%no_human$V1)

master_no_human_2$human="FALSE"
master_in_human_2$human="TRUE"
master_table_2=rbind(master_no_human_2,master_in_human_2)



##insertion bias supplements

##Fig S1A

ggplot(subset(all_combined_frame,chr!="CP003834"))+geom_histogram(aes(x=coord))+facet_wrap(~chr,scales="free_x",ncol=7)+theme_bw()
##end Fig S1A

##for metagene plots
centered_data=data.frame(Gene=NULL,adj_freq=NULL,adj_coord=NULL,set=NULL)
length_frame=data.frame(length=NULL)
for (gene in master_table_2$Gene){
  
  cur_gene=GFF[grep(gene,GFF$V9), ] 
  
  cur_min=min(cur_gene$V4)
  cur_max=max(cur_gene$V5)
  cur_chr=head(cur_gene$V1,1)
  cur_strand=head(cur_gene$V7,1)
  cur_name=substr(head(cur_gene$V9,1),4,13)
  length_frame=rbind(length_frame,data.frame(length=(cur_max-cur_min)))
  if (cur_strand == "+"){
  
  cur_gene_up = subset(all_combined_frame, chr==cur_chr & coord < (cur_min) & coord > (cur_min-2000))
  cur_gene_up$coord=cur_gene_up$coord-cur_min
  cur_gene_up$ave_freq_T0=cur_gene_up$ave_freq_T0/sum(cur_gene_up$ave_freq_T0)
  
  cur_gene_down = subset(all_combined_frame, chr==cur_chr & coord > (cur_max) & coord < (cur_max+2000))
  cur_gene_down$coord=(cur_gene_down$coord-cur_max)+2000
  cur_gene_down$ave_freq_T0=cur_gene_down$ave_freq_T0/sum(cur_gene_down$ave_freq_T0)
  
  cur_gene_body = subset(all_combined_frame, chr==cur_chr & coord >= (cur_min) & coord <= cur_max)
  cur_gene_body$coord=((cur_gene_body$coord-cur_min)/(cur_max-cur_min))*2000
  
  cur_gene=rbind(cur_gene_down,cur_gene_up,cur_gene_body)
  
  
  }
  
  
  if (cur_strand == "-"){
    
    cur_gene_up = subset(all_combined_frame, chr==cur_chr & coord > (cur_max) & coord < (cur_max+2000))
    cur_gene_up$coord=(cur_gene_up$coord-cur_max)*-1
    cur_gene_up$ave_freq_T0=cur_gene_up$ave_freq_T0/sum(cur_gene_up$ave_freq_T0)
    
    cur_gene_down = subset(all_combined_frame, chr==cur_chr & coord < (cur_min) & coord > (cur_min-2000))
    cur_gene_down$coord=((cur_gene_down$coord-cur_min)*-1)+2000
    cur_gene_down$ave_freq_T0=cur_gene_down$ave_freq_T0/sum(cur_gene_down$ave_freq_T0)
    
    cur_gene_body = subset(all_combined_frame, chr==cur_chr & coord >= (cur_min) & coord <= cur_max)
    cur_gene_body$coord=(((cur_gene_body$coord-cur_max)*-1)/(cur_max-cur_min))*2000
    
    cur_gene=rbind(cur_gene_down,cur_gene_up,cur_gene_body)
    
    
  }
  
  
  
  if(nrow(cur_gene)>0){
    centered_data=rbind(centered_data,data.frame(Gene=gene,adj_freq=cur_gene$ave_freq_T0,adj_coord=cur_gene$coord,set="ALL"))
  }
  }
  



for (gene in subset(master_table_2,ess=="ESS")$Gene){
  
  cur_gene=GFF[grep(gene,GFF$V9), ] 
  
  cur_min=min(cur_gene$V4)
  cur_max=max(cur_gene$V5)
  cur_chr=head(cur_gene$V1,1)
  cur_strand=head(cur_gene$V7,1)
  cur_name=substr(head(cur_gene$V9,1),4,13)
  
  if (cur_strand == "+"){
    
    cur_gene_up = subset(all_combined_frame, chr==cur_chr & coord < (cur_min) & coord > (cur_min-2000))
    cur_gene_up$coord=cur_gene_up$coord-cur_min
    cur_gene_up$ave_freq_T0=cur_gene_up$ave_freq_T0/sum(cur_gene_up$ave_freq_T0)
    
    cur_gene_down = subset(all_combined_frame, chr==cur_chr & coord > (cur_max) & coord < (cur_max+2000))
    cur_gene_down$coord=(cur_gene_down$coord-cur_max)+2000
    cur_gene_down$ave_freq_T0=cur_gene_down$ave_freq_T0/sum(cur_gene_down$ave_freq_T0)
    
    cur_gene_body = subset(all_combined_frame, chr==cur_chr & coord >= (cur_min) & coord <= cur_max)
    cur_gene_body$coord=((cur_gene_body$coord-cur_min)/(cur_max-cur_min))*2000
    
    cur_gene=rbind(cur_gene_down,cur_gene_up,cur_gene_body)
    
    
  }
  
  
  if (cur_strand == "-"){
    
    cur_gene_up = subset(all_combined_frame, chr==cur_chr & coord > (cur_max) & coord < (cur_max+2000))
    cur_gene_up$coord=(cur_gene_up$coord-cur_max)*-1
    cur_gene_up$ave_freq_T0=cur_gene_up$ave_freq_T0/sum(cur_gene_up$ave_freq_T0)
    
    cur_gene_down = subset(all_combined_frame, chr==cur_chr & coord < (cur_min) & coord > (cur_min-2000))
    cur_gene_down$coord=((cur_gene_down$coord-cur_min)*-1)+2000
    cur_gene_down$ave_freq_T0=cur_gene_down$ave_freq_T0/sum(cur_gene_down$ave_freq_T0)
    
    cur_gene_body = subset(all_combined_frame, chr==cur_chr & coord >= (cur_min) & coord <= cur_max)
    cur_gene_body$coord=(((cur_gene_body$coord-cur_max)*-1)/(cur_max-cur_min))*2000
    
    cur_gene=rbind(cur_gene_down,cur_gene_up,cur_gene_body)
    
    
  }
  
  
  
  
  if(nrow(cur_gene)>0){
    centered_data=rbind(centered_data,data.frame(Gene=gene,adj_freq=cur_gene$ave_freq_T0,adj_coord=cur_gene$coord,set="Essential"))
    
  }
}

for (gene in subset(master_table_2,ess=="NESS")$Gene){
  
  cur_gene=GFF[grep(gene,GFF$V9), ] 
  
  cur_min=min(cur_gene$V4)
  cur_max=max(cur_gene$V5)
  cur_chr=head(cur_gene$V1,1)
  cur_strand=head(cur_gene$V7,1)
  cur_name=substr(head(cur_gene$V9,1),4,13)
 
  if (cur_strand == "+"){
    
    cur_gene_up = subset(all_combined_frame, chr==cur_chr & coord < (cur_min) & coord > (cur_min-2000))
    cur_gene_up$coord=cur_gene_up$coord-cur_min
    cur_gene_up$ave_freq_T0=cur_gene_up$ave_freq_T0/sum(cur_gene_up$ave_freq_T0)
    
    cur_gene_down = subset(all_combined_frame, chr==cur_chr & coord > (cur_max) & coord < (cur_max+2000))
    cur_gene_down$coord=(cur_gene_down$coord-cur_max)+2000
    cur_gene_down$ave_freq_T0=cur_gene_down$ave_freq_T0/sum(cur_gene_down$ave_freq_T0)
    
    cur_gene_body = subset(all_combined_frame, chr==cur_chr & coord >= (cur_min) & coord <= cur_max)
    cur_gene_body$coord=((cur_gene_body$coord-cur_min)/(cur_max-cur_min))*2000
    
    cur_gene=rbind(cur_gene_down,cur_gene_up,cur_gene_body)
    
    
  }
  
  
  if (cur_strand == "-"){
    
    cur_gene_up = subset(all_combined_frame, chr==cur_chr & coord > (cur_max) & coord < (cur_max+2000))
    cur_gene_up$coord=(cur_gene_up$coord-cur_max)*-1
    cur_gene_up$ave_freq_T0=cur_gene_up$ave_freq_T0/sum(cur_gene_up$ave_freq_T0)
    
    cur_gene_down = subset(all_combined_frame, chr==cur_chr & coord < (cur_min) & coord > (cur_min-2000))
    cur_gene_down$coord=((cur_gene_down$coord-cur_min)*-1)+2000
    cur_gene_down$ave_freq_T0=cur_gene_down$ave_freq_T0/sum(cur_gene_down$ave_freq_T0)
    
    cur_gene_body = subset(all_combined_frame, chr==cur_chr & coord >= (cur_min) & coord <= cur_max)
    cur_gene_body$coord=(((cur_gene_body$coord-cur_max)*-1)/(cur_max-cur_min))*2000
    
    cur_gene=rbind(cur_gene_down,cur_gene_up,cur_gene_body)
    
    
  }
  if(nrow(cur_gene)>0){
    centered_data=rbind(centered_data,data.frame(Gene=gene,adj_freq=cur_gene$ave_freq_T0,adj_coord=cur_gene$coord,set="Nonessential"))
    
  }
}


##Fig S1B
heat=ggplot(centered_data)+geom_histogram(aes(x=adj_coord,fill=set),binwidth=50)+theme_bw()+facet_wrap(~set,ncol=1,scales="free_y")
heat

##end Fig S1B

##Fig S1C
gene_stats_frame_2=gene_stats_frame_2[order(gene_stats_frame_2$Gene),]
master_table_2$Num_Inserts_per_Length=gene_stats_frame_2$Num_Inserts_per_Length
ggplot(master_table_2)+geom_boxplot(aes(x=Deleted,y=Num_Inserts_per_Length))+ scale_y_continuous(trans = "log10")

##end Fig S1C


###figure 2
###figure 2A
erg11_gene_starts= c(120725,124079,126495)
erg11_gene_ends= c(123282,126188,128651)
erg11_names= c("ZFC6","ERG11","CNAG_07308")
erg11_region= subset(all_combined_frame,chr =="CP003820" & coord < 128651 & coord > 120000)

erg11_gene_frame= data.frame(erg11_gene_starts,erg11_gene_ends,erg11_names)

##some double data points where inserts are at same location in opposite directions
##just going to drop the double points and call them out in figure legend
p1=ggplot(erg11_region)+ geom_bar(aes(x=coord,y=0.5),stat="identity")+ geom_rect(data=erg11_gene_frame,aes(ymin=-1,ymax=0,xmin=erg11_gene_starts,xmax=erg11_gene_ends,fill=erg11_names))+theme_bw()
p1



### Figure 2D
##PROC

prec_recall_means$sum=prec_recall_means$Mean_Precision_Cons+prec_recall_means$Mean_Recall_Cons
prec_recall_means <- na.omit(prec_recall_means)
max_thresh=subset(prec_recall_means_2, Threshold == 0.5)

ggplot(prec_recall_means_2)  + geom_point(aes(x=Mean_Recall_Cons,y=Mean_Precision_Cons))+geom_hline(aes(yintercept=baseline),linetype="dashed",color="blue")  + theme_bw()+
  geom_point(data=max_thresh, aes( x=Mean_Recall_Cons, y=Mean_Precision_Cons), color= "blue") + ylim(c(0,1))+geom_point(data=prec_recall_means,aes(x=Mean_Recall_Cons,y=Mean_Precision_Cons))



##Figure 2E
parameters_score_table_2$Parameter<-factor(parameters_score_table_2$Parameter, levels=c("Frequency", "Frequency_per_Length","Num_Inserts","Num_Inserts_per_Length","Middle_Frequency","Middle_Frequency_per_Length","Largest_gap","Largest_gap_perc","Chromosome","Freq_Local_100k","Length"))

ggplot(parameters_score_table_2)+geom_pointrange(aes(x=Parameter,y=Accuracy,ymin=Accuracy-SD,ymax=Accuracy+SD),color="red")+theme_bw()


##Figure 2F
histo=ggplot(essential_score_table_2)+geom_histogram(aes(x=Essential,fill=ess),alpha =0.75) + theme_bw()+geom_vline(aes(xintercept=0.5),linetype="dashed",color="orange")
histo

##essentials_human
no_human=read.table(file="../Crypto_no_human_orthos.txt",header=TRUE,sep="\t")
no_human=no_human[!duplicated(no_human[,c('Gene.ID','Ortholog.count','Paralog.count')]),]

master_no_human=subset(master_table, Gene %in%no_human$Gene.ID)
master_in_human=subset(master_table, !Gene %in%no_human$Gene.ID)

master_no_human$human="FALSE"
master_in_human$human="TRUE"
master_table=rbind(master_no_human,master_in_human)

##fig 3A

histo=ggplot(master_table_2)+geom_histogram(aes(x=Essential,fill=ess),alpha =0.75) + theme_bw()+geom_vline(aes(xintercept=0.5),linetype="dashed",color="orange")+facet_wrap(~human,ncol=1)
histo


##Figure 4B

filt_master_table=subset(master_table_2, LenSite>4)
correct_size=nrow(filt_master_table)

volc=ggplot(filt_master_table)+geom_point(aes(x=Mean50_DMSO,y=-log10(P50_DMSO*correct_size),color=(P50_DMSO*correct_size)<.05))+ scale_color_manual(values=c("#009E73","#E69F00"))+xlab("Relative proportion (Fluconazole_50/DMSO)")+ylab("-log (Corrected P-value) ")+theme_bw() 
volc

target_frame=subset(filt_master_table, Gene=="CNAG_02091" | Gene =="CNAG_01720"| Gene=="CNAG_03488" | Gene=="CNAG_04863" | Gene=="CNAG_01583"| Gene=="CNAG_03582"| Gene=="CNAG_02205"| Gene=="CNAG_05601"| Gene=="CNAG_05431")

##fcys
target_frame=subset(filt_master_table, Gene=="CNAG_01681" | Gene =="CNAG_00613")

##afr
target_frame=subset(filt_master_table, Gene=="CNAG_00730")

volc = volc +geom_point(data=target_frame, aes(x=Mean50_DMSO,y=-log10(P50_DMSO*correct_size)))+ 
  xlab("Relative proportion (Fluconazole_50/DMSO)")+ylab("-log (Corrected P-value) ")+theme_bw() + scale_color_manual(values=c("#009E73","#E69F00","red"))

volc

##Figure 4C

##gene acc numbers
## nap1   = CNAG_02091
## vps23  = CNAG_01720
## vps25  = CNAG_04863 
## snf7   = CNAG_01583
## rim20  = CNAG_03582
## rim23  = CNAG_02205
## rim13  = CNAG_05601
## rim101 = CNAG_05431
## rra1   = CNAG_03488
## afr1   = CNAG_00730

fig4c_frame=subset(all_combined_frame, Gene=="CNAG_02091" | Gene =="CNAG_01720"| Gene=="CNAG_03488" | Gene=="CNAG_04863" | Gene=="CNAG_01583"| Gene=="CNAG_03582"| Gene=="CNAG_02205"| Gene=="CNAG_05601"| Gene=="CNAG_05431" | Gene=="CNAG_00730"| Gene=="Intergenic")

ggplot(fig4c_frame)+geom_boxplot(aes(x=Gene,y=log10(ave_freq_IC50/ave_freq_DMSO))) + theme_bw()

#figure 5A

cur_gene=GFF[grep(acc_num,GFF$V9), ] 

cur_min=min(cur_gene$V4)
cur_max=max(cur_gene$V5)
cur_chr=head(cur_gene$V1,1)
cur_name=substr(head(cur_gene$V9,1),4,13)


tar_gene_frame= subset(all_combined_frame, chr==cur_chr & abs(coord)<(cur_max+300) & abs(coord)>(cur_min-300))


cds_coords=read.csv("crypto_h99_genes_whole_cds.gff",sep="\t",header=FALSE)
colnames(cds_coords)<-c("Chr","Source","Type","Start","End"," ","Strand","  ","Acc")

tar_gene = subset(cds_coords,Acc==acc_num)



coding_min= min(tar_gene$Start)
coding_max= max(tar_gene$End)

mean_pre=mean(subset(tar_gene_frame,abs(coord)<coding_max & abs(coord)>coding_min)$ave_freq_T0)
mean_post=mean(subset(tar_gene_frame,abs(coord)<coding_max & abs(coord)>coding_min)$ave_freq_IC50)

max_freq= max(c(tar_gene_frame$ave_freq_T0,tar_gene_frame$ave_freq_IC50))


p1= ggplot(tar_gene_frame)+geom_bar(aes(x=abs(coord),y=(ave_freq_T0)),stat="identity",color="#56B4E9",fill="#56B4E9")+theme_bw()+
  theme(panel.grid.minor = element_blank())+theme(legend.position = "none")+
  labs(alt="Post-Sex")+
  ylab("starting insert frequency") + theme(axis.text.x.bottom = element_blank())+theme(axis.ticks.x = element_blank())+
  xlim(c(tar_gene$Start-300,tar_gene$End+300))+ ylim(c(0,max_freq*1.1)) +theme(axis.title.x = element_blank())+
  geom_segment(x=coding_min,xend=coding_max,y=mean_pre,yend=mean_pre,color ="black",linetype=2 )+
  ggtitle("Transposon insert frequencies across transcript")+theme(plot.title = element_text(hjust = 0.5))



p2= ggplot(cur_gene)+geom_hline(yintercept=0)+geom_rect(aes(xmin=V4,xmax=V5,ymin=-1,ymax=1),fill="#009e73ff") +
  xlim(c(tar_gene$Start-300,tar_gene$End+300))+  theme_bw() +theme(axis.text.y.left = element_blank())+theme(axis.ticks = element_blank())+theme(panel.grid=element_blank())+
  theme(axis.text.x.bottom = element_blank())+theme(panel.border=element_blank())



p3=ggplot(tar_gene_frame)+geom_bar(aes(x=abs(coord),y=(-ave_freq_IC50)),stat="identity",color="#E69F00",fill="#E69F00")+theme_bw()+
  theme(panel.grid.minor = element_blank())+theme(legend.position = "none")+
  labs(alt="Post-Sex")+xlab("Genome Coordinate")+
  ylab("IC50 insert frequency") +
  xlim(c(tar_gene$Start-300,tar_gene$End+300))+  ylim(c(-max_freq*1.1,0))+
  geom_segment(x=coding_min,xend=coding_max,y=-mean_post,yend=-mean_post,color ="black",linetype=2 )

plot_grid(p1,p2,p3,ncol=1,align="v",rel_heights=c(9,2,9))






#Figure 5B

erg11_5prime=subset(all_combined_frame, chr=="CP003820" & coord>123875 & coord < 124071)
erg11_3prime=subset(all_combined_frame, chr=="CP003820" & coord>126194 & coord < 126491)

erg11_5prime$T0Norm=log10(erg11_5prime$ave_freq_T0/erg11_5prime$ave_freq_T0)
erg11_5prime$DMSONorm=log10(erg11_5prime$ave_freq_DMSO/erg11_5prime$ave_freq_T0)
erg11_5prime$IC50Norm=log10(erg11_5prime$ave_freq_IC50/erg11_5prime$ave_freq_T0)

erg11_3prime$T0Norm=log10(erg11_3prime$ave_freq_T0/erg11_3prime$ave_freq_T0)
erg11_3prime$DMSONorm=log10(erg11_3prime$ave_freq_DMSO/erg11_3prime$ave_freq_T0)
erg11_3prime$IC50Norm=log10(erg11_3prime$ave_freq_IC50/erg11_3prime$ave_freq_T0)

Initial_freq_5= data.frame(erg11_5prime$T0Norm,"Initial","5prime", erg11_5prime$strand)
DMSO_freq_5=data.frame(erg11_5prime$DMSONorm,"DMSO","5prime", erg11_5prime$strand)
IC50_freq_5=data.frame(erg11_5prime$IC50Norm,"IC50","5prime", erg11_5prime$strand)

Initial_freq_3= data.frame(erg11_3prime$T0Norm,"Initial","3prime", erg11_3prime$strand)
DMSO_freq_3=data.frame(erg11_3prime$DMSONorm,"DMSO","3prime", erg11_3prime$strand)
IC50_freq_3=data.frame(erg11_3prime$IC50Norm,"IC50","3prime", erg11_3prime$strand)


colnames(Initial_freq_5)=c("Frequency","Source","Side","Strand")
colnames(DMSO_freq_5)=c("Frequency","Source","Side","Strand")
colnames(IC50_freq_5)=c("Frequency","Source","Side","Strand")

colnames(Initial_freq_3)=c("Frequency","Source","Side","Strand")
colnames(DMSO_freq_3)=c("Frequency","Source","Side","Strand")
colnames(IC50_freq_3)=c("Frequency","Source","Side","Strand")


reg_frame=rbind(DMSO_freq_5,IC50_freq_5,DMSO_freq_3,IC50_freq_3)

reg_frame$Source=as.factor(reg_frame$Source)
reg_frame$Source=factor(reg_frame$Source,levels=c("Initial","DMSO","IC50"))


ggplot(reg_frame)+geom_boxplot(aes(x=Source,y=Frequency,fill=Side))+theme_bw()

##orientation check
ggplot(reg_frame)+geom_boxplot(aes(x=Source,y=Frequency,fill=Side))+theme_bw()+facet_wrap(~Strand,ncol=2)

##maybe add this to orientation check supp fig to show the two data types correlate pretty well
ggplot(subset(master_table_2,LenSite>20 & Len_5 >20))+geom_point(aes(x=Mean50_DMSO,y=Mean50_DMSO_5))+theme_bw()+stat_smooth(method = lm,aes(x=Mean50_DMSO,y=Mean50_DMSO_5))

model <- lm(Mean50_DMSO ~ Mean50_DMSO_5, data = subset(master_table_2,LenSite>20 & Len_5 >20))
model

summary(model)

##Figure 5C

#just essentials
filt_master_table=subset(master_table_2, Len_5>4 & ess=="ESS")
correct_size=nrow(filt_master_table)

volc=ggplot(filt_master_table)+geom_point(aes(x=Mean50_DMSO_5,y=-log10(P50_DMSO_5*correct_size),color=(P50_DMSO_5*correct_size)<.05))+ 
  scale_color_manual(values=c("#009E73","#E69F00"))+xlab("Relative proportion (Fluconazole_50/DMSO)")+ylab("-log (Corrected P-value) ")+
  theme_bw()

volc



target_frame=subset(filt_master_table, Gene== "CNAG_01137"|
                      Gene== "CNAG_01204"|
                      Gene== "CNAG_01253"|
                      Gene== "CNAG_01287"|
                      Gene== "CNAG_01345"|
                      Gene== "CNAG_01470"|
                      Gene== "CNAG_01639"|
                      Gene== "CNAG_01650"|
                      Gene== "CNAG_01789"|
                      Gene== "CNAG_01797"|
                      Gene== "CNAG_01799"|
                      Gene== "CNAG_01991"|
                      Gene== "CNAG_02791"|
                      Gene== "CNAG_02842"|
                      Gene== "CNAG_03098"|
                      Gene== "CNAG_03106"|
                      Gene== "CNAG_03165"|
                      Gene== "CNAG_03260"|
                      Gene== "CNAG_03263"|
                      Gene== "CNAG_03353"|
                      Gene== "CNAG_03544"|
                      Gene== "CNAG_03629"|
                      Gene== "CNAG_03793"|
                      Gene== "CNAG_03986"|
                      Gene== "CNAG_04000"|
                      Gene== "CNAG_04032"|
                      Gene== "CNAG_04822"|
                      Gene== "CNAG_05199"|
                      Gene== "CNAG_05909"|
                      Gene== "CNAG_06630"|
                      Gene== "CNAG_06699"| 
                      Gene== "CNAG_07363"|
                      Gene== "CNAG_00554"|
                      Gene== "CNAG_01436"|
                      Gene== "CNAG_03225"|
                      Gene== "CNAG_03358"| 
                      Gene== "CNAG_05635"|
                      Gene== "CNAG_05267"
                    
                    
)



volc = volc +geom_point(data=target_frame, aes(x=Mean50_DMSO_5,y=-log10(P50_DMSO_5*correct_size)))+ 
  xlab("Relative proportion (Fluconazole_50/DMSO)")+ylab("-log (Corrected P-value) ")+theme_bw() + scale_color_manual(values=c("#009E73","#E69F00","red"))

volc

##end Figure 5C


##Figure S3B
competition_experiment <- read_excel("G:/My Drive/Papers/Crypto TNseq paper/competition_experiment.xlsx")
ggplot(competition_experiment)+geom_segment(aes(x=1,xend=2,y=1-start,yend=1-end,color=strain))+theme_bw()+facet_wrap(~strain)+ylim(c(0,1))
##end Figure S3B

##Figure S4
flu_curve=read.table(file="../flu_conc_series.csv",sep=",")
colnames(flu_curve)=c("conc","growth")

ggplot(flu_curve)+geom_point(aes(x=conc,y=growth))+theme_bw()+geom_hline(aes(yintercept=0.5))+
  scale_x_continuous(trans = "log2")
##end Figure S4












