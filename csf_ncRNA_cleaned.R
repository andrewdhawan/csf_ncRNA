#CSF Non-coding RNA project
#Andrew Dhawan
#May 4 2020

#load necessary libraries
library(miRBaseConverter)
library(reshape2)
library(YuGene)
library(umap)
library(sva)
library(ComplexHeatmap)
library(circlize)
library(gPCA)
library(RColorBrewer)
library(pals)
library(ggridges)
library(ggpubr)

#load the datasets that we have cleaned and harmonized
# setwd("~/Google Drive/ResearchProjects/CSF_ncRNA")
load('datasets_harmonized.rda') #contains datasets_list (which has the harmonized rows with updated miR names to v22)
#load the corresponding clinical data
load('clinical_data_all.rda') 

#let's limit what we are sampling to just the common samples

#check the column names 
for(dName in names(clinical_data_list)){
	print(dName)
	common_columns = intersect(colnames(clinical_data_list[[dName]]),colnames(datasets_list[[dName]]))
	datasets_list[[dName]] = datasets_list[[dName]][,common_columns]
	clinical_data_list[[dName]] = clinical_data_list[[dName]][,common_columns]
}

#the following is to obtain key statistics about each dataset that we use for description (age, gender, etc.)
all_ages <- c()
for(dName in names(clinical_data_list)){
	print(dName) #will tell us the age distributions of the patients
	age_vec = as.numeric(as.character(clinical_data_list[[dName]]['Age',]))
	all_ages <- c(all_ages,age_vec)
	print(sum(age_vec < 20))
	print(sum((age_vec >= 20) & (age_vec < 40)))
	print(sum((age_vec >= 40) & (age_vec < 60)))
	print(sum((age_vec >= 60) & (age_vec < 80)))
	print(sum(age_vec >= 80))
	print('============')
}
all_ages = all_ages[which(!is.na(all_ages))]
print('====Overall========')
print(sum(all_ages < 20))
print(sum((all_ages >= 20) & (all_ages < 40)))
print(sum((all_ages >= 40) & (all_ages < 60)))
print(sum((all_ages >= 60) & (all_ages < 80)))
print(sum(all_ages >= 80))
print('============')


for(dName in names(clinical_data_list)){
	print(dName) #will tell us the age distributions of the patients
	gender_vec = as.character(clinical_data_list[[dName]]['Sex',])
	print(table(gender_vec))
	print('============')
}
#====================================================================================================================

#first, we select out the miRs that are well-expressed in each dataset
good_rows_all <- list() #we'll store the well-expressed miRs for each dataset here.
cutoff_value = 0.9 #means that, at worst, 90% of the values can be 0 or NA
for(dName in names(datasets_list)){
	na_vals = sum(is.na(datasets_list[[dName]])) 
	if(na_vals > 0){
		good_rows  = which(rowSums(is.na(datasets_list[[dName]])) <= (dim(datasets_list[[dName]])[2] * cutoff_value))
	}else{
		good_rows  = which(rowSums(datasets_list[[dName]]==0) <= (dim(datasets_list[[dName]])[2] * cutoff_value))
		if(dName == 'Reed_2018'){ #the Reed dataset seems to have a different 'zero' value, though does not affect analysis regardless
			good_rows  = which(rowSums(datasets_list[[dName]]==min(datasets_list[[dName]])) <= (dim(datasets_list[[dName]])[2] * cutoff_value))
		}
	}
	print(dName)
	print(length(good_rows))
	print(dim(datasets_list[[dName]])[1])
	print(length(good_rows)/dim(datasets_list[[dName]])[1]) #proportion of rows being kept after removal
	good_rows_all[[dName]] <- names(good_rows)
}

#let's limit only to the rows that are well-expressed in each sample and ensure every dataset is log-transformed
datasets_list_filtered_transformed = datasets_list
for(dName in names(datasets_list)){
	datasets_list_filtered_transformed[[dName]] = datasets_list[[dName]][good_rows_all[[dName]],]

	#if the max is higher than 50 in a given dataset, then we need to apply log transform
	if(max(datasets_list[[dName]],na.rm=T) >  50){
		datasets_list_filtered_transformed[[dName]] = log2(datasets_list[[dName]] + 1) 
	}

	#we need to remove the NA values - set these to minimum of the dataset
	datasets_list_filtered_transformed[[dName]][which(is.na(datasets_list_filtered_transformed[[dName]]),arr.ind=T)] = min(datasets_list_filtered_transformed[[dName]],na.rm=T)
}

#now we are ready to perform YuGene transform for each dataset
yugene_transformed = list()
for(dName in names(datasets_list_filtered_transformed)){
	 tmp = YuGene(data.prop = datasets_list_filtered_transformed[[dName]])
	 class(tmp) = 'matrix'
	 yugene_transformed[[dName]] = tmp
}

#set a threshold for the miRNA we are considering (this will include miRNA that are found in all datasets)
most_expressed_miRs = names(which(sort(table(melt(good_rows_all)[,1]),decreasing = T) == length(names(yugene_transformed))))

#we need a vector of all the diseases that patients have, as well as their IDs to match
all_disease_status_names <- c()
all_sample_names <- c()
for(dName in names(clinical_data_list)){
	all_sample_names = c(all_sample_names,colnames(clinical_data_list[[dName]]))
	all_disease_status_names <- c(all_disease_status_names,clinical_data_list[[dName]]['Disease',])
}
all_sample_names = melt(all_disease_status_names)[,2]
all_disease_status_names = melt(all_disease_status_names)[,1]
names(all_disease_status_names) =  all_sample_names

#the pheno variable is a construct used to help with the batch correction
pheno = matrix(0,nrow=length(all_disease_status_names),ncol=2)
colnames(pheno) = c('disease_status','batch')
row.names(pheno) = names(all_disease_status_names)
pheno[,'disease_status'] = as.character(all_disease_status_names)
counter=1
for(q in 1:length(yugene_transformed)){
	pheno[counter:(counter+(dim(yugene_transformed[[q]])[2])-1),'batch'] = as.character(q)
	counter=counter+(dim(yugene_transformed[[q]])[2])
}

pheno=as.data.frame(pheno)

combined_data <- c() 
for(n in 1:length(yugene_transformed)){
		combined_data = cbind(combined_data,yugene_transformed[[n]][most_expressed_miRs,]) #subset to only the common miRNA
}

#combat batch correction
modcombat = model.matrix(~1, data=pheno)
combat_edata = ComBat(dat=combined_data, batch=pheno$batch, mod=modcombat, par.prior=TRUE)

#examine gPCA delta value to assure that we've correctly accounted for batch effect
out_pre_correction<-gPCA.batchdetect(x=t(combined_data),batch=pheno$batch,center=F, scaleY=F,filt=NULL,nperm=1000,seed=NULL)
out_post_correction<-gPCA.batchdetect(x=t(combat_edata),batch=pheno$batch,center=F, scaleY=F,filt=NULL,nperm=1000,seed=NULL)

print(out_pre_correction$delta)
print(out_pre_correction$p.val)

print(out_post_correction$delta)
print(out_post_correction$p.val)


#make a UMAP plot of the dataset for clustering
df = t(combat_edata)
custom.config = umap.defaults
custom.config$random_state = 1234567890
custom.config$n_epochs = 2500
custom.config$n_neighbors = 30

umap.data = umap(df,custom.config)
plot.umap(x=umap.data,labels=all_disease_status_names, fName='all_patients_umap',col=stepped(n=18))

all_disease_status_names_simplified <- as.character(all_disease_status_names)
all_disease_status_names_simplified[which(all_disease_status_names_simplified %in% c('IVH','SAH'))] <- 'Hemorrhage'
all_disease_status_names_simplified[which(all_disease_status_names_simplified %in% c('AD','HD','PD','PDD','ALS'))] <- 'Degenerative'
all_disease_status_names_simplified[which(all_disease_status_names_simplified %in% c('Glioma','Ependymoma','Meningioma','Glioblastoma','Lymphoma','Medulloblastoma'))] <- 'Malignancy'
all_disease_status_names_simplified[which(all_disease_status_names_simplified %in% c('Lung_met','Breast_met'))] <- 'Metastasis'
names(all_disease_status_names_simplified) = names(all_disease_status_names)
plot.umap(x=umap.data,labels=as.factor(all_disease_status_names_simplified), fName='all_patients_umap_simplified',col=brewer.pal(n = 7, name = "Set2"))

#make a heatmap of overall expression
cols = colorRamp2(breaks=c(-0.5,0.25,1),colors=c('blue','white','red'))
disease_cols <- stepped(n=18)
dataset_cols <- brewer.pal(n = length(unique(pheno$batch)), name = "Set3")

names(disease_cols) = unique(pheno$disease_status)
names(dataset_cols) = unique(pheno$batch)

top_annotation = HeatmapAnnotation(Disease = pheno$disease_status,
show_annotation_name=c(F),
annotation_name_gp = gpar(fontsize=8),
col=list(Disease=disease_cols),
annotation_legend_param = list(title_gp=gpar(fontsize=8),
labels_gp=gpar(fontsize=7),direction = "horizontal"))
row_labels <- rownames(combat_edata)

split_cols = all_disease_status_names_simplified

dev.new()
h = Heatmap(combat_edata[,rownames(pheno)] ,
	column_title=NULL, col=cols,row_names_gp = gpar(fontsize = 7),row_names_side = "left",row_labels=row_labels, column_split=split_cols,
row_title_gp = gpar(fontsize=8),row_title_rot = 0,top_annotation=top_annotation,show_row_names=T,show_column_names=F,
 show_row_dend=F,cluster_columns= T,show_column_dend=F,cluster_rows=T,show_heatmap_legend=T,
 heatmap_legend_param = list(title = "Expr",title_position = "topleft",labels_gp=gpar(fontsize=7),
 title_gp=gpar(fontsize=8),direction = "vertical"))
draw(h,merge_legend=T,heatmap_legend_side = "right",  annotation_legend_side = "right")
dev.copy(pdf,file='overall_summary_matrix.pdf',width=8,height=4)
dev.off()
#========================================================================================================


#now, let's perform differential miRNA expression analysis using Wilcoxon across diseases
healthy_pts = rownames(pheno)[which(pheno$disease_status == 'Healthy')]

all_p_vals_diseases <- list()
for(disease_name in unique(pheno$disease_status)){
	print(disease_name)
	if(disease_name!="Healthy"){
		all_p_vals <- c()
		disease_samples = rownames(pheno)[which(pheno$disease_status == disease_name)]
		for(miR_name in rownames(combat_edata)){
			tmp_greater = wilcox.test(combat_edata[miR_name,disease_samples], combat_edata[miR_name,healthy_pts],alternative='greater')$p.value
			tmp_less = wilcox.test(combat_edata[miR_name,disease_samples], combat_edata[miR_name,healthy_pts],alternative='less')$p.value
			if(tmp_greater < tmp_less){
				all_p_vals <- c(all_p_vals,tmp_greater)		
			}else{
				all_p_vals <- c(all_p_vals,(-1) * tmp_less)	
			}
		}
		names(all_p_vals) = rownames(combat_edata)
		all_p_vals = sort(all_p_vals)
		print(names(all_p_vals)[which(abs(all_p_vals) < (0.05 / length(all_p_vals)))])
		all_p_vals_diseases[[disease_name]] = all_p_vals
	}
}


all_p_vals_diseases_simplified <- list()
for(disease_name in unique(all_disease_status_names_simplified)){
	print(disease_name)
	if(disease_name!="Healthy"){
		all_p_vals <- c()
		disease_samples = names(all_disease_status_names_simplified)[which(all_disease_status_names_simplified == disease_name)]
		for(miR_name in rownames(combat_edata)){
			tmp_greater = wilcox.test(combat_edata[miR_name,disease_samples], combat_edata[miR_name,healthy_pts],alternative='greater')$p.value
			tmp_less = wilcox.test(combat_edata[miR_name,disease_samples], combat_edata[miR_name,healthy_pts],alternative='less')$p.value
			if(tmp_greater < tmp_less){
				all_p_vals <- c(all_p_vals,tmp_greater)		
			}else{
				all_p_vals <- c(all_p_vals,(-1) * tmp_less)	
			}
		}
		names(all_p_vals) = rownames(combat_edata)
		all_p_vals = sort(all_p_vals)
		print(names(all_p_vals)[which(abs(all_p_vals) < (0.05 / length(all_p_vals)))])
		all_p_vals_diseases_simplified[[disease_name]] = all_p_vals
	}
}


melted_df <- melt(combat_edata)
melted_df$disease_name = pheno$disease_status[melted_df$Var2]
melted_df$disease_name_simplified = all_disease_status_names_simplified[melted_df$Var2]
dev.new()
ggboxplot(melted_df[which(melted_df$Var1 == 'hsa-miR-767-5p' & (melted_df$disease_name %in% c('AD','Healthy'))),],x='disease_name',y='value',color='black',palette='npj',add='jitter',add.params=list(alpha=0.2),xlab='',ylab='Normalized Expr',alpha=0.2,title='hsa-miR-767-5p') + rotate_x_text(angle=30)+  stat_compare_means() 
dev.copy(pdf,file='healthy_ad_mir767.pdf',width=3,height=4)
dev.off()

dev.new()
ggboxplot(melted_df[which(melted_df$Var1 == 'hsa-miR-361-3p' & (melted_df$disease_name %in% c('HD','Healthy'))),],x='disease_name',y='value',color='black',palette='npj',add='jitter',add.params=list(alpha=0.2),xlab='',ylab='',alpha=0.2,title='hsa-miR-361-3p') + rotate_x_text(angle=30)+  stat_compare_means()  
dev.copy(pdf,file='healthy_hd_mir361.pdf',width=3,height=4)
dev.off()

dev.new()
ggboxplot(melted_df[which(melted_df$Var1 == 'hsa-miR-885-5p' & (melted_df$disease_name %in% c('HD-Pre-Low','Healthy'))),],x='disease_name',y='value',color='black',palette='npj',add='jitter',add.params=list(alpha=0.2),xlab='',ylab='',alpha=0.2,title='hsa-miR-885-5p') +  rotate_x_text(angle=30)+ stat_compare_means()  
dev.copy(pdf,file='healthy_hd_pre_low_mir885.pdf',width=3,height=4)
dev.off()

dev.new()
ggboxplot(melted_df[which(melted_df$Var1 == 'hsa-miR-142-5p' & (melted_df$disease_name %in% c('ALS','Healthy'))),],x='disease_name',y='value',color='black',palette='npj',add='jitter',add.params=list(alpha=0.2),xlab='',ylab='',alpha=0.2,title='hsa-miR-142-5p') + rotate_x_text(angle=30)+  stat_compare_means()  
dev.copy(pdf,file='healthy_als_mir142.pdf',width=3,height=4)
dev.off()


dev.new()
p1 = ggboxplot(melted_df[which(melted_df$Var1 == 'hsa-miR-142-5p' & (melted_df$disease_name_simplified %in% c('Degenerative','Healthy'))),],x='disease_name_simplified',y='value',color='black',palette='npj',add='jitter',add.params=list(alpha=0.2),xlab='',ylab='',alpha=0.2,title='hsa-miR-142-5p') + rotate_x_text(angle=30)+  stat_compare_means()  
p2 = ggboxplot(melted_df[which(melted_df$Var1 == 'hsa-miR-361-3p' & (melted_df$disease_name_simplified %in% c('Degenerative','Healthy'))),],x='disease_name_simplified',y='value',color='black',palette='npj',add='jitter',add.params=list(alpha=0.2),xlab='',ylab='',alpha=0.2,title='hsa-miR-361-3p') + rotate_x_text(angle=30)+  stat_compare_means()  
ggarrange(p1,p2,nrow=1)
dev.copy(pdf,file='healthy_degenerative_mir142_mir361.pdf',width=6,height=4)
dev.off()


#=====================================================================================================

#-----Here, we identify the miRNA that are highly expressed consistently in the healthy state.

#Rank Product computation
RP.out <- RP(combat_edata[,healthy_pts],rep(1,length(healthy_pts)))
RP_out_values <- RP.out
rank_prod_tables <- topGene(RP.out,cutoff = 0.05,method="pfp",gene.names=rownames(combat_edata))

#Plot the results
plot_df = melted_df[which((melted_df$disease_name == 'Healthy') & (melted_df$Var1 %in% rownames(rank_prod_tables$Table2))),]
plot_df$text <-  -log10(rank_prod_tables$Table2[as.character(plot_df$Var1),5])
plot_df$text = format(plot_df$text,digits=2)
plot_df$Var1 = ordered(plot_df$Var1,levels=rev(rownames(rank_prod_tables$Table2)))
dev.new()
ggplot(plot_df, aes(x = value, y = Var1, fill = Var1)) +
scale_fill_brewer(palette="Pastel2") + 
geom_text(aes(label =text, x = -0.4),fontface = "plain",family = c("sans"),size=4,
    # position = position_dodge(0.9),
    vjust = 0) + 
  geom_density_ridges() +
  theme_ridges() + 
  labs(x=NULL,y=NULL)+
  theme(legend.position = "none")	
  dev.copy(pdf,file='joyplot_miRNA_expr_normal.pdf',width=5,height=4)
  dev.off()

#==================================================================================================

# Create Donut plot for the samples types present (Descriptive statistic)
plot_data <- data.frame(
  category=names(table(all_disease_status_names)),
  count=as.vector(table(all_disease_status_names))
)
 
# Compute percentages
plot_data$fraction <- plot_data$count / sum(plot_data$count)
# Compute the cumulative percentages (top of each rectangle)
plot_data$ymax <- cumsum(plot_data$fraction)
# Compute the bottom of each rectangle
plot_data$ymin <- c(0, head(plot_data$ymax, n=-1))
# Compute label position
plot_data$labelPosition <- (plot_data$ymax + plot_data$ymin) / 2
# Compute a good label
plot_data$label <- paste0(plot_data$category, "\n value: ", plot_data$count)
dev.new()
# Make the plot
ggplot(plot_data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_text( x=2, aes(y=labelPosition, label=label, color=category), size=6) + # x here controls label position (inner / outer)
  scale_fill_brewer(palette='Set3') +
  scale_color_brewer(palette='Set3') +
  coord_polar(theta="y") +
  xlim(c(0, 4)) +
  theme_void() +
  theme(legend.position = "none")
dev.copy(pdf,file='donut_plot_all_samples_2.pdf',height=4,width=4)
dev.off()
#============================================================================================


#---------Now, we examine co-expression networks among subgroups:
all_cor_mats <- list() #list of matrices of correlation coefficients
all_cor_mats_p_values <- list() #list of matries for the p values of the correlation coefficients
for(disease_name in unique(all_disease_status_names_simplified)){
	tmp_mat <- matrix(NA,nrow=length(rownames(combat_edata)),ncol=length(rownames(combat_edata)))
	row.names(tmp_mat) <- rownames(combat_edata)
	colnames(tmp_mat) <- rownames(combat_edata)
	tmp_mat_p_val <- tmp_mat
	desired_columns = names(all_disease_status_names_simplified)[which(all_disease_status_names_simplified == disease_name)]

	for(miR_name_1 in rownames(combat_edata)){
		for(miR_name_2 in rownames(combat_edata)){
			cor_val = cor.test(as.numeric(combat_edata[miR_name_1,desired_columns]),as.numeric(combat_edata[miR_name_2,desired_columns]),method='spearman')
			tmp_mat[miR_name_1,miR_name_2] = cor_val$estimate
			tmp_mat_p_val[miR_name_1,miR_name_2] = cor_val$p.value
		}
	}
	all_cor_mats[[disease_name]] <- tmp_mat
	all_cor_mats_p_values[[disease_name]] <- tmp_mat_p_val
}



all_cor_mats_filtered <- list()
cols = colorRamp2(breaks=c(-1,0,1),colors=c('blue','white','red'))
hmaps <- c()
counter = 1

#this is just to get the plotting information from Hemorrhage heatmap
disease_name='Hemorrhage'
tmp_mat <- all_cor_mats[[disease_name]]
tmp_mat[which(is.na(tmp_mat),arr.ind=T)] <- 0
tmp_mat_p_val <- all_cor_mats_p_values[[disease_name]]
tmp_mat_p_val[which(is.na(tmp_mat_p_val),arr.ind=T)] <- 1
tmp_mat[which(tmp_mat_p_val >= 0.05,arr.ind=T)] <- 0
all_cor_mats_filtered[[disease_name]] <- tmp_mat
print(head(tmp_mat))
show_leg = (disease_name == 'Malignancy')
tmp_hmap = Heatmap(tmp_mat,col=cols,show_row_dend=F,
	show_column_dend=F, row_names_gp=gpar(fontsize=5),column_title=disease_name,column_order=col_order,
	row_names_side = "left",show_column_names=F, cluster_columns=T,
	show_heatmap_legend=show_leg,
	heatmap_legend_param = list(title = "Spearman's rho",
	title_position = "topleft",labels_gp=gpar(fontsize=7),
	title_gp=gpar(fontsize=8),direction = "vertical"))
col_order=column_order(tmp_hmap) 
#================---------------------======================


for(disease_name in unique(all_disease_status_names_simplified)){
	if(table(all_disease_status_names_simplified)[disease_name]> 24){
		tmp_mat <- all_cor_mats[[disease_name]]
		tmp_mat[which(is.na(tmp_mat),arr.ind=T)] <- 0
		tmp_mat_p_val <- all_cor_mats_p_values[[disease_name]]
		tmp_mat_p_val[which(is.na(tmp_mat_p_val),arr.ind=T)] <- 1
		tmp_mat[which(tmp_mat_p_val >= 0.05,arr.ind=T)] <- 0
		all_cor_mats_filtered[[disease_name]] <- tmp_mat
		print(head(tmp_mat))
		show_leg = (disease_name == 'Malignancy')

		hmaps <- hmaps + Heatmap(tmp_mat,col=cols,show_row_dend=F,
			show_column_dend=F, row_names_gp=gpar(fontsize=6),column_title=disease_name,column_order=col_order,row_order=col_order,
			row_names_side = "left",show_column_names=F, cluster_columns=F,cluster_rows=F,
			show_heatmap_legend=show_leg,
			heatmap_legend_param = list(title = "Spearman's rho",
			title_position = "topleft",labels_gp=gpar(fontsize=7),
			title_gp=gpar(fontsize=8),direction = "vertical"))

		if(disease_name == 'Malignancy'){
			plot_counter = counter
		}
		counter = counter + 1
	}	
}

dev.new()
draw(hmaps,main_heatmap = plot_counter)
dev.copy(pdf,file=paste0('heatmap_correlation_all','.pdf'),height=3,width=12) #actual correlation heatmaps being plotted
dev.off()	

#===========================================================================
#-----Now, identify those miRNA that are correlated in opposing directions across the disease groups

inversely_correlated_miRs <- c()
for(disease_name_1 in unique(all_disease_status_names_simplified)){
	for(disease_name_2 in unique(all_disease_status_names_simplified)){
		if(disease_name_1 != disease_name_2){

			inversely_correlated_cases <- which((all_cor_mats_filtered[[disease_name_1]] *  all_cor_mats_filtered[[disease_name_2]]) < 0 ,arr.ind=T)
			if(length(inversely_correlated_cases) > 0){
				values_to_rank_by = abs(all_cor_mats_filtered[[disease_name_1]][inversely_correlated_cases]) + abs(all_cor_mats_filtered[[disease_name_2]][inversely_correlated_cases])
				inversely_correlated_miRs$miR_1 = c(inversely_correlated_miRs$miR_1,rownames(all_cor_mats_filtered[[disease_name_1]])[inversely_correlated_cases[,1]])
				inversely_correlated_miRs$miR_2 = c(inversely_correlated_miRs$miR_2,colnames(all_cor_mats_filtered[[disease_name_1]])[inversely_correlated_cases[,2]])
				inversely_correlated_miRs$disease_1 = c(inversely_correlated_miRs$disease_1,rep(disease_name_1,times=dim(inversely_correlated_cases)[1]))
				inversely_correlated_miRs$disease_2 = c(inversely_correlated_miRs$disease_2,rep(disease_name_2,times=dim(inversely_correlated_cases)[1]))
				inversely_correlated_miRs$values = c(inversely_correlated_miRs$values, values_to_rank_by) 
			}
		}
	}
}

inversely_correlated_miRs = as.data.frame(inversely_correlated_miRs)
switch_rows <- which(as.character(inversely_correlated_miRs$miR_1) < as.character(inversely_correlated_miRs$miR_2))
tmp = inversely_correlated_miRs$miR_2[switch_rows]
inversely_correlated_miRs$miR_2[switch_rows] = inversely_correlated_miRs$miR_1[switch_rows]
inversely_correlated_miRs$miR_1[switch_rows] = tmp

switch_rows <- which(as.character(inversely_correlated_miRs$disease_1) < as.character(inversely_correlated_miRs$disease_2))
tmp = inversely_correlated_miRs$disease_2[switch_rows]
inversely_correlated_miRs$disease_2[switch_rows] = inversely_correlated_miRs$disease_1[switch_rows]
inversely_correlated_miRs$disease_1[switch_rows] = tmp

inversely_correlated_miRs <- unique(inversely_correlated_miRs)
inversely_correlated_miRs[order(inversely_correlated_miRs$values,decreasing=T),]

#now plotting commands for a specific case

malig_cols = names(all_disease_status_names_simplified)[which(all_disease_status_names_simplified == 'Malignancy')]
hem_cols = names(all_disease_status_names_simplified)[which(all_disease_status_names_simplified == 'Hemorrhage')]
deg_cols = names(all_disease_status_names_simplified)[which(all_disease_status_names_simplified == 'Degenerative')]
healthy_cols = names(all_disease_status_names_simplified)[which(all_disease_status_names_simplified == 'Healthy')]

data_11 <- combat_edata['hsa-miR-769-5p',malig_cols]
data_12 <- combat_edata['hsa-miR-127-3p',malig_cols]
all_cor_mats_filtered[['Malignancy']]['hsa-miR-769-5p','hsa-miR-127-3p']
all_cor_mats_p_values[['Malignancy']]['hsa-miR-769-5p','hsa-miR-127-3p']

data_21 <- combat_edata['hsa-miR-769-5p',hem_cols]
data_22 <- combat_edata['hsa-miR-127-3p',hem_cols]
all_cor_mats_filtered[['Hemorrhage']]['hsa-miR-769-5p','hsa-miR-127-3p']
all_cor_mats_p_values[['Hemorrhage']]['hsa-miR-769-5p','hsa-miR-127-3p']


data_31 <- combat_edata['hsa-miR-769-5p',deg_cols]
data_32 <- combat_edata['hsa-miR-127-3p',deg_cols]
all_cor_mats_filtered[['Degenerative']]['hsa-miR-769-5p','hsa-miR-127-3p']
all_cor_mats_p_values[['Degenerative']]['hsa-miR-769-5p','hsa-miR-127-3p']


data_41 <- combat_edata['hsa-miR-769-5p',healthy_cols]
data_42 <- combat_edata['hsa-miR-127-3p',healthy_cols]
all_cor_mats_filtered[['Healthy']]['hsa-miR-769-5p','hsa-miR-127-3p']
all_cor_mats_p_values[['Healthy']]['hsa-miR-769-5p','hsa-miR-127-3p']


dev.new()
plot(data_11,data_12,col=alpha('black',0.2),pch=20,xlab='',ylab='',xlim=c(-0.5,1.5),ylim=c(-0.5,1.5),main=NULL,frame.plot=F)
dev.copy(pdf,file='malignant_769_127.pdf',width=3,height=3)
dev.off()

dev.new()
plot(data_21,data_22,col=alpha('black',0.2),pch=20,xlab='',ylab='',xlim=c(-0.5,1.5),ylim=c(-0.5,1.5),main=NULL,frame.plot=F)
dev.copy(pdf,file='hem_769_127.pdf',width=3,height=3)
dev.off()

dev.new()
plot(data_31,data_32,col=alpha('black',0.2),pch=20,xlab='',ylab='',xlim=c(-0.5,1.5),ylim=c(-0.5,1.5),main=NULL,frame.plot=F)
dev.copy(pdf,file='deg_769_127.pdf',width=3,height=3)
dev.off()

dev.new()
plot(data_41,data_42,col=alpha('black',0.2),pch=20,xlab='',ylab='',xlim=c(-0.5,1.5),ylim=c(-0.5,1.5),main=NULL,frame.plot=F)
dev.copy(pdf,file='healthy_769_127.pdf',width=3,height=3)
dev.off()

#======================================================================================================================


#-----helper functions----

plot.umap <- function(x, labels, fName,
         main="",
         colors=rainbow(20),
         pad=0.1, cex=0.65, pch=19, add=FALSE, show_leg=T,legend.suffix="",
         cex.main=1, cex.legend=1) {
#adapted from plot.iris

  layout = x
  if (is(x, "umap")) {
    layout = x$layout
  } 
  
  xylim = range(layout)
  xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
	dev.new()
  if (!add) {
    par(mar=c(0.2,0.9,1.2,0.9), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F,main=main)
  }
  points(layout[,1], layout[,2], col=colors[as.integer(labels)],
         cex=cex, pch=pch)

  labels.u = unique(labels)
  legend.pos = "topright"
  legend.text = as.character(labels.u)
  if (add) {
    legend.pos = "bottomright"
    legend.text = paste(as.character(labels.u), legend.suffix)
  }
  if(show_leg){
	  legend(legend.pos, legend=legend.text,
         col=colors[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)

  }
  dev.copy(pdf,file=paste0(fName,'.pdf'),height=7,width=7)
  dev.off()
}


#----Repeat the above analysis but with a dataset removed------

all_disease_status_names_all <- all_disease_status_names
all_disease_status_names_simplified_all <- all_disease_status_names_simplified

for(dName2 in names(yugene_transformed)){
	new_names <- setdiff(names(yugene_transformed),dName2)
	yugene_transformed_one_removed <- yugene_transformed[new_names]
	clinical_data_one_removed <- clinical_data_list[new_names]

	#redo the analysis here
	good_rows_all_one_removed <- good_rows_all[new_names]
	most_expressed_miRs = names(which(sort(table(melt(good_rows_all_one_removed)[,1]),decreasing = T) == length(names(yugene_transformed_one_removed))))

	#we need a vector of all the diseases that patients have, as well as their IDs to match
	all_disease_status_names <- c()
	all_sample_names <- c()
	for(dName in names(clinical_data_one_removed)){
		all_sample_names = c(all_sample_names,colnames(clinical_data_one_removed[[dName]]))
		all_disease_status_names <- c(all_disease_status_names,clinical_data_one_removed[[dName]]['Disease',])
	}
	all_sample_names = melt(all_disease_status_names)[,2]
	all_disease_status_names = melt(all_disease_status_names)[,1]
	names(all_disease_status_names) =  all_sample_names

	pheno = matrix(0,nrow=length(all_disease_status_names),ncol=2)
	colnames(pheno) = c('disease_status','batch')
	row.names(pheno) = names(all_disease_status_names)
	pheno[,'disease_status'] = as.character(all_disease_status_names)
	counter=1
	for(q in 1:length(yugene_transformed_one_removed)){
		pheno[counter:(counter+(dim(yugene_transformed_one_removed[[q]])[2])-1),'batch'] = as.character(q)
		counter=counter+(dim(yugene_transformed_one_removed[[q]])[2])
	}

	pheno=as.data.frame(pheno)

	combined_data <- c()
	for(n in 1:length(yugene_transformed_one_removed)){
		combined_data = cbind(combined_data,yugene_transformed_one_removed[[n]][most_expressed_miRs,])
	}

	modcombat = model.matrix(~1, data=pheno)
	combat_edata = ComBat(dat=combined_data, batch=pheno$batch, mod=modcombat, par.prior=TRUE)

	df = t(combat_edata)
	custom.config = umap.defaults
	custom.config$random_state = 1234567890
	custom.config$n_epochs = 2500
	custom.config$n_neighbors = 30

	umap.data = umap(df,custom.config)
	plot.umap(x=umap.data,labels=all_disease_status_names_all, fName=paste0('all_patients_umap_',dName2,'_removed'),col=stepped(n=18),main=paste0(dName2,' Removed'),show_leg=F)
	plot.umap(x=umap.data,labels=as.factor(all_disease_status_names_simplified_all), fName=paste0('all_patients_umap_simplified_',dName2,'_removed'),col=brewer.pal(n = 7, name = "Set2"),main=paste0(dName2,' Removed'),show_leg=F)
}