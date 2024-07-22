update.packages(checkBuilt=TRUE, ask=FALSE)

setwd("/Volumes/mum2023/home/z6/Documents/workingdir/Final-nasalmicrob/2022/")
{#library 
library("knitr")
library("ggplot2")
library("microbiomeSeq")
library("vegan")
library("plyr")
library("rstatix")
library("dplyr")
library("qiime2R")
library("phyloseq")
library("microbiome")
library("ggpubr")
}
library(devtools)
install_github("umerijaz/microbiomeSeq") 
library(microbiomeSeq)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("AnnotationDbi", "DESeq2", "GO.db", "impute", "phyloseq", "preprocessCore"))

#datainfo
summarize_phyloseq(ps1_unq)
df <- as.data.frame(sample_sums(ps1_unq))
sd(df$`sample_sums(ps1_unq)`)

{#####read table of seq variants/ import data qiime to r
ps<-qza_to_phyloseq(
  features = "ubms-v4fwd-fmitochlo_table.qza", 
  tree="ubms-v4fwd_insertion-tree.qza",
  taxonomy="ubms-v4fwd_taxonomy.qza",
  #otu table
  metadata = "metadata.tsv"
)
#####filtering and QC
table(tax_table(ps)[, "Phylum"], exclude = NULL) #check uncharacterized phylum
#filter out samples with less than 10,000 ASV
ps0 <- prune_samples(sample_sums(ps) >= 10000, ps) #reads less than 10,000
ps1 = subset_taxa(ps0, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))#filter un-characterized phylum 


#unique samples
meta(ps1_unq)
ps1_unq <- subset_samples(ps1, unique == "Y")
filter <- phyloseq::genefilter_sample(ps1_unq, filterfun_sample(function(x) x >= 100), A = 2)
psprev_unq <- prune_taxa(filter,ps1_unq)
psprev_unq
}

################KEY TAXA INTERESTED= 130
taxaname_psprev_unq <- psprev_unq
taxa_names(taxaname_psprev_unq) <- paste0("ASV",seq(ntaxa(taxaname_psprev_unq))) #add asv id
overall130 <- as.data.frame(tax_table(taxaname_psprev_unq))
write.table(overall130,"figures/overall130.csv")

ps1_prev_id <- taxa_names(psprev_unq) #get taxa ID the taxa interested
check <- aggregate_taxa(ps, "Genus")
table = data.frame(tax_table(psprev_unq),
                   Abundance = abundances(psprev_unq)*100)

##sparsity of microbiome data (%of zeros)
sum(otu_table(ps1_unq) == 0) /(6141*87)#dim(otu_table(ps1_unq))
#ps1_prev <- filter_taxa(ps1, function(x){sum(x > 0) > 5}, prune = TRUE) #prevalence in 5samples
#filter <- phyloseq::genefilter_sample(ps1, filterfun_sample(function(x) x >= 50), A = 5)
rarecurve(t(otu_table(psprev_unq)),step=500, cex=.75,las=1,xlim=c(-50,100000))


##############################################
#######################Compositional
#compare between OA groups
psprev_unq_phyla_JH <- prune_samples(psprev_unq, group == "Jehai_FC")
psprev_unq_phyla_JH <- tax_glom(psprev_unq_phyla_JH, taxrank = 'Phylum') # agglomerate taxa
phyla.sum = tapply(taxa_sums(psprev_unq_phyla_JH), tax_table(psprev_unq_phyla_JH)[, "Phylum"], sum)
sort(phyla.sum, decreasing=TRUE)
Majorcomp <- data.frame(tax_table(psprev_unq_phyla_JH),Prevalence = prevalence(psprev_unq_phyla_JH,count=TRUE), # Add taxonomy and total read counts to this data.frame
                        TotalAbundance = taxa_sums(psprev_unq_phyla_JH), 
                        Percent=(taxa_sums(psprev_unq_phyla_JH)/1632905)*100
                        )
sum(Majorcomp$TotalAbundance)
#OVERALL
psprev_unq_phyla <- tax_glom(psprev_unq, taxrank = 'Phylum') # agglomerate taxa
phyla.sum = tapply(taxa_sums(psprev_unq_phyla), tax_table(psprev_unq_phyla)[, "Phylum"], sum)
sort(phyla.sum, decreasing=TRUE)
Majorcomp <- data.frame(tax_table(psprev_unq_phyla),Prevalence = prevalence(psprev_unq_phyla,count=TRUE), # Add taxonomy and total read counts to this data.frame
                        TotalAbundance = taxa_sums(psprev_unq_phyla), 
                        Percent=(taxa_sums(psprev_unq_phyla)/sum(taxa_sums(psprev_unq_phyla)))*100)
Majorcomp
write.table(Majorcomp,"figures/majorcomp.csv")
psprev_unq_phylabar <- psmelt(psprev_unq_phyla)
psprev_unq_phylabar$Phylum[psprev_unq_phylabar$Phylum 
                        %in% c( "Halobacterota", "Campilobacterota", "Patescibacteria", "Fusobacteriota", "Basidiomycota")] <- "Other Phyla"
topphy = c("Actinobacteriota","Firmicutes","Proteobacteria","Bacteroidota","Other Phyla")
psprev_unq_phylabar$Phylum<- factor(psprev_unq_phylabar$Phylum, levels = topphy)
p <- ggplot(data=psprev_unq_phylabar, aes(x=participant, y=Abundance, fill= Phylum))
phylum <- c("darkgreen", "deeppink","white" ,"yellow","black","grey","orange","yellow",
                   "blue", "darkslategray1","purple", "red","black","cornsilk",
                   "black","#679289","darkolivegreen1","black","#F9ECCC" )
p + geom_bar(stat="identity", position = "fill", color= "black", width=1.0) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(ncol = 5))+
  scale_fill_manual(name="Phylum",values = phylum)+
  facet_grid(cols = vars(group), space="free_x", scales="free_x", switch="x", drop = FALSE)+
  theme(legend.position = "bottom", axis.text.y.left = element_text(size = 14), 
        axis.text.x = element_text(size =11, angle=90, vjust = 0.5),
        strip.text.x = element_text(size=14),legend.text = element_text(size=12),
        axis.title.y = element_text(size = 14),axis.title.x = element_blank())
ggsave("figures/taxaplot_overall.jpeg", width=40, height=25, units = "cm")


######BARPLOT BASED ON CLASSES
psprev_unq_glom <- tax_glom(psprev_unq, taxrank = 'Class') # agglomerate taxa
class.sum = tapply(taxa_sums(psprev_unq_glom), tax_table(psprev_unq_glom)[, "Class"], sum)
sort(class.sum, decreasing=TRUE)
psprev_unq_tmglom <- psmelt(psprev_unq_glom) # create dataframe from phyloseq object
#rename genera with < 1% abundance
psprev_unq_tmglom$Class[psprev_unq_tmglom$Class %in% c("Gracilibacteria","Fusobacteriia","Malasseziomycetes","Saccharimonadia")] <- "Classes < 0.05% abund."
topclass = c("Actinobacteria","Bacilli","Gammaproteobacteria","Clostridia","Bacteroidia",
             "Negativicutes","Halobacteria","Campylobacteria","Alphaproteobacteria","Classes < 0.05% abund.")
psprev_unq_tmglom$Class<- factor(psprev_unq_tmglom$Class, levels = topclass)
#plot
p <- ggplot(data=psprev_unq_tmglom, aes(x=participant, y=Abundance, fill=Class))
phylum_colors <- c("aquamarine3", "lightpink", "dark grey","orange","yellow",
                   "blue", "darkslategray1","purple", "red","black","cornsilk",
                   "black","#679289","darkolivegreen1","black","#F9ECCC" )
p + geom_bar(stat="identity", position = "fill", color= "black", width=1.0) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(ncol = 5))+
  scale_fill_manual(name="Class",values = phylum_colors)+
  facet_grid(cols = vars(group), space="free_x", scales="free_x", switch="x", drop = FALSE)+
  theme(legend.position = "bottom", axis.text.y.left = element_text(size = 14), 
        axis.text.x = element_text(size =11, angle=90, vjust = 0.5),
        strip.text.x = element_text(size=14),legend.text = element_text(size=12),
        axis.title.y = element_text(size = 14),axis.title.x = element_blank())
ggsave("figures/taxaplot_group.jpeg", width=40, height=25, units = "cm")

#plot based on gender
p <- ggplot(data=psprev_unq_tmglom, aes(x=participant, y=Abundance, fill=Class))
p + geom_bar(stat="identity", position = "fill", color= "black", width=1.0) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(ncol = 5))+
  scale_fill_manual(name="Class",values = phylum_colors)+
  facet_grid(cols = vars(gender), space="free_x", scales="free_x", switch="x", drop = FALSE)+
  theme(legend.position = "bottom", axis.text.y.left = element_text(size = 14), 
        axis.text.x = element_text(size =11, angle=90, vjust = 0.5),
        strip.text.x = element_text(size=14),legend.text = element_text(size=12),
        axis.title.y = element_text(size = 14),axis.title.x = element_blank())
ggsave("figures/taxaplot_gender.jpeg", width=40, height=25, units = "cm")
#plot based on waist
p <- ggplot(data=psprev_unq_tmglom, aes(x=participant, y=Abundance, fill=Class))
p + geom_bar(stat="identity", position = "fill", color= "black", width=1.0) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(ncol = 5))+
  scale_fill_manual(name="Class",values = phylum_colors)+
  facet_grid(cols = vars(waistCM), space="free_x", scales="free_x", switch="x", drop = FALSE)+
  theme(legend.position = "bottom", axis.text.y.left = element_text(size = 14), 
        axis.text.x = element_text(size =11, angle=90, vjust = 0.5),
        strip.text.x = element_text(size=14),legend.text = element_text(size=12),
        axis.title.y = element_text(size = 14),axis.title.x = element_blank())
ggsave("figures/taxaplot_waist.jpeg", width=40, height=25, units = "cm")

###master taxa name, rare classes
taxa_names(taxaname_psprev_unq) #get taxa ID the taxa interested
taxaname_psprev_unq <- subset_taxa(taxaname_psprev_unq, Class %in% c("Gracilibacteria","Fusobacteriia","Malasseziomycetes","Saccharimonadia"))
rare <- data.frame(tax_table(taxaname_psprev_unq),Prevalence = prevalence(taxaname_psprev_unq,count=TRUE), # Add taxonomy and total read counts to this data.frame
                   TotalAbundance = taxa_sums(taxaname_psprev_unq))
write.csv(rare,"figures/rareclasses.csv")


#######################################################
#checking specific taxa
taxa_test <- tax_glom(psprev_unq, taxrank = 'Family') 
taxa_test <- microbiome::transform(taxa_test,"compositional")
test <- subset_taxa(taxa_test, Family == "Neisseriaceae")
title = "plot_bar; Bacteroidetes-only"
plot_bar(test, "participant", "Abundance", title=title,fill="Family")
#checking phylum
tax_table(test)
lawsonella_test <- tax_glom(psprev_unq, taxrank = 'Kingdom') 
test = subset_taxa(lawsonella_test, Kingdom == "d__Archaea")
title = "plot_bar; d__Archaea-only"
plot_bar(test, "participant", "Abundance", fill= "Kingdom",title=title)+
  facet_grid(cols = vars(group), space="free_x", scales="free_x", switch="x", drop = FALSE)
#######################################################

#taxa_names(ps_comp) <- paste0("ASV",seq(ntaxa(ps_comp))) #add asv id
!!!!!!!!!!!!!!check tax glom
tax_glom(psprev_unq, taxrank = 'Family') 
######Core with compositionals and Heatmap
ps_comp <-  microbiome::transform(psprev_unq, "compositional") #transform to relative abundance
ps_comp <- merge_taxa(ps_comp,"Genus")
pscore_id <- core_members(ps_comp,detection=0.1/100, prevalence=50/100) #more than 0.1%relative abundance in >50% sample
ps_core.taxa <- core(ps_comp, detection = 0.1/100, prevalence = 50/100, sort=TRUE)
pscore90_id <- core_members(ps_comp, detection = 0.1/100, prevalence = 90/100) #higher prevalence
ps_core.taxa
write.csv(core11,"figures/core11.csv")
all <- data.frame(tax_table(taxaname_psprev_unq),Prevalence = prevalence(taxaname_psprev_unq,count=TRUE), # Add taxonomy and total read counts to this data.frame
                   TotalAbundance = taxa_sums(taxaname_psprev_unq))
taxonomy <- as.data.frame(tax_table(ps_core.taxa)) #get taxa id
core_taxa_id <- subset(taxonomy,rownames(taxonomy) %in% pscore_id)
hcore_taxa_id <- subset(taxonomy,rownames(taxonomy) %in% pscore90_id)
#check
table_coreabund <- sample_sums(core(ps_comp, detection = 0.1/100, prevalence = 50/100)) 
table <- as.data.frame(prevalence(ps_comp, detection=0.1/100, sort=TRUE, count=TRUE))
#heatmap
prevalences <- seq(.05, 0.5, .1)
detections <- round(10^seq(log10(1e-3), log10(.2), length = 9),3)
library(RColorBrewer)
#labelled using genera
taxa_names(ps_comp) <- paste0("ASV",seq(ntaxa(ps_comp)),":g_",tax_core$Genus) #assign prefix and numbers
tax_core <- as.data.frame(tax_table(ps_comp))
tax_core$Genus
#ps_comp_fig <- microbiomeutilities::format_to_besthit(ps_comp,prefix = '') #Naming to the genus/species
plot_core(ps_comp, plot.type = "heatmap", 
          colours = rev(brewer.pal(5, "RdGy")), prevalences = prevalences, detections=detections, 
          #min.prevalence=prevalence(ps_comp_figure, sort=TRUE)[77])+
          min.prevalence = 0.5)+
  xlab("Detection Threshold (Relative Abundance (%))")+
  theme_bw(base_size=20)+theme(legend.key.size = unit(0.8,'cm'),legend.text = element_text(size=14))+
  guides(fill=guide_legend(label.vjust = 0.5))



####################################
##### DIVERSITY INDEX ALPHA ###
library(ggpubr)
library(rstatix)
{tab <- microbiome::alpha(psprev_unq, index ="all")
head(tab)
ps_meta <- meta(psprev_unq) #combined meta files and metrics
ps_meta$Shannon <- tab$diversity_shannon
ps_meta$Observed <- tab$observed
ps_meta$Pielou <- tab$evenness_pielou
}
#not sure whether to test for residuals or data
#test for normality ####choose "diversity_shannon","observed","evenness_pielou"
metrics <- c("evenness_pielou") 
shapiro.test(tab[,({{metrics}})]) #p more than 0.05-normal
plot(density(tab[,({{metrics}})]))
ggqqplot(tab[,({{metrics}})])
hist(tab[,({{metrics}})])
#COLOR =JEHAI"#339900", TEMIAR "#FF6600",TEMUAN"#CC33CC"

#####################
{#####kruskal wallis adjusted- GROUP
  ps_meta <- ps_meta %>%  reorder_levels(group, order = c("Jehai_FC",  "Temiar_RC","Temuan_UC")) #reorder categ
  #Observed #interpretation of eff #0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
  observed.kruskal <- ps_meta %>% kruskal_test(Observed ~ group)
  ps_meta %>% kruskal_effsize(Observed ~ group) #eta squared based on H stats, test effect size
  #pairwise using wilcoxon's test
  pwc.k_observed <- ps_meta %>% wilcox_test(Observed ~ group, p.adjust.method = "BH")
  #pwc.k_observed <- ps_meta %>% dunn_test(Observed ~ group, p.adjust.method = "bonferroni")
  pwc.k_observed <- pwc.k_observed %>% add_xy_position(x = "group")
  Observed <- ggboxplot(ps_meta, "group", "Observed", color="black", alpha=0.3,
                        fill= "group",palette=c("#339900", "#FF6600","#CC33CC","black"),outlier.shape = NA) +
    geom_jitter(aes(color = group),width=0.15,size=2) +
    stat_pvalue_manual(pwc.k_observed) +
    theme(axis.title.x = element_blank(),panel.background = element_rect(color = "gray", size = 0.5),
          plot.caption = element_text(size=10))+
    theme_bw(base_size = 24)+labs(x="",  y="Value") +
    labs(caption = get_test_label(observed.kruskal, detailed = FALSE),
         subtitle = "Observed Richness")
  Observed
  #Pielou
  #The interpretation of eff #0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
  Pielou.kruskal <- ps_meta %>% kruskal_test(Pielou ~ group)
  ps_meta %>% kruskal_effsize(Pielou ~ group) #eta squared based on H stats, test effect size
  #pairwise using wilcoxon's test
  pwc.k_Pielou <- ps_meta %>% wilcox_test(Pielou ~ group, p.adjust.method = "BH")
  ps_meta %>% dunn_test(Pielou ~ group, p.adjust.method = "BH")
  pwc.k_Pielou <- pwc.k_Pielou %>% add_xy_position(x = "group")
  Pielou <- ggboxplot(ps_meta, "group", "Pielou", color="black", alpha=0.3,
                      fill= "group",palette=c("#339900", "#FF6600","#CC33CC","black"),outlier.shape = NA) +
    geom_jitter(aes(color = group),width=0.15,size=2) +
    stat_pvalue_manual(pwc.k_Pielou) +
    theme(axis.title.x = element_blank(),panel.background = element_rect(color = "gray", size = 0.5),
          plot.caption = element_text(size=10))+
    theme_bw(base_size = 24)+labs(x="",  y="") +
    labs(caption = get_test_label(Pielou.kruskal, detailed = FALSE),
         subtitle = "Pielou")
  Pielou
  #Shannon
  Shannon.kruskal <- ps_meta %>% kruskal_test(Shannon ~ group)
  ps_meta %>% kruskal_effsize(Shannon ~ group) #eta squared based on H stats, test effect size
  #pairwise using wilcoxon's test
  pwc.k_Shannon <- ps_meta %>% wilcox_test(Shannon ~ group, p.adjust.method = "BH")
  pwc.k_Shannon <- pwc.k_Shannon %>% add_xy_position(x = "group")
  Shannon <- ggboxplot(ps_meta, "group", "Shannon", color="black", alpha=0.3,
                       fill= "group",palette=c("#339900", "#FF6600","#CC33CC","black"),outlier.shape = NA) +
    geom_jitter(aes(color = group),width=0.15,size=2) +
    stat_pvalue_manual(pwc.k_Shannon) +
    theme(axis.title.x = element_blank(),panel.background = element_rect(color = "gray", size = 0.5),
          plot.caption = element_text(size=10))+
    theme_bw(base_size = 24)+labs(x="",  y="") +
    labs(caption = get_test_label(Shannon.kruskal, detailed = FALSE),
         subtitle = "Shannon")
  Shannon
  alpha_group <- ggarrange(Observed,Pielou, Shannon,ncol =3, nrow = 1,
                           #labels = c("A.", "B.","C."),  
                           common.legend = TRUE, legend="top")
  alpha_group
}
ggsave("figures/alpha-s_groupOA.jpeg", width=50, height=25, units = "cm", bg = "white")

{#####kruskal wallis adjusted- waistCM
  ps_metawaist <- filter(ps_meta, !waistCM == "")
  ps_metawaist <- ps_metawaist %>%  reorder_levels(waistCM, order = c("Healthy","Risk")) #reorder categ
  #Observed
  observed.kruskal_waist <- ps_metawaist %>% kruskal_test(Observed ~ waistCM)
  ps_metawaist %>% kruskal_effsize(Observed ~ waistCM) #eta squared based on H stats, test effect size
  #pairwise using wilcoxon's test
  pwc.k_observed_waist <- ps_metawaist %>% wilcox_test(Observed ~ waistCM, p.adjust.method = "BH")
  pwc.k_observed_waist <- pwc.k_observed_waist %>% add_xy_position(x = "waistCM")
  Observed <- ggboxplot(ps_metawaist, "waistCM", "Observed", color="black", alpha=0.4,
                        fill= "waistCM",palette=c("deepskyblue1", "cadetblue"),outlier.shape = NA) +
    geom_jitter(aes(color = waistCM),width=0.15,size=2) +
    stat_compare_means(method="wilcox.test",comparisons=list(c("Healthy","Risk")),label = "p.signif",size=5)+
    theme(axis.title.x = element_blank(),panel.background = element_rect(color = "gray", size = 0.5),
          plot.caption = element_text(size=6))+
    theme_bw(base_size = 24)+labs(x="",  y="Value") +
    labs(caption = get_test_label(observed.kruskal_waist, detailed = FALSE),
         subtitle = "Observed Rich.")
  Observed
  #Pielou
  pielou.kruskal_waist <- ps_metawaist %>% kruskal_test(Pielou ~ waistCM)
  ps_metawaist %>% kruskal_effsize(Pielou ~ waistCM) #eta squared based on H stats, test effect size
  #pairwise using wilcoxon's test
  pwc.k_pielou_waist <- ps_metawaist %>% wilcox_test(Pielou ~ waistCM, p.adjust.method = "BH")
  pwc.k_pielou_waist <- pwc.k_pielou_waist %>% add_xy_position(x = "waistCM")
  Pielou <- ggboxplot(ps_metawaist, "waistCM", "Pielou", color="black", alpha=0.4,
                      fill= "waistCM",palette=c("deepskyblue1", "cadetblue"),outlier.shape = NA) +
    geom_jitter(aes(color = waistCM),width=0.15,size=2) +
    stat_compare_means(method="wilcox.test",comparisons=list(c("Healthy","Risk")),label = "p.signif",size=5)+
    theme(axis.title.x = element_blank(),panel.background = element_rect(color = "gray", size = 0.5),
          plot.caption = element_text(size=6))+
    theme_bw(base_size = 24)+labs(x="",  y="") +
    labs(caption = get_test_label(pielou.kruskal_waist, detailed = FALSE),
         subtitle = "Pielou")
  Pielou
  #Shannon
  Shannon.kruskal_waist <- ps_metawaist %>% kruskal_test(Shannon ~ waistCM)
  ps_metawaist %>% kruskal_effsize(Shannon ~ waistCM) #eta squared based on H stats, test effect size
  #pairwise using wilcoxon's test
  pwc.k_shannon_waist <- ps_metawaist %>% wilcox_test(Shannon ~ waistCM, p.adjust.method = "BH")
  pwc.k_shannon_waist <- pwc.k_shannon_waist %>% add_xy_position(x = "waistCM")
  Shannon <- ggboxplot(ps_metawaist, "waistCM", "Shannon", color="black", alpha=0.4,
                       fill= "waistCM",palette=c("deepskyblue1", "cadetblue"),outlier.shape = NA) +
    geom_jitter(aes(color = waistCM),width=0.15,size=2) +
    stat_compare_means(method="wilcox.test",comparisons=list(c("Risk","Healthy")),label = "p.signif",size=5)+
    theme(axis.title.x = element_blank(),panel.background = element_rect(color = "gray", size = 0.5),
          plot.caption = element_text(size=6))+
    theme_bw(base_size = 24)+labs(x="",  y="") +
    labs(caption = get_test_label(Shannon.kruskal_waist, detailed = FALSE),
         subtitle = "Shannon")
  Shannon
  alpha_waist <- ggarrange(Observed,Pielou, Shannon,ncol =3, nrow = 1,
                           #labels = c("A.", "B.","C."),  
                           common.legend = TRUE, legend="top")
}
alpha_waist
ggarrange(alpha_group, NULL, alpha_waist,ncol =3, nrow = 1,labels = c("A.","", "B."),  widths = c(1.2,0.1,1),
          common.legend = TRUE, legend="top")
ggsave("figures/alpha-s_waist.jpeg", width=50, height=25, units = "cm", bg = "white")


{#####kruskal wallis adjusted- GENDER
ps_meta <- ps_meta %>%  reorder_levels(gender, order = c("F", "M")) #reorder categ
#Observed
observed.kruskal_gender <- ps_meta %>% kruskal_test(Observed ~ gender)
ps_meta %>% kruskal_effsize(Observed ~ gender) #eta squared based on H stats, test effect size
#pairwise using wilcoxon's test
pwc.k_observed_gender <- ps_meta %>% wilcox_test(Observed ~ gender, p.adjust.method = "BH")
pwc.k_observed_gender <- pwc.k_observed_gender %>% add_xy_position(x = "gender")
Observed <- ggboxplot(ps_meta, "gender", "Observed", color="black", alpha=0.4,
                      fill= "gender",palette=c("deeppink", "deepskyblue4"),outlier.shape = NA) +
  geom_jitter(aes(color = gender),width=0.15,size=2) +
  #stat_compare_means(method="wilcox.test",comparisons=list(c("F","M")),label = "p.signif",size=5)+
  theme(axis.title.x = element_blank(),panel.background = element_rect(color = "gray", size = 0.5),
        plot.caption = element_text(size=6))+
  theme_bw(base_size = 20)+labs(x="",  y="Value") +
  labs(caption = get_test_label(observed.kruskal_gender, detailed = FALSE),
       subtitle = "Observed Rich.")
#Pielou
pielou.kruskal_gender <- ps_meta %>% kruskal_test(Pielou ~ gender)
ps_meta %>% kruskal_effsize(Pielou ~ gender) #eta squared based on H stats, test effect size
#pairwise using wilcoxon's test
pwc.k_pielou_gender <- ps_meta %>% wilcox_test(Pielou ~ gender, p.adjust.method = "BH")
pwc.k_pielou_gender <- pwc.k_pielou_gender %>% add_xy_position(x = "gender")
Pielou <- ggboxplot(ps_meta, "gender", "Pielou", color="black", alpha=0.4,
                      fill= "gender",palette=c("deeppink", "deepskyblue4"),outlier.shape = NA) +
  geom_jitter(aes(color = gender),width=0.15,size=2) +
  #stat_compare_means(method="wilcox.test",comparisons=list(c("F","M")),label = "p.signif",size=5)+
  theme(axis.title.x = element_blank(),panel.background = element_rect(color = "gray", size = 0.5),
        plot.caption = element_text(size=6))+
  theme_bw(base_size = 20)+labs(x="",  y="") +
  labs(caption = get_test_label(pielou.kruskal_gender, detailed = FALSE),
       subtitle = "Pielou")
Pielou
#Shannon
Shannon.kruskal_gender <- ps_meta %>% kruskal_test(Shannon ~ gender)
ps_meta %>% kruskal_effsize(Shannon ~ gender) #eta squared based on H stats, test effect size
#pairwise using wilcoxon's test
pwc.k_shannon_gender <- ps_meta %>% wilcox_test(Shannon ~ gender, p.adjust.method = "BH")
pwc.k_shannon_gender <- pwc.k_shannon_gender %>% add_xy_position(x = "gender")
Shannon <- ggboxplot(ps_meta, "gender", "Shannon", color="black", alpha=0.4,
                      fill= "gender",palette=c("deeppink", "deepskyblue4"),outlier.shape = NA) +
  geom_jitter(aes(color = gender),width=0.15,size=2) +
  #stat_compare_means(method="wilcox.test",comparisons=list(c("F","M")),label = "p.signif",size=5)+
  theme(axis.title.x = element_blank(),panel.background = element_rect(color = "gray", size = 0.5),
        plot.caption = element_text(size=6))+
  theme_bw(base_size = 20)+labs(x="",  y="") +
  labs(caption = get_test_label(Shannon.kruskal_gender, detailed = FALSE),
       subtitle = "Shannon")
alpha_gender <- ggarrange(Observed,Pielou, Shannon,ncol =3, nrow = 1,
                         #labels = c("A.", "B.","C."),  
                         common.legend = TRUE, legend="top")
alpha_gender

#####kruskal wallis adjusted- SMOKING STATUS
ps_meta <- ps_meta %>%  reorder_levels(group, order = c("Never","Former","Smoker")) #reorder categ
ps_metasmoke <- filter(ps_meta, !smoke == "N/A")
#Observed
observed.kruskal_smoke <- ps_metasmoke %>% kruskal_test(Observed ~ smoke)
ps_metasmoke %>% kruskal_effsize(Observed ~ smoke) #eta squared based on H stats, test effect size
#pairwise using wilcoxon's test
pwc.k_observed_smoke <- ps_metasmoke %>% wilcox_test(Observed ~ smoke, p.adjust.method = "BH")
pwc.k_observed_smoke <- pwc.k_observed_smoke %>% add_xy_position(x = "smoke")
Observed <- ggboxplot(ps_metasmoke, "smoke", "Observed", color="black", alpha=0.3,
                      fill= "smoke",palette=c("darkgoldenrod1", "brown", "brown1"),outlier.shape = NA) +
  geom_jitter(aes(color = smoke),width=0.15,size=2) +
  #stat_pvalue_manual(pwc.k_observed_smoke) +
  theme(axis.title.x = element_blank(),panel.background = element_rect(color = "gray", size = 0.5),
        plot.caption = element_text(size=6))+
  theme_bw(base_size = 20)+labs(x="",  y="Value") +
  labs(caption = get_test_label(observed.kruskal_smoke, detailed = FALSE),
       subtitle = "Observed Rich.")
#Pielou
Pielou.kruskal_smoke <- ps_metasmoke %>% kruskal_test(Pielou ~ smoke)
ps_metasmoke %>% kruskal_effsize(Pielou ~ smoke) #eta squared based on H stats, test effect size
#pairwise using wilcoxon's test
pwc.k_Pielou_smoke <- ps_metasmoke %>% wilcox_test(Pielou ~ smoke, p.adjust.method = "BH")
pwc.k_Pielou_smoke <- pwc.k_Pielou_smoke %>% add_xy_position(x = "smoke")
Pielou <- ggboxplot(ps_metasmoke, "smoke", "Pielou", color="black", alpha=0.3,
                    fill= "smoke",palette=c("darkgoldenrod1", "brown", "brown1"),outlier.shape = NA) +
  geom_jitter(aes(color = smoke),width=0.15,size=2) +
  #stat_pvalue_manual(pwc.k_Pielou_smoke) +
  theme(axis.title.x = element_blank(),panel.background = element_rect(color = "gray", size = 0.5),
        plot.caption = element_text(size=6))+
  theme_bw(base_size = 20)+labs(x="",  y="") +
  labs(caption = get_test_label(Pielou.kruskal_smoke, detailed = FALSE),
       subtitle = "Pielou")
Pielou
#Shannon
#The interpretation of eff #0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).
Shannon.kruskal_smoke <- ps_metasmoke %>% kruskal_test(Shannon ~ smoke)
ps_metasmoke %>% kruskal_effsize(Shannon ~ smoke) #eta squared based on H stats, test effect size
#pairwise using wilcoxon's test
pwc.k_Shannon_smoke <- ps_metasmoke %>% wilcox_test(Shannon ~ smoke, p.adjust.method = "bonferroni")
pwc.k_Shannon_smoke <- pwc.k_Shannon_smoke %>% add_xy_position(x = "smoke")
Shannon <- ggboxplot(ps_metasmoke, "smoke", "Shannon", color="black", alpha=0.3,
                     fill= "smoke",palette=c("darkgoldenrod1", "brown", "brown1"),outlier.shape = NA) +
  geom_jitter(aes(color = smoke),width=0.15,size=2) +
  #stat_pvalue_manual(pwc.k_Shannon_smoke) +
  theme(axis.title.x = element_blank(),panel.background = element_rect(color = "gray", size = 0.5),
        plot.caption = element_text(size=6))+
  theme_bw(base_size = 20)+labs(x="",  y="") +
  labs(caption = get_test_label(Shannon.kruskal_smoke, detailed = FALSE),
       subtitle = "Shannon")
Shannon
alpha_smoke <- ggarrange(Observed,Pielou, Shannon,ncol =3, nrow = 1,
                         #labels = c("A.", "B.","C."),  
                         common.legend = TRUE, legend="top")
alpha_smoke
ggarrange(alpha_gender, NULL, alpha_smoke,ncol =3, nrow = 1,labels = c("A.","", "B."),  widths = c(1,0.1,1.2),
          common.legend = TRUE, legend="top")
}
ggsave("figures/alphadiv_ns.jpeg.jpeg", width=50, height=25, units = "cm",bg="white")



#############BETA DIVERSITY
microbiome::summarize_phyloseq(psprev_unq)
#PCOA (Traditional)
{##### DIVERSITY INDEX BETA ###
psprev_unq_comp = microbiome::transform(psprev_unq, 'compositional')
wunifrac=ordinate(psprev_unq_comp,"PCoA", "unifrac",weighted=TRUE)
wuni = plot_ordination(psprev_unq_comp,wunifrac,color="group")+geom_point(size=6)+
  scale_color_manual(values=c("#339900", "#FF6600","#CC33CC","black"))+
  #stat_ellipse(aes(group=group), type="norm" ,level =.95, linetype = 1)+ #multivariate t distribution and multi normal
  stat_ellipse(aes(group=group), level =.95, linetype = 2)+ #multivariate normal less conservative
    theme_bw(base_size = 18)+ggtitle("Weighted Unifrac")+scale_y_reverse()+scale_x_reverse()
wuni
unifrac=ordinate(psprev_unq,"PCoA", "unifrac",weighted=FALSE)
uni = plot_ordination(psprev_unq,unifrac,color="group")+geom_point(size=6)+
  scale_color_manual(values=c("#339900", "#FF6600","#CC33CC","black"))+
  stat_ellipse(aes(group=group), level =.95, linetype = 2)+
  theme_bw(base_size = 18)+ggtitle("Unweighted Unifrac")
ggarrange(wuni,uni,common.legend = TRUE,legend = "top")
ggsave("figures/betadiv-unifrac.jpeg", width=45, height=30, units = "cm",bg="white")
###bray
bray=ordinate(psprev_unq,"PCoA", "bray")
bray=plot_ordination(psprev_unq,bray,color="group")+geom_point(size=6)+
  scale_color_manual(values=c("#339900", "#FF6600","#CC33CC","black"))+
  stat_ellipse(aes(group=group), level =.95, linetype = 2)+scale_x_reverse()+
  theme_bw(base_size = 18)+ggtitle("Bray-Curtis")
#PCOA after clr transformation
library(vegan)
#psprev_unq_glom <- tax_glom(psprev_unq, taxrank = 'Family') # agglomerate taxa
ps_clr = microbiome::transform(psprev_unq, 'clr')
gp_euc_pcoa = ordinate(ps_clr, "PCoA", "euclidean")
PCoA_clr.euc = plot_ordination(ps_clr, gp_euc_pcoa, "sample", color="waistCM",shape="waistCM")+
  scale_color_manual(values=c("#339900", "#FF6600","#CC33CC","black"))+
  scale_shape_manual(values=c(19,19,19,17,19,17))+ scale_y_reverse()+
  ggtitle("PCoA clr_euclidean")+theme_bw()+ theme(legend.position = "right")+theme_bw(base_size = 18)+
  #stat_ellipse(type = "norm", linetype = 2)+ 
  stat_ellipse(aes(group=waistCM), level =.95, linetype = 2)+
  geom_point(size=4.5)
PCoA_clr.euc
ggarrange(bray,PCoA_clr.euc,common.legend = TRUE,legend = "top")
}
ggsave("figures/betadiv.jpeg", width=45, height=30, units = "cm",bg="white")




##########################PERMANOVA: test location effect (~mean), assume equal dispersions and linear response to covariate
ps_clr = microbiome::transform(psprev_unq, 'clr')
#metadf <- data.frame(sample_data(ps_clr)) #metadata
clr_dist_matrix <- phyloseq::distance(ps_clr,method = "euclidean")#Generate distance matrix

adonis2(clr_dist_matrix ~ group+gender+smoke+waistCM, data=meta(ps_clr), permutations = 999, method = "euclidean")
adonis2(clr_dist_matrix ~ group, data=meta(ps_clr), permutations = 5000, method = "euclidean")
#adonis2(clr_dist_matrix ~ group+gender+smoke+waistCM, data=meta(ps_clr), permutations = 999)
#adonis2(clr_dist_matrix ~ group+gender+smoke+waistCM, data = metadf, permutations = 999, method = "euclidean")
#different results using margin vs term

#controlled for other factor # control factor (strata)
adonis2(clr_dist_matrix ~ group, strata = meta(ps_clr)$waistCM, data = meta(ps_clr), permutations = 999)
adonis2(clr_dist_matrix ~ group*gender*smoke*waistCM, strata=meta(ps_clr)$waistCM ,data = meta(ps_clr), permutations = 999)
anosim(clr_dist_matrix,metadf$group, permutations = 1000)#alternative to adonis
####################pairwise adonis for group
library(pairwiseAdonis)
padonis = pairwise.adonis(clr_dist_matrix,factors=meta(ps_clr)$group,sim.function='vegdist',
                          sim.method='euclidian',p.adjust.m='BH', perm= 999)
padonis

#permutest(dispr,pairwise = TRUE) #ALTERNATIVE:reject null hypothesis: means difference in dispersion
#Tukey: apply exactly to balanced designs where there are the same number of observations made at each level of the factor.

###betadisper
#group or gender
dispr <- betadisper(clr_dist_matrix, meta(ps_clr)$group, bias.adjust = TRUE)
plot(dispr,main="Multivariate Dispersion: OA Groups",ylab = "PCoA2",xlab = "PCoA1", 
     col=c("#339900", "#FF6600","#CC33CC"), label = FALSE,pch = c(0,1,2),cex=1)
points(dispr$centroids[,1:2], pch=c(15,16,17), cex=1.5,col=c("#339900", "#FF6600","#CC33CC"))
legend("topright", c("Jehai_FC", "Temiar_RC","Temuan_UC"), pch = c(0,1,2), col = c("#339900", "#FF6600","#CC33CC"))
#boxplot(dispr)
#betadisper pairwise CAN use EITHER anova or permutest
anova(dispr) #perform test
thsd <- TukeyHSD(dispr)
plot(thsd)

permutest(dispr,pairwise = TRUE, permutations = 999) #reject null hypothesis: means difference in dispersion
#Pr(>F) overall significance level

#gender
dispr <- betadisper(clr_dist_matrix, meta(ps_clr)$gender, bias.adjust = TRUE)
plot(dispr,main="Multivariate Dispersion: Gender",ylab = "PCoA2",xlab = "PCoA1", 
     col=c("deeppink", "deepskyblue4"), label = FALSE,pch = c(1,2),cex=1)
points(dispr$centroids[,1:2], pch=c(16,17), cex=1.5,col=c("deeppink", "deepskyblue4"))
legend("topright", c("Female", "Male"), pch = c(1,2), col = c("deeppink", "deepskyblue4"))


#WAISTCM
pc_clr_waist <- subset_samples(psprev_unq, !waistCM == "")
ps_clr_waist = microbiome::transform(pc_clr_waist, 'clr')
wclr_dist_matrix <- phyloseq::distance(ps_clr_waist,method = "euclidean")#Generate distance matrix
#wmetadf <- data.frame(sample_data(ps_clr_waist)) #metadata
adonis2(wclr_dist_matrix ~ waistCM, data=meta(ps_clr_waist), permutations = 999, method = "euclidean")
dispr <- betadisper(wclr_dist_matrix, meta(ps_clr_waist)$waistCM, bias.adjust = TRUE)
plot(dispr,main="Multivariate Dispersion: WaistCM",ylab = "PCoA2",xlab = "PCoA1", 
     col=c("deepskyblue1", "cadetblue"), label = FALSE,pch = c(1,2),cex=1)
points(dispr$centroids[,1:2], pch=c(16,17), cex=1.5,col=c("deepskyblue1", "cadetblue"))
legend("topright", c("Healthy", "At risk"), pch = c(1,2), col = c("deepskyblue1", "cadetblue"))
padonis = pairwise.adonis(wclr_dist_matrix,factors=meta(ps_clr_waist)$waistCM,sim.function='vegdist',
                          sim.method='euclidian',p.adjust.m='BH', perm= 9999)
padonis
#smoke
pc_clr_smoke <- subset_samples(psprev_unq, !smoke == "N/A")
pc_clr_smoke = microbiome::transform(pc_clr_smoke, 'clr')
sclr_dist_matrix <- phyloseq::distance(pc_clr_smoke,method = "euclidean")#Generate distance matrix
adonis2(sclr_dist_matrix ~ smoke, data=meta(pc_clr_smoke), permutations = 999, method = "euclidean")
dispr <- betadisper(sclr_dist_matrix, meta(pc_clr_smoke)$smoke, bias.adjust = TRUE)
plot(dispr,main="Multivariate Dispersion: Smoke",ylab = "PCoA2",xlab = "PCoA1", 
     col=c("darkgoldenrod1", "brown", "brown1"), label = FALSE,pch = c(0,1,2),cex=1)
points(dispr$centroids[,1:2], pch=c(15,16,17), cex=1.5,col=c("darkgoldenrod1", "brown", "brown1"))
legend("topright", c("Former", "Smoking","Never"), pch = c(0,1,2), col = c("darkgoldenrod1", "brown", "brown1"))
     
     

adonis2(clr_dist_matrix ~ smoke, data=meta(ps_clr), permutations = 999, method = "euclidean")
padonis = pairwise.adonis(clr_dist_matrix,factors=metadf$waistCM,sim.function='vegdist',
                          sim.method='euclidian',p.adjust.m='BH', perm= 9999)

update.packages(checkBuilt = TRUE, ask = FALSE)

#### ALDEX ####
#Comparison between 3 groups
library("ALDEx2")
#psprev_unq_aldex <- merge_taxa(psprev_unq,"Genus")
#psprev_unq_aldex <- aggregate_taxa(psprev_unq,"Genus")
psprev_unq_aldex<- tax_glom(psprev_unq,"Genus")
otu_table(psprev_unq_aldex)
aldex_taxatable <- as.data.frame(tax_table(psprev_unq_aldex))
taxa_names(psprev_unq_aldex) <- paste0("Genus", seq(ntaxa(psprev_unq_aldex)))
psprev_unq_aldex <- microbiomeutilities::format_to_besthit(psprev_unq_aldex)
taxa_names(psprev_unq_aldex)

TRK_JH = subset_samples(psprev_unq_aldex, group  %in% c("Temiar_RC","Jehai_FC"))
TM_JH = subset_samples(psprev_unq_aldex, group  %in% c("Temuan_UC","Jehai_FC"))
TRK_TM = subset_samples(psprev_unq_aldex, group  %in% c( "Temiar_RC","Temuan_UC"))

#Group 1
#A=sample_data(TRK_JH)
#group1 =A$group
#group1
cond_TRK_JH <- c(rep("Temiar_RC", 24),rep("Jehai_FC", 34))
TRK_JH_df = otu_table(TRK_JH)
TRK_JH.all <- aldex((TRK_JH_df), cond_TRK_JH, mc.samples = 1000, test="t", effect=TRUE, include.sample.summary=TRUE, 
                    denom="all", verbose=TRUE)
#TRK_JH.all <- aldex((TRK_JH_df), group1, mc.samples = 1000, test="t", effect=TRUE, include.sample.summary=TRUE, 
#                    denom="all", verbose=TRUE)
par(mfrow=c(1,2))
aldex.plot(TRK_JH.all, type="MA", test="welch", xlab="Log-ratio abundance", ylab="Difference",cutoff.pval=0.05,
           all.cex=2,called.cex=2,rare.cex=2,all.pch = 20,called.pch = 20,rare.pch = 20)
aldex.plot(TRK_JH.all, type="MW", test="welch", xlab="Dispersion", ylab="Difference",cutoff.pval=0.05,
           all.cex=2,called.cex=2,rare.cex=2,all.pch = 20,called.pch = 20,rare.pch = 20)
tab.TRK_JH.all <- TRK_JH.all %>% filter(TRK_JH.all$we.eBH <= 0.05 & TRK_JH.all$wi.eBH <= 0.05) %>% 
  dplyr::select(effect,wi.eBH,we.eBH)%>%tibble::rownames_to_column("Name")
tab.TRK_JH.all$Enriched <- ifelse(tab.TRK_JH.all$effect > 0, "Temiar_RC","Jehai_FC")
a <- ggplot(data=tab.TRK_JH.all, aes(x=reorder(Name, effect), y=effect, fill=Enriched)) +
  geom_bar(stat="identity")+coord_flip()+
  theme(axis.title.y = element_blank(), panel.background = element_rect(fill = "darkgray"), panel.grid.major = element_line(color = "lightgrey", size = 0.5),axis.title.x =  element_blank(),
        panel.grid.minor.x=element_line(color = "lightgrey",size=0.5), axis.line = element_line(colour = "lightgrey"),
        text = element_text(size = 18))+scale_fill_manual(values = c("#339900","#FF6600"))+
        geom_hline(yintercept=c(1,-1), linetype="dotted", color = "black",size=1.5)+  ylim(-1.8, 2.8)
a
#tabeff.JH_TRK.all <- tab.JH_TRK.all %>% filter(tab.JH_TRK.all$we.eBH <= 0.05 & tab.JH_TRK.all$wi.eBH <= 0.05 &abs(effect) >=1) %>% 
#  dplyr::select(effect,wi.eBH,we.eBH)%>%tibble::rownames_to_column("Name")

#Group 2
cond_TM_JH <- c(rep("Temuan_UC", 29),rep("Jehai_FC", 34))
TM_JH_df = otu_table(TM_JH)
TM_JH.all <- aldex((TM_JH_df), cond_TM_JH, mc.samples = 1000, test="t", effect=TRUE, include.sample.summary=TRUE, 
                    denom="all", verbose=TRUE)
aldex.plot(TM_JH.all, type="MA", test="welch", xlab="Log-ratio abundance", ylab="Difference",cutoff.pval=0.05,
           all.cex=2,called.cex=2,rare.cex=2,all.pch = 20,called.pch = 20,rare.pch = 20)
aldex.plot(TM_JH.all, type="MW", test="welch", xlab="Dispersion", ylab="Difference",cutoff.pval=0.05,
           all.cex=2,called.cex=2,rare.cex=2,all.pch = 20,called.pch = 20,rare.pch = 20)
tab.TM_JH.all <- TM_JH.all %>% filter(TM_JH.all$we.eBH <= 0.05 & TM_JH.all$wi.eBH <= 0.05) %>% 
  dplyr::select(effect,wi.eBH,we.eBH)%>%tibble::rownames_to_column("Name")
tab.TM_JH.all$Enriched <- ifelse(tab.TM_JH.all$effect > 0, "Temuan_UC","Jehai_FC")
b <- ggplot(data=tab.TM_JH.all, aes(x=reorder(Name, effect), y=effect, fill=Enriched)) +
  geom_bar(stat="identity")+coord_flip()+
  theme(axis.title.y = element_blank(), panel.background = element_rect(fill = "darkgray"), panel.grid.major = element_line(color = "lightgrey", size = 0.5),axis.title.x =  element_blank(),
        panel.grid.minor.x=element_line(color = "lightgrey",size=0.5), axis.line = element_line(colour = "lightgrey"),
        text = element_text(size = 18))+scale_fill_manual(values = c("#339900","#CC33CC"))+
  geom_hline(yintercept=c(1,-1), linetype="dotted", color = "black",size=1.5)+  ylim(-1.8, 2.8)
b
#group3
cond_TRK_TM <- c(rep("Temiar_RC", 24),rep("Temuan_UC", 29))
TRK_TM_df = otu_table(TRK_TM)
TRK_TM.all <- aldex((TRK_TM_df), cond_TRK_TM, mc.samples = 1000, test="t", effect=TRUE, include.sample.summary=TRUE, 
                   denom="all", verbose=TRUE)
aldex.plot(TRK_TM.all, type="MA", test="welch", xlab="Log-ratio abundance", ylab="Difference",cutoff.pval=0.05,
           all.cex=2,called.cex=2,rare.cex=2,all.pch = 20,called.pch = 20,rare.pch = 20)
aldex.plot(TRK_TM.all, type="MW", test="welch", xlab="Dispersion", ylab="Difference", cutoff.pval=0.05,
           all.cex=2,called.cex=2,rare.cex=2,all.pch = 20,called.pch = 20,rare.pch = 20)
tab.TRK_TM.all <- TRK_TM.all %>% filter(TRK_TM.all$we.eBH <= 0.05 & TRK_TM.all$wi.eBH <= 0.05) %>% 
  dplyr::select(effect,wi.eBH,we.eBH)%>%tibble::rownames_to_column("Name")
tabeff.TRK_TM.all <- TRK_TM.all %>% filter(TRK_TM.all$we.eBH <= 0.05 & TRK_TM.all$wi.eBH <= 0.05 & abs(TRK_TM.all$effect) >=1) %>% 
  dplyr::select(effect,wi.eBH,we.eBH)%>%tibble::rownames_to_column("Name")
tab.TRK_TM.all$Enriched <- ifelse(tab.TRK_TM.all$effect > 0, "Temiar_RC","Temuan_UC")
c <- ggplot(data=tab.TRK_TM.all, aes(x=reorder(Name, effect), y=effect, fill=Enriched)) +
  geom_bar(stat="identity")+coord_flip()+
  theme(axis.title.y = element_blank(), panel.background = element_rect(fill = "darkgray"), panel.grid.major = element_line(color = "lightgrey", size = 0.5),
        panel.grid.minor.x=element_line(color = "lightgrey",size=0.5), axis.line = element_line(colour = "lightgrey"),
        text = element_text(size = 18))+scale_fill_manual(values = c("#FF6600","#CC33CC"))+
  geom_hline(yintercept=c(1,-1), linetype="dotted", color = "black",size=1.5)+  ylim(-1.8, 2.8)
c
ggarrange(NULL,a,NULL,b,NULL,c, nrow = 6, align ="v",legend = "right",  heights =  c(0.08,1.1,0.08,.6,0.08,0.4),label.y = 1.15,label.x=-0.05,
          labels = c("Jehai_FC vs Temiar_RC:","","","Jehai_FC vs Temuan_UC:", "","Temuan_UC vs Temiar_RC:"))

#merge table
tab.TRK_TM.all$compare <- "TRK_TM"
tab.TM_JH.all$compare <- "TM_JH"
tab.TRK_JH.all$compare <- "TRK_JH"
combine <- rbind(tab.TRK_TM.all,tab.TM_JH.all,tab.TRK_JH.all)
write.csv(quote=FALSE, combine,"figures/aldex_table",row.names = FALSE)



#### ALDEX ####
#Comparison between 3 groups
library("ALDEx2")
psprev_unq_aldex<- tax_glom(psprev_unq,"Genus")
tax_table(psprev_unq_aldex)
taxa_names(psprev_unq_aldex) <- paste0("Genus", seq(ntaxa(psprev_unq_aldex)))
psprev_unq_aldex <- microbiomeutilities::format_to_besthit(psprev_unq_aldex)
taxa_names(psprev_unq_aldex)

waist_H = subset_samples(psprev_unq_aldex, waistCM  %in% c("Healthy"))
waist_R = subset_samples(psprev_unq_aldex, waistCM  %in% c("Risk"))
waist_HR <- cbind(otu_table(waist_H),otu_table(waist_R))
cond_waist <- c(rep("Healthy", 39),rep("Risk", 46))
sample_data(waist_HR)
TM_JH = subset_samples(psprev_unq_aldex, group  %in% c("Temuan_UC","Jehai_FC"))
TRK_TM = subset_samples(psprev_unq_aldex, group  %in% c( "Temiar_RC","Temuan_UC"))
colnames(waist_HR)
#group1 =A$group
#group1
cond_TRK_JH <- c(rep("Temiar_RC", 24),rep("Jehai_FC", 34))
TRK_JH_df = otu_table(TRK_JH)
waist_HR.all <- aldex((waist_HR), cond_waist, mc.samples = 1000, test="t", effect=TRUE, include.sample.summary=TRUE, 
                    denom="all", verbose=TRUE)
aldex.plot(waist_HR.all, type="MA", test="welch", xlab="Log-ratio abundance", ylab="Difference",cutoff.pval=0.05,
           all.cex=2,called.cex=2,rare.cex=2,all.pch = 20,called.pch = 20,rare.pch = 20)
aldex.plot(waist_HR.all, type="MW", test="welch", xlab="Dispersion", ylab="Difference",cutoff.pval=0.05,
           all.cex=2,called.cex=2,rare.cex=2,all.pch = 20,called.pch = 20,rare.pch = 20)
tab.waist_HR.all <- waist_HR.all %>% filter(waist_HR.all$we.eBH <= 0.05 & TRK_JH.all$wi.eBH <= 0.05) %>% 
  dplyr::select(effect,wi.eBH,we.eBH)%>%tibble::rownames_to_column("Name")
tab.TRK_JH.all$group <- ifelse(tab.TRK_JH.all$effect > 0, "Temiar_RC","Jehai_FC")


#gender
gender_M = subset_samples(psprev_unq_aldex, gender  %in% c("M"))
gender_F = subset_samples(psprev_unq_aldex, gender  %in% c("F"))
gender_MF <- cbind(otu_table(gender_M),otu_table(gender_F))
cond_gender <- c(rep("Male", 25),rep("Female", 62))
sample_data(gender_MF)
colnames(gender_MF)
#group1 =A$group
#group1
gender_MF.all <- aldex((gender_MF), cond_gender, mc.samples = 1000, test="t", effect=TRUE, include.sample.summary=TRUE, 
                    denom="all", verbose=TRUE)
aldex.plot(gender_MF.all, type="MA", test="welch", xlab="Log-ratio abundance", ylab="Difference",cutoff.pval=0.05,
           all.cex=2,called.cex=2,rare.cex=2,all.pch = 20,called.pch = 20,rare.pch = 20)
aldex.plot(gender_MF.all, type="MW", test="welch", xlab="Dispersion", ylab="Difference",cutoff.pval=0.05,
           all.cex=2,called.cex=2,rare.cex=2,all.pch = 20,called.pch = 20,rare.pch = 20)
tab.gender_MF.all <- gender_MF.all %>% filter(TRK_JH.all$we.eBH <= 0.05 & TRK_JH.all$wi.eBH <= 0.05) %>% 
  dplyr::select(effect,wi.eBH,we.eBH)%>%tibble::rownames_to_column("Name")




