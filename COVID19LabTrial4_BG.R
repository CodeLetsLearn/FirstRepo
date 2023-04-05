# Introduction ----
# The goal of this short script to it make sure you are able to read the count data and the study design for the COVID hackdash into your R environment

# We only need a single package
library(tidyverse)
# Begin by reading in study design that includes all info for both human and ferret samples
targets.all <- read_tsv("covid_metadata.txt")
# then read in the human covid data and convert to a matrix with gene symbols as rownames
#human_covid_data <- read_tsv("GSE147507_RawReadCounts_Human.tsv")
#human_covid_data <- as.matrix(column_to_rownames(human_covid_data, "...1"))
# repeat for the ferret covid data
ferret_covid_data <- read_tsv("GSE147507_RawReadCounts_Ferret.tsv")
ferret_covid_data <- as.matrix(column_to_rownames(ferret_covid_data, "...1"))

# Now proceed with your exploration and analysis of the data!

#Ferret trechea difference in presentation of SARS COV, and other viruses, in ferret trachea at 3 days

# Question -----
# What is the difference in presentation of IAV and 
# SARS-CoV-2 in ferret trachea at 3 days?

#Part 1 of Jackie's ferret code-----
# start processing the targets object for ferret data
targets <- targets.all %>%
  filter(host_species == "Mustela_putorius_furo",
         tissue_cell_type == "trachea")

targets$sample <- c("Ctl_d3_1",      
                    "Ctl_d3_2",      
                    "Ctl_d3_3",       
                    "Ctl_d3_4",       
                    "SARSCoV2_d3_1",
                    "SARSCoV2_d3_2",
                    "SARSCoV2_d3_3",
                    "SARSCoV2_d3_4",
                    "IAV_d3_1",       
                    "IAV_d3_2",       
                    "IAV_d3_3",       
                    "IAV_d3_4",       
                    "IAV_d3_5",       
                    "IAV_d3_6")

targets$group <- c("Mock",
                   "Mock",
                   "Mock",
                   "Mock",
                   "SARS",
                   "SARS",
                   "SARS",
                   "SARS",
                   "IAV",
                   "IAV",
                   "IAV",
                   "IAV",
                   "IAV",
                   "IAV")

Txi_gene <- ferret_covid_data[,19:32]

colnames(Txi_gene)[11] <- "SARSCoV2_d3_1"
colnames(Txi_gene)[12] <- "SARSCoV2_d3_2"
colnames(Txi_gene)[13] <- "SARSCoV2_d3_3"
colnames(Txi_gene)[14] <- "SARSCoV2_d3_4"


#Part 2 of Jackie's ferret code/optimized for all Data wrangling (Step2 script) ----
library(tidyverse)
library(edgeR)
library(matrixStats)
library(cowplot)

sampleLabels <- c("Ctl_d3_1", 
                  "Ctl_d3_2",
                  "Ctl_d3_3",
                  "Ctl_d3_4",
                  "IAV_d3_1",      
                  "IAV_d3_2",      
                  "IAV_d3_3",       
                  "IAV_d3_4", 
                  "IAV_d3_5",
                  "IAV_d3_6",
                  "SARSCoV2_d3_1",
                  "SARSCoV2_d3_2",
                  "SARSCoV2_d3_3",
                  "SARSCoV2_d3_4")

myDGEList <- DGEList(Txi_gene)

log2.cpm <- cpm(myDGEList, log=TRUE)

#Question: Not sure if I should call these "Transcript ID" or "gene ID", I think I'm working on gene ID but the original script from Jackie says "transcript ID"
log2.cpm.df <- as_tibble(log2.cpm, rownames = "transcriptID") 
colnames(log2.cpm.df) <- c("transcriptID", sampleLabels)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = -1, # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)

p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 16, 
               size = 2, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()
p1 <- p1 + coord_flip()

cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1)>=4 #user defined
myDGEList.filtered <- myDGEList[keepers,]

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "transcriptID")
colnames(log2.cpm.filtered.df) <- c("transcriptID", sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                           cols = -1, # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)
p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 16, 
               size = 2, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()
p2 <- p2 + coord_flip()

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("transcriptID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = -1, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)
p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 16, 
               size = 2, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()
p3 <- p3 + coord_flip()

plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)

#Part 3 of Jackie's ferret code/optimized for SARS vs. Flu comparison (Step3 script) ----
library(tidyverse)
library(DT)
library(gt)
library(plotly)

group <-targets$group 
group <- factor(group)

#PCA of PC1 vs. PC2
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
pc.var <- pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1) 
pca.res.df <- as_tibble(pca.res$x)
pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color = group) +
  geom_point(size=2) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  # coord_fixed() +
  theme_bw()

ggplotly(pca.plot)

#PCA of PC2 vs. PC3

pca.plot <- ggplot(pca.res.df) +
  aes(x=PC2, y=PC3, label=sampleLabels, color = group) +
  geom_point(size=2) +
  stat_ellipse() +
  xlab(paste0("PC2 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC3 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  # coord_fixed() +
  theme_bw()

ggplotly(pca.plot)


#PCA of PC1 vs. PC3

pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC3, label=sampleLabels, color = group) +
  geom_point(size=2) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC3 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  # coord_fixed() +
  theme_bw()

ggplotly(pca.plot)

#Log fold change of transcript counts between ferrets infected with SARS-CoV-2 and IAV treatment. 
mydata.df <- mutate(log2.cpm.filtered.norm.df,
                    #Mock.AVG = (Ctl_d3_1 +
                    #              Ctl_d3_2 +
                    #              Ctl_d3_3 +
                    #               Ctl_d3_4)/4, 
                    SARS.CoV.2.AVG = (SARSCoV2_d3_1 +
                                        SARSCoV2_d3_2 +
                                        SARSCoV2_d3_3 +
                                        SARSCoV2_d3_4)/4,
                     IAV.AVG = (IAV_d3_1 +       
                               IAV_d3_2 +      
                               IAV_d3_3 +      
                               IAV_d3_4 +       
                               IAV_d3_5 +       
                               IAV_d3_6)/6,
                    #now make columns comparing each of the averages above that you're interested in
                    LogFC = (SARS.CoV.2.AVG - IAV.AVG),
                    # LogFC.IAV = (IAV.AVG - mock.AVG)
) %>% #note that this is the first time you've seen the 'pipe' operator
  mutate_if(is.numeric, round, 2)

datatable(mydata.df[,c(1,16:18)], #note: the Avg for my COVID data, and IAV (flu), and logFoldChange data are in column 16-18 respectively, that's why I specified those columns here
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         #dom = "Blfrtip", 
                         #buttons = c("copy", "csv", "excel"),
                         lengthMenu = c("10", "25", "50", "100")))


#Part 4 of Jackie's ferret code/optimized for SARS vs. Flu comparison (DiffGenes (Step 5 script)) ----
library(tidyverse)
library(limma)
library(edgeR)
library(gt)
library(DT)
library(plotly)

group <- factor(group)
design <- model.matrix(~0 + group)

colnames(design) <- levels(group)

v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)
fit <- lmFit(v.DEGList.filtered.norm, design)
contrast.matrix <- makeContrasts(infectionSARS = SARS - IAV,
                                 levels=design)

fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "transcriptID")

#Plot volcano plot comparison of SARS-COV2 vs. IAV
vplot <- ggplot(myTopHits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", transcriptID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(title="Ferret trachea infections with SARS-CoV-2 versus IAV",
       subtitle = "Volcano Plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

ggplotly(vplot)

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=1)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
diffGenes.df <- as_tibble(diffGenes, rownames = "transcriptID")
datatable(diffGenes.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: Differential transcript counts in Ferret trachea SARS-CoV-2 infection',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:9), digits=2)

#Note: there were no differentially expressed genes between SARS vs. IAV!
#Note: your adjusted p-values are determined on your sample number and your number of genes, 
#So if you're interested only a subset of genes, you can just focus on those genes (ie: focus only on immune genes)


#Part 5 of Jackie's ferret code/optimized for SARS vs. Flu comparison (Modules (Step 6 script)) ----
library(gplots)
library(RColorBrewer)
myheatcolors <- rev(brewer.pal(name="RdBu", n=11))
clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") #cluster rows by pearson correlation
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
heatmap.2(diffGenes, 
          Rowv=as.dendrogram(clustRows), 
          # Colv=as.dendrogram(clustColumns),
          Colv = NA,
          RowSideColors=module.color,
          col=myheatcolors, scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20))

modulePick <- 2 
myModule_up <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
hrsub_up <- hclust(as.dist(1-cor(t(myModule_up), method="pearson")), method="complete") 

heatmap.2(myModule_up, 
          Rowv=as.dendrogram(hrsub_up), 
          Colv=NA, 
          labRow = NA,
          col=myheatcolors, scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.color[module.assign%in%modulePick], margins=c(8,20))

modulePick <- 1 
myModule_down <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
hrsub_down <- hclust(as.dist(1-cor(t(myModule_down), method="pearson")), method="complete") 

heatmap.2(myModule_down, 
          Rowv=as.dendrogram(hrsub_down), 
          Colv=NA, 
          labRow = NA,
          col=myheatcolors, scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.color[module.assign%in%modulePick], margins=c(8,20))

#Notes from Lab discussion-----
#Deconvolution: Check out this website to get immune deconvolution:
#https://omnideconv.org/immunedeconv/articles/immunedeconv.html
#You can adjust "Scale by row" for the the heatmap; (each values are compared to average for each row)
#overrepresented vs. underrepresented instead of upregulated vs. downregulated-Ask 
