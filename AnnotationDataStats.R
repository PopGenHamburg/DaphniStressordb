Sys.setenv(LANG = "en")

##########################################
### Literature Database statistics
##########################################

#1) Data Availabilty Statistics:

# How many papers available for each search strategy: 

search.strategy <- data.frame(c("FunctionSearch","FunctionSearch","FunctionSearch","Mendeley","Mendeley","Mendeley","Twitter","Twitter","Twitter"))
category <- data.frame(c("Obtained", "Daphnia_GE", "Retrievable_info","Obtained", "Daphnia_GE", "Retrievable_info","Obtained", "Daphnia_GE", "Retrievable_info"))
numbers.obtained <- data.frame(c(47,24,21,74,45,43,37,31,26))
data.avail <- data.frame(cbind(search.strategy, category, numbers.obtained))
colnames(data.avail) <- c("SearchStrategy", "Category", "numbers")
data.avail

#install.packages("ggplot2")
library(ggplot2)

data.avail$Category <- factor(data.avail$Category, levels=c("Obtained", "Daphnia_GE", "Retrievable_info"))
## here we change the order in which we want the bar plot to appear, by default, R alphabetically orders the category levels and hence we can change it to desired order using the above command

ggplot(data.avail, aes(SearchStrategy, numbers)) + geom_bar(aes(fill=Category), width=0.4, position=position_dodge(width=0.5), stat="identity") +  
  theme(legend.position="top", legend.title = element_blank(),axis.title.x=element_blank(), axis.title.y=element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size=35)) + scale_fill_grey()

ggsave("Literature_stats.eps", plot=last_plot(), device="eps", path="D:/PhD/Work/Project2_Annotation/AnnotationDataStats/OUTput/new_plots/final_lit_stats_plots/")

########################################
### Literature information
########################################

lit.data.input <- read.table("D:/PhD/Work/Project2_Annotation/AnnotationDataStats/INput/LitData_input.tab", header=T, sep="\t")
head(lit.data.input)


#------------------ Number of papers for each species

lit.data.yes <- subset(lit.data.input, Data_Availability =='Yes')

lit.species.yes <- lit.data.yes$species
head(lit.species.yes)

species.yes.tab <- table(lit.species.yes)
species.yes.tab

#install.packages("plotrix")
library('plotrix')

setEPS()
postscript("D:/PhD/Work/Project2_Annotation/AnnotationDataStats/OUTput/new_plots/final_lit_stats_plots/Species_info.eps", width=50, height=30)
par(oma=c(10,10,10,20))
pie(species.yes.tab, labels=c(expression(italic("")),expression(italic("D.carinata")(4)),expression(italic("Daphnia sp.")(2)), expression(italic("D.galeata")(1)),expression(italic("D.magna")(52)),expression(italic("")),expression(italic("D.pulex")(30)), expression(italic("D.pulex+ D.pulicaria")(1))), 
    col = c("gray", "gray0", "lavenderblush4", "lightblue4", "lightyellow4", "palegreen3", "sandybrown", "wheat3", "pink"),
    main="Species distribution", cex=3, cex.main=6)
dev.off()

#------------------ Timeline of papers

lit.data.yes <- subset(lit.data.input, Data_Availability =='Yes')

lit.year.all <- lit.data.input$Year
lit.year.yes <- lit.data.yes$Year
lit.year.tab <- table(lit.year.yes)

setEPS()
postscript("D:/PhD/Work/Project2_Annotation/AnnotationDataStats/OUTput/new_plots/final_lit_stats_plots/Year_info.eps", width=50, height=50)
par(oma=c(20,10,0,0))
barplot(lit.year.tab, ylim=c(0,30), col="grey", cex=8, cex.axis=7, cex.names=5, width=1, axisnames = FALSE)
dev.off()


#------------------ Methods used

lit.data.yes <- subset(lit.data.input, Data_Availability =='Yes')

lit.methods.yes <- lit.data.yes$methods

lit.methods.yes.tab <- table(lit.methods.yes)
lit.methods.yes.tab

lit.methods.val <- c(42,30,7,7,1,1,1)
names(lit.methods.val) <- c("qPCR", "Microarray", "qPCR\n+ Microarray", "RNA-Seq", "Assay", "Comparative\n Genomics", "Protein 2D-Gel")

setEPS()
postscript("D:/PhD/Work/Project2_Annotation/AnnotationDataStats/OUTput/new_plots/final_lit_stats_plots/methods_info.eps", width=30, height=30)
par(oma=c(20,30,10,10))
barplot(lit.methods.val, space=1, xlim = c(0,50), cex.lab=3, cex.axis=3, cex.names=3, las=1, horiz = T)
mtext(text="Methods used", side=2, outer=TRUE, cex=4, line=15)
mtext(text="Number of studies", side=1, outer=TRUE, cex=4, line=2)
dev.off()

##-------------------------- Venn diagram for BLAST results

##########################################
### Finding overlaps between BLASTx and tBLASTx
##########################################

blastx.out <- read.table("D:/PhD/Work/Project2_Annotation/AnnotationDataStats/BLAST_new/blastx_more_50_filt.tab", header=T, sep="\t")
tblastx.out <- read.table("D:/PhD/Work/Project2_Annotation/AnnotationDataStats/BLAST_new/tblastx_more_50_filt.tab", header=T, sep="\t")
head(blastx.out)
head(tblastx.out)

dgal.blastx <- blastx.out$QuerySequence
dgal.tblastx <- tblastx.out$QuerySequence

dgal.blast.overlap <- list(dgal.blastx, dgal.tblastx)
sort.dgal.blast.overlap <- sort(unique(unlist(dgal.blast.overlap)))

dgal.blast.mat <- matrix(rep(0, length(sort.dgal.blast.overlap)*length(dgal.blast.overlap)), ncol=2)
colnames(dgal.blast.mat) <- c("BLASTx", "tBLASTx")
rownames(dgal.blast.mat) <- sort.dgal.blast.overlap
lapply(seq_along(dgal.blast.overlap), function(i){
  dgal.blast.mat[sort.dgal.blast.overlap %in% dgal.blast.overlap[[i]], i] <<- table(dgal.blast.overlap[[i]])
})

dgal.blast.mat

#source("http://www.bioconductor.org/biocLite.R")
#biocLite("limma")
library(limma)
install.packages("Cairo")

setEPS()
cairo_ps("D:/PhD/Work/Project2_Annotation/AnnotationDataStats/BLAST_new/BLAST_summary.eps")
vennDiagram(dgal.blast.mat, circle.col = c("blue3", "deeppink"))
dev.off()

#vennDiagram(mat.blast, names=c("tblastx", "blastx"), main="Number of transcripts overlapping", cex=c(1.5,1.5), lwd=2, circle.col = c("blue", "pink"), counts.col = c("black"))
#dev.off()


#------------ pie charts

## dappu: there are 6452 out of 30939 genes with stressor as on 2nd Feb 2018 => 20.85%
### dapma: there are 7096 out of 29127 genes with stressor as on 2nd Feb 2018 => 24.36%

## data input
#total=100
#dpul.stress=19.61 # Percentage of genes that contain an associated stressor from literature
#dpul.no.stress=total-dpul.stress  #  Percentage of genes that DONT contain an associated stressor from literature
dpul.slices <- c(20.85, 79.15)
dpul.labels <- c("Stressor", "NoStressor")
setEPS()
postscript("D:/PhD/Work/Project2_Annotation/AnnotationDataStats/OUTput/new_plots/final_lit_stats_plots/Dpul_stress_pie.eps")
pie(dpul.slices, labels=dpul.labels, cex=2, cex.main=1, col = c("black", "grey"),  main=substitute(paste("Number of ", italic(D.pulex), " genes with an associated stressor in literature")))
dev.off()

dmag.slices <- c(24.36, 75.64)
dmag.labels <- c("Stressor", "NoStressor")
setEPS()
postscript("D:/PhD/Work/Project2_Annotation/AnnotationDataStats/OUTput/new_plots/final_lit_stats_plots/Dmag_stress_pie.eps")
pie(dmag.slices, labels=dmag.labels, cex=2, cex.main=1, col = c("black", "grey"),  main=substitute(paste("Number of ", italic(D.magna), " genes with an associated stressor in literature")))
dev.off()


dgal.slices <- c(14.23, 85.76)
dgal.labels <- c("Stressor", "NoStressor")
setEPS()
postscript("D:/PhD/Work/Project2_Annotation/AnnotationDataStats/OUTput/new_plots/final_lit_stats_plots/Dgal_stress_pie.eps")
pie(dgal.slices, labels=dgal.labels, cex=2, cex.main=4, col = c("black", "grey"),  main=substitute(paste("Number of ", italic(D.galeata), " transcripts: Comparative approach")))
dev.off()

