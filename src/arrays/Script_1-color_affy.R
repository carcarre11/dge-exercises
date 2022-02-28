#Master in Bioinformatics 2022
#Normalization and diferential expression for Affymetrix arrays
# Autor: Gonzalo Gómez. ggomez@cnio.es

#setwd("mypath") ### Establecer directorio de trabajo

##1. Load libraries (do not forget to install them before)
library("affy")
library("limma")
library("genefilter")


#2. Import targets.txt file
targets <- readTargets("targets.txt", row.names="FileName")


#3. Import .CEL files
data <- ReadAffy(filenames=targets$FileName) 		#Import intensities from Affymetrix arrays (.CEL)
													# data is an object of AffyBatch class

#4. Normalize with RMA 
#generates object eset (class ExprSet), 
#expresso function provides intensities in log scale
eset <- expresso(data,
 				bg.correct = TRUE, 
                bgcorrect.method="rma",
                normalize = TRUE, 
                normalize.method="quantiles", 
                pmcorrect.method="pmonly", 
                summary.method="medianpolish",
                verbose = TRUE,
				 ) 




#5. Generate BOXPLOTS before and after normalization

#boxplot for raw data
boxplot(data,
		 main="Boxplot Before Normalization",
		 col = "lightgrey")
		

#boxplot for normalized data
exprseset <- as.data.frame(exprs(eset))		
boxplot(data.frame(exprseset),
		 main="Boxplot After Normalization (log scale)",
		 col = "white")


#6.Data filtering using IQR.
esetIQR <- varFilter(eset, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)




#######Differential expression analysis.#######
#7. Design matrix.
design<-cbind(Control=c(1,1,1,1,1,0,0,0,0,0), Nanog_RNAi=c(0,0,0,0,0,1,1,1,1,1))
rownames(design)<-targets$FileName

#8. Contrasts matrix.
cont.matrix<-makeContrasts(Nanog_RNAivsControl=Nanog_RNAi-Control,levels=design) 


#9. Obtaining differentially expressed genes (DEGs)
#Linear model and eBayes 
fit<-lmFit(esetIQR,design)  ##getting DEGs from IQR 
fit2<-contrasts.fit(fit, cont.matrix)
fit2<-eBayes(fit2)

#Table with DEGs results
toptableIQR<-topTable(fit2, number=dim(exprs(esetIQR))[1], adjust.method="BH", sort.by="p")


##10. Save results
save(toptableIQR,file="MyResults.RData")






