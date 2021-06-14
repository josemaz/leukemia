require(TCGAbiolinks)
require(SummarizedExperiment)
require(NOISeq)
require(EDASeq)

query.rna <- GDCquery(project = "TARGET-ALL-P1",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - Counts")
GDCdownload(query.rna)
rnas <- GDCprepare(query = query.rna, summarizedExperiment = TRUE)

# View(as.data.frame(colData(rnas))) # See clinical info

rnas.PB <- rnas[,rnas$sample_type == "Primary Blood Derived Cancer - Peripheral Blood"]
rnas.BM <- rnas[,rnas$sample_type == "Primary Blood Derived Cancer - Bone Marrow"]

# Filtro
m1 <- assay(rnas)
dim(m1)
threshold <- round(dim(rnas)[2]/2)
print(paste0("Threshold: ",threshold))
m1<- m1[rowSums(m1 == 0) <= threshold, ]
print(paste0("Rows After Zeros (dim): ",dim(m1)[1]))
m1 <- m1[rowMeans(m1) >= 10, ]
print(paste0("Rows After means (dim): ",dim(m1)[1]))
rnas <- rnas[rownames(rnas) %in% rownames(m1),]

# TAREA 1: Descargar biomart y anotar
rowData(rnas) # tiene que tener las columnas de biomart (gene_length y gc content, chrom, etc.)
# Tarea 2: Correr Normalizaciones

# Normalization
fac <- data.frame(sampleType= rnas$sample_type, 
                  row.names=colnames(rnas))
ln.data <- withinLaneNormalization(assay(rnas), 
                                   rowData(rnas)$geneLength, which = "full")