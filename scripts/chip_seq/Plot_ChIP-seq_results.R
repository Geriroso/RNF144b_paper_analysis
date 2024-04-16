### Plot ChIP-seq results

### Libraries
suppressPackageStartupMessages(library(Gviz))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(Mus.musculus))
suppressPackageStartupMessages(library(Homo.sapiens))
suppressPackageStartupMessages(library(here))

### Create TxDb object for mm10
#### Check the table info from mm10
#supportedUCSCtables(genome="mm10")
#### Obtain the desired gene identifier
#mm10_txdb <- makeTxDbFromUCSC(genome="mm10", tablename="refGene")
#### Saving the object
#saveDb(mm10_txdb, file=here("data","mm10_NCBI_RefSeq_Genes.sqlite"))

#### Obtain the desired gene identifier
#hg38_txdb <- makeTxDbFromUCSC(genome="hg38", tablename="refGene")
#saveDb(hg38_txdb, file=here("data","hg38_NCBI_RefSeq_Genes.sqlite"))

### Load data
hg38_txdb <- loadDb(here("data","hg38_NCBI_RefSeq_Genes.sqlite"))
mm10_txdb <- loadDb(here("data","mm10_NCBI_RefSeq_Genes.sqlite"))
P53_WT_BCELL_MOCK_bw = here("results", "Nextflow","bwa","mergedLibrary","bigwig","P53_WT_BCELL_MOCK.bigWig")
P53_WT_BCELL_IR_bw = here("results", "Nextflow","bwa","mergedLibrary","bigwig","P53_WT_BCELL_IR.bigWig")
P53_NULL_IR_bw = here("results", "Nextflow","bwa","mergedLibrary","bigwig","P53_NULL_SPLEEN_IR.bigWig")
P53_WT_NONBCELL_MOCK_bw = here("results", "Nextflow","bwa","mergedLibrary","bigwig","P53_WT_NONBCELL_MOCK.bigWig")
P53_WT_NONBCELL_IR_bw = here("results", "Nextflow","bwa","mergedLibrary","bigwig","P53_WT_NONBCELL_IR.bigWig")
alternative = here("results", "Nextflow","bwa","mergedLibrary","bigwig","SRX769393.bw")
## GSE46240
MEF_P53_NULL_DOX_bw = here("results","Nextflow_GSE46240","MEF_P53_NULL_DOX.bigWig")
MEF_P53_WT_DOX_bw = here("results","Nextflow_GSE46240","MEF_P53_WT_DOX.bigWig")
MEF_P53_NULL_DOX_bed = here("results","Nextflow_GSE46240","MEF_P53_NULL_DOX_Rnf144b.bed")
MEF_P53_WT_DOX_bed = here("results","Nextflow_GSE46240","MEF_P53_WT_DOX_Rnf144b.bed")

## GSE55727
MEF_P53_WT_DOX_bw_GSE55727 = here("results","Nextflow_MEF_GSE55727","MEF_1_P53_WT_DOX.bigWig")
MEF_P53_WT_DOX_bed_GSE55727 = here("results","Nextflow_MEF_GSE55727","MEF_1_P53_WT_DOX_Rnf144b.bed")

## GSE55727_FSF
INPUT_FSF_GM00011_WT_DOX_bw_GSE55727 = here("results","Nextflow_FSF_GSE55727","INPUT_FSF_GM00011_WT_DOX.bigWig")
INPUT_FSF_GM06170_WT_DOX_bw_GSE55727 = here("results","Nextflow_FSF_GSE55727","INPUT_FSF_GM06170_WT_DOX.bigWig")
FSF_GM00011_P53_WT_DOX_bw_GSE55727 = here("results","Nextflow_FSF_GSE55727","FSF_GM00011_P53_WT_DOX.bigWig")
FSF_GM06170_P53_WT_DOX_bw_GSE55727 = here("results","Nextflow_FSF_GSE55727","FSF_GM06170_P53_WT_DOX.bigWig")
FSF_GM00011_P53_WT_DOX_bed_GSE55727 = here("results","Nextflow_FSF_GSE55727","FSF_GM00011_P53_WT_DOX_Rnf144b.bed")
FSF_GM06170_P53_WT_DOX_bed_GSE55727 = here("results","Nextflow_FSF_GSE55727","FSF_GM06170_P53_WT_DOX_Rnf144b.bed")


### Called_peaks in bed file after running broadPeak2bed.sh
P53_WT_BCELL_IR_bed_path = here("results", "Nextflow","bwa","mergedLibrary","macs2","broadPeak","P53_WT_BCELL_IR_Rnf144b.bed")
P53_WT_NONBCELL_IR_bed_path = here("results", "Nextflow","bwa","mergedLibrary","macs2","broadPeak","P53_WT_NONBCELL_IR_Rnf144b.bed")
### Convert bed in GRanges object
P53_WT_BCELL_IR_bed_df = read.table(P53_WT_BCELL_IR_bed_path)
P53_WT_NONBCELL_IR_bed_df = read.table(P53_WT_NONBCELL_IR_bed_path)
P53_WT_BCELL_IR_GR = GRanges(seqnames = P53_WT_BCELL_IR_bed_df$V1,
				ranges = IRanges(start = P53_WT_BCELL_IR_bed_df$V2,
												 end = P53_WT_BCELL_IR_bed_df$V3),
				strand = "*")
P53_WT_NONBCELL_IR_GR = GRanges(seqnames = P53_WT_NONBCELL_IR_bed_df$V1,
														 ranges = IRanges(start = P53_WT_NONBCELL_IR_bed_df$V2,
														 								 end = P53_WT_NONBCELL_IR_bed_df$V3),
														 strand = "*")

### Convert GSE55727 human to GRanges
FSF_GM00011_P53_WT_DOX_bed_GSE55727_df = read.table(FSF_GM00011_P53_WT_DOX_bed_GSE55727)
FSF_GM00011_P53_WT_DOX_GR_GSE55727 = GRanges(seqnames = FSF_GM00011_P53_WT_DOX_bed_GSE55727_df$V1,
				ranges = IRanges(start = FSF_GM00011_P53_WT_DOX_bed_GSE55727_df$V2,
												 end = FSF_GM00011_P53_WT_DOX_bed_GSE55727_df$V3),
				strand = "*")
FSF_GM06170_P53_WT_DOX_bed_GSE55727_df = read.table(FSF_GM06170_P53_WT_DOX_bed_GSE55727)
FSF_GM06170_P53_WT_DOX_GR_GSE55727 = GRanges(seqnames = FSF_GM06170_P53_WT_DOX_bed_GSE55727_df$V1,
																						 ranges = IRanges(start = FSF_GM06170_P53_WT_DOX_bed_GSE55727_df$V2,
																						 								 end = FSF_GM06170_P53_WT_DOX_bed_GSE55727_df$V3),
																						 strand = "*")

### Convert GSE55727 mouse to GRanges
MEF_P53_WT_DOX_bed_GSE55727_df = read.table(MEF_P53_WT_DOX_bed_GSE55727)
MEF_P53_WT_DOX_GR_GSE55727 = GRanges(seqnames = MEF_P53_WT_DOX_bed_GSE55727_df$V1,
																		 ranges = IRanges(start = MEF_P53_WT_DOX_bed_GSE55727_df$V2,
																		 								 end = MEF_P53_WT_DOX_bed_GSE55727_df$V3),
																		 strand = "*")

### Convert GSE46240 to GRanges
MEF_P53_WT_DOX_bed_df = read.table(MEF_P53_WT_DOX_bed)
MEF_P53_WT_DOX_GR = GRanges(seqnames = MEF_P53_WT_DOX_bed_df$V1,
														ranges = IRanges(start = MEF_P53_WT_DOX_bed_df$V2,
																						 end = MEF_P53_WT_DOX_bed_df$V3),
														strand = "*")


### Plotting BCELL
gtrack = GenomeAxisTrack()
itrack = IdeogramTrack(genome = "mm10", chromosome = "chr13")
txTr = GeneRegionTrack(mm10_txdb, transcriptAnnotation = "symbol", showId = TRUE, geneSymbol = TRUE, name = "Gene", col.axis = "black",background.title = "white", col.title = "black")
#### Get gene symbol from Transcript Name
z = ranges(txTr)
z$symbol <- mapIds(Mus.musculus, z$gene, "SYMBOL", "GENEID")
ranges(txTr) = z
dtrack_bcell_mock_bw = DataTrack(range = P53_WT_BCELL_MOCK_bw, genome = "mm10", type = "l", name = "BCELL_mock_bW", ylim = c(0,2), col.axis = "black",background.title = "white", col.title = "black") 
dtrack_bcell_ir_bw = DataTrack(range = P53_WT_BCELL_IR_bw, genome = "mm10", type = "l", name = "BCELL_IR_bW", ylim = c(0,2), col.axis = "black",background.title = "white", col.title = "black") 
dtrack_p53_null_bw = DataTrack(range = P53_NULL_IR_bw, genome = "mm10", type = "l", name = "p53_null_bW", ylim = c(0,2), col.axis = "black",background.title = "white", col.title = "black") 
dtrack_alternative = DataTrack(range = alternative, genome = "mm10", type = "l", name = "altern", ylim = c(0,2))

bptrack_bcell_ir = AnnotationTrack(P53_WT_BCELL_IR_GR, name = "Called Peaks", genome = "mm10")
plotTracks(list(itrack,gtrack,dtrack_bcell_mock_bw,dtrack_bcell_ir_bw,dtrack_p53_null_bw,bptrack_bcell_ir,txTr), from = 47080000, to = 47165000, showBandId = TRUE, col = NULL, extend.left = 600, extend.right = 500, ylim = c(0,1.5))

# Mouse_BCELL_ylim1.5_long.png
pdf(file = here("results","BCELL_ChIP_seq.pdf"), width = 15, height = 10)
plotTracks(list(itrack,gtrack,dtrack_bcell_mock_bw,dtrack_bcell_ir_bw,dtrack_p53_null_bw,txTr), from = 47107000, to = 47250000, showBandId = TRUE, col = NULL, extend.left = 600, extend.right = 500, ylim = c(0,1.5))
dev.off()

## Alternative
plotTracks(list(itrack,gtrack,dtrack_alternative,txTr), from = 47120000, to = 47250000, showBandId = TRUE, col = NULL, extend.left = 6000, extend.right = 5000)

### Plotting NONBCELL
dtrack_nonbcell_mock_bw = DataTrack(range = P53_WT_NONBCELL_MOCK_bw, genome = "mm10", type = "l", name = "NONBCELL_mock_bW", ylim = c(0,2), col.axis = "black",background.title = "white", col.title = "black") 
dtrack_nonbcell_ir_bw = DataTrack(range = P53_WT_NONBCELL_IR_bw, genome = "mm10", type = "l", name = "NONBCELL_IR_bW", ylim = c(0,2),col.axis = "black",background.title = "white", col.title = "black") 
dtrack_p53_null_bw = DataTrack(range = P53_NULL_IR_bw, genome = "mm10", type = "l", name = "p53_null_bW", ylim = c(0,2),col.axis = "black",background.title = "white", col.title = "black") 
bptrack_nonbcell_ir = AnnotationTrack(P53_WT_NONBCELL_IR_GR, name = "Called Peaks", genome = "mm10")

pdf(file = here("results","NONBCELL_ChIP_seq.pdf"), width = 15, height = 10)
plotTracks(list(itrack,gtrack,dtrack_nonbcell_mock_bw,dtrack_nonbcell_ir_bw,dtrack_p53_null_bw,txTr), from = 47107000, to = 47250000, showBandId = TRUE, col = NULL, extend.left = 600, extend.right = 500, ylim = c(0,1.5))
dev.off()


### Plotting GSE46240
dtrack_MEF_P53_WT_DOX_bw = DataTrack(range = MEF_P53_WT_DOX_bw, genome = "mm10", type = "l", name = "MEF_P53_WT_DOX_bW", ylim = c(0,2))
dtrack_MEF_P53_NULL_DOX_bw = DataTrack(range = MEF_P53_NULL_DOX_bw, genome = "mm10", type = "l", name = "MEF_P53_NULL_DOX_bW", ylim = c(0,2))
bptrack_MEF_P53_WT_DOX = AnnotationTrack(MEF_P53_WT_DOX_GR, name = "Called Peaks", genome = "mm10")
plotTracks(list(itrack,gtrack,dtrack_MEF_P53_WT_DOX_bw,dtrack_MEF_P53_NULL_DOX_bw,bptrack_MEF_P53_WT_DOX,txTr), from = 47120000, to = 47250000, showBandId = TRUE, col = NULL, extend.left = 6000, extend.right = 5000)

### Plotting GSE55727
dtrack_MEF_P53_WT_DOX_GSE55727_bw = DataTrack(range = MEF_P53_WT_DOX_bw_GSE55727, genome = "mm10", type = "l", name = "MEF_P53_WT_DOX_bW", ylim = c(0,2))
bptrack_MEF_P53_WT_DOX_GSE55727 = AnnotationTrack(MEF_P53_WT_DOX_GR_GSE55727, name = "Called Peaks", genome = "mm10")
plotTracks(list(itrack,gtrack,dtrack_MEF_P53_WT_DOX_GSE55727_bw,bptrack_MEF_P53_WT_DOX_GSE55727,txTr), from = 47120000, to = 47250000, showBandId = TRUE, col = NULL, extend.left = 6000, extend.right = 5000)


### Plotting GSE55727 human
gtrack_hg38 = GenomeAxisTrack()
itrack_hg38 = IdeogramTrack(genome = "hg38", chromosome = "chr6")
txTr_hg38 = GeneRegionTrack(hg38_txdb, transcriptAnnotation = "symbol", showId = TRUE, geneSymbol = TRUE, name = "Gene",col.axis = "black",background.title = "white", col.title = "black")
#### Get gene symbol from Transcript Name
z_hg38 = ranges(txTr_hg38)
z_hg38$symbol <- mapIds(Homo.sapiens, z_hg38$gene, "SYMBOL", "GENEID")
ranges(txTr_hg38) = z_hg38
dtrack_FSF_GM00011_p53_WT_DOX_bw = DataTrack(range = FSF_GM00011_P53_WT_DOX_bw_GSE55727, genome = "hg38", type = "l", name = "FSF_GM00011_WT_DOX", ylim = c(0,2), col.axis = "black",background.title = "white", col.title = "black") 
dtrack_FSF_GM06170_p53_WT_DOX_bw = DataTrack(range = FSF_GM06170_P53_WT_DOX_bw_GSE55727, genome = "hg38", type = "l", name = "FSF_GM06170_WT_DOX", ylim = c(0,2), col.axis = "black",background.title = "white", col.title = "black") 
dtrack_INPUT_FSF_GM00011_p53_WT_DOX_bw = DataTrack(range = INPUT_FSF_GM00011_WT_DOX_bw_GSE55727, genome = "hg38", type = "l", name = "INPUT_FSF_GM00011_WT_DOX", ylim = c(0,2), col.axis = "black",background.title = "white", col.title = "black") 
dtrack_INPUT_FSF_GM06170_p53_WT_DOX_bw = DataTrack(range = INPUT_FSF_GM06170_WT_DOX_bw_GSE55727, genome = "hg38", type = "l", name = "INPUT_FSF_GM06170_WT_DOX", ylim = c(0,2), col.axis = "black",background.title = "white", col.title = "black") 
bptrack_FSF_GM00011_p53_WT_DOX = AnnotationTrack(FSF_GM00011_P53_WT_DOX_GR_GSE55727, name = "Called Peaks", genome = "hg38")
bptrack_FSF_GM06170_p53_WT_DOX = AnnotationTrack(FSF_GM06170_P53_WT_DOX_GR_GSE55727, name = "Called Peaks", genome = "hg38")
p = plotTracks(list(itrack_hg38,gtrack_hg38,dtrack_FSF_GM00011_p53_WT_DOX_bw,dtrack_INPUT_FSF_GM00011_p53_WT_DOX_bw,bptrack_FSF_GM00011_p53_WT_DOX,dtrack_FSF_GM06170_p53_WT_DOX_bw, dtrack_INPUT_FSF_GM06170_p53_WT_DOX_bw,bptrack_FSF_GM06170_p53_WT_DOX,txTr_hg38), from = 18386200, to = 18468870, showBandId = TRUE, col = NULL, extend.left = 6000, extend.right = 5000)

pdf(file = here("results","GSE55727_FSF_P53_WT_DOX.pdf"), width = 15, height = 10)
plotTracks(list(itrack_hg38,gtrack_hg38,dtrack_FSF_GM00011_p53_WT_DOX_bw,dtrack_INPUT_FSF_GM00011_p53_WT_DOX_bw,dtrack_FSF_GM06170_p53_WT_DOX_bw,dtrack_INPUT_FSF_GM06170_p53_WT_DOX_bw,txTr_hg38), from = 18376200, to = 18468870, showBandId = TRUE, col = NULL, extend.left = 600, extend.right = 500, ylim=c(0,2))
dev.off()

pdf(file = here("results","GSE55727_FSF_GM00011_P53_WT_DOX.pdf"), width = 15, height = 10)
plotTracks(list(itrack_hg38,gtrack_hg38,dtrack_FSF_GM00011_p53_WT_DOX_bw,dtrack_INPUT_FSF_GM00011_p53_WT_DOX_bw,txTr_hg38), from = 18376200, to = 18468870, showBandId = TRUE, col = NULL, extend.left = 600, extend.right = 500, ylim = c(0,2))
dev.off()

pdf(file = here("results","GSE55727_FSF_GM06170_P53_WT_DOX.pdf"), width = 15, height = 10)
plotTracks(list(itrack_hg38,gtrack_hg38,dtrack_FSF_GM06170_p53_WT_DOX_bw,dtrack_INPUT_FSF_GM06170_p53_WT_DOX_bw,txTr_hg38), from = 18376200, to = 18468870, showBandId = TRUE, col = NULL, extend.left = 600, extend.right = 500, ylim=c(0,2))
dev.off()

