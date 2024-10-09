#!/bin/Rscript

#install.packages(c('ape', 'dplyr', 'graphics', 'grDevices', 
#   'magrittr', 'phytools', 'shape', 'stats', 'stringr', 'tibble', 
#   'tidyr', 'utils', 'viridis'), dep=T)
#library(remotes)
#install_github("ChrispinChaguza/RCandy", build_vignettes = FALSE)



library(RCandy)
library(ape)
library(tidyr)
library(RcmdrMisc)	# for Hist() and numSummary()
library(abind)
library(e1071)

verbose=FALSE

# set input file names
tree.file <- "zyg_gubbins.node_labelled.final_tree.tre"
tree.file <- "zyg_gubbins.final_tree.tre"
metadata.file <- "zyg_metadata.tsv"
gubbins.gff.file <- "zyg_gubbins.recombination_predictions.gff"
ref.genome.gff.file <- "EcoliK12.gff3"
outdir <- "zyg_R"

outf <- function(f) {
    return( paste(outdir, f, sep="/") )
}

min_max_normalize <- function(x) {
    (x - min(x)) / (max(x) - min(x))
}

printf <- function(...) cat(sprintf(...))

make.backup <- function(file) {
    i = 0
    while (i < 1000) {
        # since we just add a number, any path will also be preserved
        bck <- sprintf("%s,%03d", file, i)
        if (! file.exists(bck) ) break;
        i <- i + 1
    }
    file.copy(file, bck)
}

say.and.continue <- function(condition) {
    cat('\n\n    *** ERROR ***\n\n')
    return()
}


# Save all output
if (verbose) {
    logfile <- outf("R.RCandy.log")
} else {
    logfile <- outf("R.RCandy.sum")
}

make.backup(logfile)
sink(logfile, split=TRUE)

cat('\n\n-------------------------------------\n')
cat(    "Reading data for statistical analysis\n")
cat(    '-------------------------------------\n\n')


cat("Reading Gubbins output\n")
if ( ! file.exists(outf("R.gubbins.GFF.Rds")) ) {
    gubbins.GFF <- load.gubbins.GFF(gubbins.gff.file, 
                   recom.input.type = "Gubbins") 
    saveRDS(gubbins.GFF, outf("R.gubbins.GFF.Rds"))
    # cannot be saved as tab-separated
} else {
    gubbins.GFF <- readRDS(outf("R.gubbins.GFF.Rds"))
}

cat("Reading recombination events\n")
if ( ! file.exists(outf("R.rec.events.Rds")) ) {
    #rec.events <- load.gubbins.GFF(gubbins.gff.file, recom.input.type = "Gubbins")
    rec.events <- gubbins.GFF
    saveRDS(rec.events, file=outf("R.rec.events.Rds"))
    # load with rec.events <- readRDS("R.rec.events.Rds")
} else {
    rec.events <- readRDS(outf("R.rec.events.Rds"))
    # cannot be saved as tab-separated
}

cat("Reading phylogenetic tree\n")
if (!  file.exists(outf("R.tree.Rds")) ) {
    #tree <- ape::read.tree.file(tree.file)
    tree <- read.tree.file(tree.file)
    saveRDS(tree, outf("R.tree.Rds"))
} else {
    tree <- readRDS(outf("R.tree.Rds"))
    # cannot be saved as tab-separated
}

cat("Reading metadata\n")
if ( ! file.exists(outf("R.metadata.df.Rds")) ) {
    metadata.df <- load.taxon.metadata(metadata.file)
    saveRDS(metadata.df, outf("R.metadata.df.Rds"))
    # cannot be saved as tab-separated
} else {
    metadata.df <- readRDS(outf("R.metadata.df.Rds"))
}
metadata <- tidyr::as_tibble(read.table(metadata.file, header=T, sep='\t', comment.char="?"))

cat("Reading reference genome\n")
if ( ! file.exists(outf("R.refgenome.gff.Rds")) ) {
    refgenome.GFF <- load.genome.GFF(ref.genome.gff.file)
    saveRDS(refgenome.GFF, outf("R.refgenome.GFF.Rds"))
    # cannot be saved as tab-separated
} else {
    refgenome.GFF <- readRDS(outf("R.refgenome.GFF.Rds"))
}

cat("Reading rec frequency per base\n")
if ( ! file.exists(outf("R.rec.freq.Rds")) ) {
    rec.freq <- count.rec.events.per.base(gubbins.gff.file, recom.input.type="Gubbins")
    saveRDS(rec.freq, file=outf("R.rec.freq.Rds"))
    write.table(rec.freq, file=outf("R.rec.freq.tab"), sep='\t', row.names=FALSE)
} else {
    rec.freq <- readRDS(outf("R.rec.freq.Rds"))
}

cat("Reading rec events per genome\n")
if ( ! file.exists(outf("R.rec.genome.Rds")) ) {
    rec.genome <- count.rec.events.per.genome(gubbins.gff.file, 
		    recom.input.type="Gubbins",
                    #taxon.names=metadata.df$ID)
                    taxon.names=tree$tip.label)
    saveRDS(rec.genome, file=outf("R.rec.genome.Rds"))
    write.table(rec.genome, file=outf("R.rec.genome.tab"), sep='\t', row.names=FALSE)
    head(rec.genome)
} else {
    rec.genome <- readRDS(outf("R.rec.genome.Rds"))
}

genome.recombinations <- data.frame(rec=rec.genome)
colnames(genome.recombinations) <- c('genome', 'rec.Freq')
head(genome.recombinations)





if(verbose) {
    # save all (screen) plots to a PDF file
    make.backup(outf("R.RCandy_plots.pdf"))
    pdf(outf("R.RCandy_plots.pdf"), paper="a4", width=8, height=11)

    # we may specify specific taxons names to include in the figure (e.g.)
    #tree1 <- ape::read.tree(tree.file)
    #subtree.taxa <- tree1$tip.label[1:50]
    # and add ", subtree.tips=subtree.taxa to the plot command
    # let us make some useful subsets:
    subtree.resistant.taxa <- metadata[ metadata$resistant=="Y", ]$ID
    subtree.sensible.taxa <- metadata[ metadata$resistant=="N", ]$ID
    # we have one repeated ID
    subtree.sensible.taxa <- subtree.sensible.taxa[subtree.sensible.taxa != 'ERR871394']

    # run RCandy
    cat("Creating simple plot with tree and categories\n")
    RCandyVis(tree.file.name = tree, #tree.file, 
	      midpoint.root = TRUE, 
	      ladderize.tree.right = TRUE, 
    	      taxon.metadata.file = metadata, #metadata.file, 
	      taxon.metadata.columns = c("resistant", "S83L", "D87N", "S83L.D87N")
             )


    cat("Creating plots with genomic recombinations\n")
    if ( ! file.exists(outf("R.RCandy_plot.pdf")) ) {
	# do a separate PDF plot annotated with the GFF genome information
	# and the recombination frequency per genome
	RCandyVis(
            # tree
            tree.file.name=tree, #tree.file, 
            midpoint.root=TRUE, 
            ladderize.tree.right=TRUE, 
            #subtree.tips=NULL,
            #color.tree.tips.by.column="resistant",
            #trait.for.ancestral.reconstr='S83L'
            #ace.model.name="ER",	# "ARD", "ER", "SYM"
            #tree.scale.length=NULL,
            #show.tip.label=FALSE,
            #align.tip.label=FALSE,
            # metadata
            taxon.metadata.file=metadata, #metadata.file, 
            taxon.metadata.columns = c("S83L", "D87N", "S83L.D87N", "resistant"), 
            taxon.id.column = "ID", 
            show.metadata.columns=TRUE,
            show.metadata.label=TRUE,
            metadata.column.label.angle=90,
            color.pallette="plasma", # plasma, cividis, viridis, magma, inferno
            # reference genome
            ref.genome.name=refgenome.GFF, #ref.genome.gff.file,
            show.gene.label=FALSE,	# here genes have multiple overwriting labels
            gene.label.angle=90,
            show.genome.ticks=TRUE,
            show.genome.axis=TRUE,
            show.genome.annot=TRUE,
            #ref.genome.length=NULL,
            #genome.start=NULL,
            #genome.end=NULL,
            # recombinations
            gubbins.gff.file=gubbins.GFF, #gubbins.gff.file, 
            recom.input.type="Gubbins",
            show.rec.events=TRUE,
            show.rec.freq.per.genome=TRUE, 
            show.rec.per.genome.scale=TRUE,
            show.rec.freq.per.base=TRUE,
            rec.events.per.base.as.heatmap=FALSE,
            #show.rec.plot.tracks=FALSE, 
            #show.rec.plot.border=FALSE,
            rec.heatmap.color=c("#FF000022", "#0000FF66"),
            #rec.heatmap.color=c("#FFFFFF00", "#0000FF66"),
            rec.plot.bg.transparency=0.0,
            #plot.width=12,
            #plot.height=9.5,
            show.fig.legend=TRUE,
            save.to.this.file=outf("R.RCandy_plot.pdf")
            #save.to.this.file="R.RCandy_plot_term.pdf"
	)

    }
    # Repeat on-screen (or on the global PDF)
    RCandyVis(
	# tree
	tree.file.name=tree, #tree.file, 
	midpoint.root=TRUE, 
	ladderize.tree.right=TRUE, 
	#subtree.tips=NULL,
	#color.tree.tips.by.column="resistant",
	#trait.for.ancestral.reconstr='S83L'
	#ace.model.name="ER",	# "ARD", "ER", "SYM"
	#tree.scale.length=NULL,
	#show.tip.label=FALSE,
	#align.tip.label=FALSE,
	# metadata
	taxon.metadata.file=metadata, #metadata.file, 
	taxon.metadata.columns = c("S83L", "D87N", "S83L.D87N"), 
	taxon.id.column = "ID", 
	show.metadata.columns=TRUE,
	show.metadata.label=TRUE,
	metadata.column.label.angle=90,
	color.pallette="plasma", # plasma, cividis, viridis, magma, inferno
	# reference genome
	ref.genome.name=refgenome.GFF, #ref.genome.gff.file,
	show.gene.label=FALSE,	# here genes have multiple overwriting labels
	gene.label.angle=90,
	show.genome.ticks=TRUE,
	show.genome.axis=TRUE,
	show.genome.annot=TRUE,
	#ref.genome.length=NULL,
	#genome.start=NULL,
	#genome.end=NULL,
	# recombinations
	gubbins.gff.file=gubbins.GFF, #gubbins.gff.file, 
	recom.input.type="Gubbins",
	show.rec.events=TRUE,
	show.rec.freq.per.genome=TRUE, 
	show.rec.per.genome.scale=TRUE,
	show.rec.freq.per.base=TRUE,
	rec.events.per.base.as.heatmap=TRUE,
	#show.rec.plot.tracks=FALSE, 
	#show.rec.plot.border=FALSE,
	rec.heatmap.color=c("#FF000022", "#0000FF66"),
	rec.plot.bg.transparency=0.0,
	#plot.width=12,
	#plot.height=9.5,
	show.fig.legend=TRUE,
	save.to.this.file=NULL
    )

    # And then do subplots for resistant and sensible strains
    RCandyVis(
	# tree
	tree.file.name=tree, #tree.file, 
	midpoint.root=TRUE, 
	ladderize.tree.right=TRUE, 
	subtree.tips=subtree.resistant.taxa,
	# metadata
	taxon.metadata.file=metadata, #metadata.file, 
	taxon.metadata.columns = c("S83L", "D87N", "S83L.D87N"), 
	taxon.id.column = "ID", 
	show.metadata.columns=TRUE,
	show.metadata.label=TRUE,
	metadata.column.label.angle=90,
	color.pallette="plasma", # plasma, cividis, viridis, magma, inferno
	# reference genome
	ref.genome.name=refgenome.GFF, #ref.genome.gff.file,
	show.gene.label=FALSE,	# here genes have multiple overwriting labels
	gene.label.angle=90,
	show.genome.ticks=TRUE,
	show.genome.axis=TRUE,
	show.genome.annot=TRUE,
	# recombinations
	gubbins.gff.file=gubbins.GFF, #gubbins.gff.file, 
	recom.input.type="Gubbins",
	show.rec.events=TRUE,
	show.rec.freq.per.genome=FALSE, 
	show.rec.per.genome.scale=FALSE,
	show.rec.freq.per.base=FALSE,
	rec.events.per.base.as.heatmap=FALSE,
	rec.heatmap.color=c("#FF000022", "#0000FF66"),
	rec.plot.bg.transparency=0.0,
	show.fig.legend=TRUE,
	save.to.this.file=NULL
    )

    RCandyVis(
	# tree
	tree.file.name=tree, #tree.file, 
	midpoint.root=TRUE, 
	ladderize.tree.right=TRUE, 
	subtree.tips=subtree.sensible.taxa,
	# metadata
	taxon.metadata.file=metadata, #metadata.file, 
	taxon.metadata.columns = c("S83L", "D87N", "S83L.D87N"), 
	taxon.id.column = "ID", 
	show.metadata.columns=TRUE,
	show.metadata.label=TRUE,
	metadata.column.label.angle=90,
	color.pallette="plasma", # plasma, cividis, viridis, magma, inferno
	# reference genome
	ref.genome.name=refgenome.GFF, #ref.genome.gff.file,
	show.gene.label=FALSE,	# here genes have multiple overwriting labels
	gene.label.angle=90,
	show.genome.ticks=TRUE,
	show.genome.axis=TRUE,
	show.genome.annot=TRUE,
	# recombinations
	gubbins.gff.file=gubbins.GFF, #gubbins.gff.file, 
	recom.input.type="Gubbins",
	show.rec.events=TRUE,
	show.rec.freq.per.genome=FALSE, 
	show.rec.per.genome.scale=FALSE,
	show.rec.freq.per.base=FALSE,
	rec.events.per.base.as.heatmap=FALSE,
	rec.heatmap.color=c("#FF000022", "#0000FF66"),
	rec.plot.bg.transparency=0.0,
	show.fig.legend=TRUE,
	save.to.this.file=NULL
    )


    cat("Creating plot of recombinations for gyrA region\n")
    #gyrA.start <- 1030000	# in ST-131
    #gyrA.end   <- 1040000
    gyrA.start <- 2336000	# 2336793 in E coli K12
    gyrA.end   <- 2340000	# 2339420 in E coli K12
    if ( ! file.exists("R.RCandy_gyrA_region.pdf") ) {
	# Do a partial plot of the GyrA region (just for the sake of it)
	RCandyVis(tree.file.name = tree, 
            midpoint.root = TRUE, 
            ladderize.tree.right = TRUE, 
            taxon.metadata.file = metadata, 
            taxon.metadata.columns = c("S83L", "D87N", "S83L.D87N"), 
            taxon.id.column = "ID", 
            genome.start = gyrA.start, genome.end = gyrA.end, 
            #genome.start = 1032000, genome.end = 1035000,
            gubbins.gff.file = gubbins.GFF, 
            ref.genome.name = refgenome.GFF, 
            rec.events.per.base.as.heatmap=TRUE,
            show.gene.label = FALSE,
            rec.heatmap.color=c("#FF000022", "#0000FF"),
            rec.plot.bg.transparency=0.0,
            color.pallette="plasma", # plasma, cividis, viridis, magma, inferno
            #color.tree.tips.by.column="S83L",
            #trait.for.ancestral.reconstr='S83L',
            save.to.this.file = outf("R.RCandy_gyrA_region.pdf")
	)
    }
    # And to the screen (or global PDF)
    RCandyVis(tree.file.name = tree, 
	midpoint.root=TRUE, 
	ladderize.tree.right = TRUE, 
	taxon.metadata.file = metadata, 
	taxon.metadata.columns = c("S83L", "D87N", "S83L.D87N"), 
	taxon.id.column = "ID", 
	genome.start = gyrA.start, genome.end = gyrA.end, 
	gubbins.gff.file = gubbins.GFF, 
	ref.genome.name = refgenome.GFF, 
	rec.events.per.base.as.heatmap = TRUE,
	show.gene.label = FALSE,
	rec.heatmap.color=c("#FF000022", "#0000FF"),
	rec.plot.bg.transparency=0.0,
	color.pallette="plasma", # plasma, cividis, viridis, magma, inferno
	#color.tree.tips.by.column="S83L",
	#trait.for.ancestral.reconstr='S83L',
	save.to.this.file = NULL
    )


    cat("Creating plot of recombinations for parC region\n")
    parC.start <- 3163000	# 3163715 in E coli K12
    parC.end   <- 3166000	# 3165973 in E coli K12
    if (! file.exists(outf("R.RCandy_parC_region.pdf")) ) {
	# Do a partial plot of the parC region (just for the sake of it)
	RCandyVis(tree.file.name = tree, 
            midpoint.root = TRUE, 
            ladderize.tree.right = TRUE, 
            taxon.metadata.file = metadata, 
            taxon.metadata.columns = c("S83L", "D87N", "S83L.D87N"), 
            taxon.id.column = "ID", 
            #parC goes from 1956486 to 1958744
            genome.start = 1920000, genome.end = 1980000, 
            gubbins.gff.file = gubbins.GFF, 
            ref.genome.name = refgenome.GFF, 
            rec.events.per.base.as.heatmap=TRUE,
            show.gene.label = FALSE,
            #save.to.this.file = outf("R.RCandy_parC_region.pdf"),
            rec.heatmap.color=c("#FF000022", "#0000FF"),
            rec.plot.bg.transparency=0.0,
            color.pallette="plasma", # plasma, cividis, viridis, magma, inferno
            #color.tree.tips.by.column="S83L",
            #trait.for.ancestral.reconstr='S83L',
            save.to.this.file = outf("R.RCandy_parC_region.pdf")
	)
    }
    # And to the screen (or global PDF)
    RCandyVis(tree.file.name = tree, 
	midpoint.root = TRUE, 
	ladderize.tree.right = TRUE, 
	taxon.metadata.file = metadata, 
	taxon.metadata.columns = c("S83L", "D87N", "S83L.D87N"), 
	taxon.id.column = "ID", 
	#parC goes from 1956486 to 1958744
	genome.start = 1902000, genome.end = 1980000, 
	gubbins.gff.file = gubbins.GFF, 
	ref.genome.name = refgenome.GFF, 
	rec.events.per.base.as.heatmap=TRUE,
	show.gene.label = FALSE,
	rec.heatmap.color=c("#FF000022", "#0000FF"),
	rec.plot.bg.transparency=0.0,
	color.pallette="plasma", # plasma, cividis, viridis, magma, inferno
	#color.tree.tips.by.column="S83L",
	#trait.for.ancestral.reconstr='S83L',
	save.to.this.file = NULL
    )

    # And add the ancestral reconstruction of S83L (which => resistance)
    RCandyVis(tree.file.name = tree, 
	midpoint.root = TRUE, 
	ladderize.tree.right = TRUE, 
	taxon.metadata.file = metadata, 
	taxon.metadata.columns = c("S83L", "D87N", "S83L.D87N"), 
	taxon.id.column = "ID", 
	#parC goes from 1956486 to 1958744
	genome.start = 1902000, genome.end = 1980000, 
	gubbins.gff.file = gubbins.GFF, 
	ref.genome.name = refgenome.GFF, 
	rec.events.per.base.as.heatmap=TRUE,
	show.gene.label = FALSE,
	rec.heatmap.color=c("#FF000022", "#0000FF"),
	rec.plot.bg.transparency=0.0,
	color.pallette="plasma", # plasma, cividis, viridis, magma, inferno
	#color.tree.tips.by.column="S83L",
	#trait.for.ancestral.reconstr='S83L',
	save.to.this.file = NULL
    )


} # endif (verbose)

# ----------------------------
# Prepare statistical analysis
# ----------------------------

cat('\n\n------------------\n')
cat(    "Preparing datasets\n")
cat(    '------------------\n\n')

# BY RESISTANT/SENSIBLE
cat("By R/S status\n")
# select metadata.df$ID that are resistant and not resistant
resistant <- metadata.df$resistant == 'Y'
# get their names so we can select matching genomes
resistant.names <- metadata.df$ID[resistant]

#genome.recombinations$genome %in% resistant.names
# use %in% to get the recombination frequencies from rec.genome
genome.recomb.resistant <- genome.recombinations$genome %in% resistant.names
names.resistant <- as.character(genome.recombinations$genome[genome.recomb.resistant])
freqs.resistant <- genome.recombinations$rec.Freq[genome.recomb.resistant]
names.not.resistant <- as.character(genome.recombinations$genome[!genome.recomb.resistant])
freqs.not.resistant <- genome.recombinations$rec.Freq[!genome.recomb.resistant]
freqs.res <- data.frame(
        genome=c(names.resistant, names.not.resistant),
	freqs=c(freqs.resistant, freqs.not.resistant), 
        ab=c(rep('Y', length(freqs.resistant)), 
             rep('N', length(freqs.not.resistant ))
            )
        )
freqs.res$ab <- as.factor(freqs.res$ab)

# BY MUTATION
cat("By S83L status\n")
s83l <- metadata.df$S83L == 'Y'
s83l.names <- metadata.df$ID[s83l]

genome.recomb.s83l <- genome.recombinations$genome %in% s83l.names
names.s83l <- as.character(genome.recombinations$genome[genome.recomb.s83l])
freqs.s83l <- genome.recombinations$rec.Freq[genome.recomb.s83l]
names.not.s83l <- as.character(genome.recombinations$genome[!genome.recomb.s83l])
freqs.not.s83l <- genome.recombinations$rec.Freq[!genome.recomb.s83l]
freqs.s83l <- data.frame(
        genome=c(names.s83l, names.not.s83l),
	freqs=c(freqs.s83l, freqs.not.s83l), 
        ab=c(rep('Y', length(freqs.s83l)), 
             rep('N', length(freqs.not.s83l ))
            )
        )
freqs.s83l$ab <- as.factor(freqs.s83l$ab)

cat("By D87N status\n")
d87n <- metadata.df$D87N == 'Y'
d87n.names <- metadata.df$ID[d87n]

genome.recomb.d87n <- genome.recombinations$genome %in% d87n.names
names.d87n <- as.character(genome.recombinations$genome[genome.recomb.d87n])
freqs.d87n <- genome.recombinations$rec.Freq[genome.recomb.d87n]
names.not.d87n <- as.character(genome.recombinations$genome[!genome.recomb.d87n])
freqs.not.d87n <- genome.recombinations$rec.Freq[!genome.recomb.d87n]
freqs.d87n <- data.frame(
        genome=c(names.d87n, names.not.d87n),
	freqs=c(freqs.d87n, freqs.not.d87n), 
        ab=c(rep('Y', length(freqs.d87n)), 
             rep('N', length(freqs.not.d87n ))
            )
        )
freqs.d87n$ab <- as.factor(freqs.d87n$ab)

cat("By S83L+D87N status\n")
s83l.d87n <- metadata.df$S83L.D87N == 'Y'
s83l.d87n.names <- metadata.df$ID[s83l.d87n]

genome.recomb.s83l.d87n <- genome.recombinations$genome %in% s83l.d87n.names
names.s83l.d87n <- as.character(genome.recombinations$genome[genome.recomb.s83l.d87n])
freqs.s83l.d87n <- genome.recombinations$rec.Freq[genome.recomb.s83l.d87n]
names.not.s83l.d87n <- as.character(genome.recombinations$genome[!genome.recomb.s83l.d87n])
freqs.not.s83l.d87n <- genome.recombinations$rec.Freq[!genome.recomb.s83l.d87n]
freqs.s83l.d87n <- data.frame(
        genome=c(names.s83l.d87n, names.not.s83l.d87n),
	freqs=c(freqs.s83l.d87n, freqs.not.s83l.d87n), 
        ab=c(rep('Y', length(freqs.s83l.d87n)), 
             rep('N', length(freqs.not.s83l.d87n ))
            )
        )
freqs.s83l.d87n$ab <- as.factor(freqs.s83l.d87n$ab)

# and now prepare all datasets in a list for use in a for loop
freqs.list <- list(
	resistance=freqs.res,
        S83L=freqs.s83l,
        D87N=freqs.d87n,
        S83L.D87N=freqs.s83l.d87n
        )

#
# ANALYZE THE DATASETS
#

for ( f in 1:length(freqs.list) ) {
    freqs <- freqs.list[[f]]
    trait <- names(freqs.list)[f]

    cat("\n\n=======================================\n")
    cat("\n    ANALYZING ", trait, "\n")
    cat("\n\n=======================================\n")

    cat('\nExpanding dataset\n')
    cat(  '-----------------\n\n')
    cat("Normalizing\n")
    # let us normalize the data between 0 and 1
    # IMPORTANT
    #  in the first two, we take advantage of the fact we know the 
    #  dataset has 'Y' first and 'N' afterwards
    #nfy <- min_max_normalize(freqs[freqs$ab=='Y', ]$freqs)
    #nfn <- min_max_normalize(freqs[freqs$ab=='N', ]$freqs)
    #freqs <- cbind(freqs, norm.freqs=c(nfy, nfn))
    #  another way
    #nf <- Tapply(freqs$freqs ~ freqs$ab, min_max_normalize)
    #freqs <- cbind(freqs, norm.freqs=c(nf$Y, nf$N)
    #   and another one, that ensures each row gets its value, and
    #   plus, it does not add a new column if it already exists
    freqs$norm.freqs[freqs$ab == 'Y'] <- min_max_normalize(freqs$freqs[freqs$ab == 'Y'])
    freqs$norm.freqs[freqs$ab == 'N'] <- min_max_normalize(freqs$freqs[freqs$ab == 'N'])

    cat("Transforming\n")
    # Apply additional transformations
    # We can try to transform both the frequency and the normalized freq data
    #freqs <- cbind(freqs, log.freqs=log1p(freqs$freqs), sqrt.freqs=sqrt(freqs$freqs))
    #freqs <- cbind(freqs, log.norm.freqs=log1p(freqs$norm.freqs), 
    #                      sqrt.norm.freqs=sqrt(freqs$norm.freqs))
    # this also works but avoids column duplication (overwrites it if it exists)
    freqs$log.freqs=log1p(freqs$freqs)
    freqs$sqrt.freqs=sqrt(freqs$freqs)
    freqs$log.norm.freqs=log1p(freqs$norm.freqs)
    freqs$sqrt.norm.freqs=sqrt(freqs$norm.freqs)
    
    cat("Saving\n")
    # save data to speed up future analyses if needed
    fname <- outf(paste("R.freqs.by", trait, "Rds", sep='.'))
    saveRDS(freqs, fname)
    fname <- outf(paste("R.freqs.by", trait, "tab", sep='.'))
    write.table(freqs, file=fname, sep='\t')

    for (col in c("freqs", "log.freqs", "sqrt.freqs",
                  "norm.freqs", "log.norm.freqs", "sqrt.norm.freqs") )
    {

        cat('\n\nAnalyzing', col, "by", trait, '\n')
        cat('----------------------------------\n\n')
        #
      if (verbose) {
        numSummary(freqs[,col, drop=FALSE], groups=freqs$ab, statistics=c("mean", "sd", "IQR", "quantiles"), quantiles=c(0,
          .25,.5,.75,1))

        # there seems to be a slight difference in the means (119 vs. 115) and
        # a large difference in standard deviations (16.7 vs. 27.1)

        cat('\n    Normality tests\n')
        cat(  '    ---------------\n\n')
        # test for normality with a battery of tests
        for ( nt in c("shapiro.test","cvm.test","sf.test","ad.test","lillie.test","pearson.test") )
        {
            cat("\n        ", nt, "of", col, "by", trait, '\n')
            normalityTest(freqs[ , col] ~ freqs$ab, test=nt)
        }

        cat('\n    Plots\n')
        cat(  '    -----\n\n')
        # all seem to agree the data is strongly normally distributed, but this
        # is something I've often seen with Pearson distributions,
        # however the histogram of resistants seems to agree more with a Pearson
        # distribution than with a normal distribution:         
        cat("        Histogram of ", col, "by", trait, '\n')
        Hist(freqs[ ,col], groups=freqs$ab, 
	     scale="frequency", breaks=40, col="darkgray", 
             xlab=col, ylab="count",
             main=paste(col, "by", trait))
        Hist(freqs[ ,col], groups=freqs$ab, 
	     scale="percent", breaks=40, col="darkgray", 
             xlab=col, ylab="percent",
             main=paste(col, "by", trait))
        # and save to a file
        png(file=paste(outf('R.histogram-pct'), trait, col, 'png', sep='.'), 
            width=1024, height=1024)
        Hist(freqs[ ,col], groups=freqs$ab, 
	     scale="percent", breaks=40, col="darkgray", 
             xlab=col, ylab="percent",
             main=paste(col, "by", trait))
        dev.off()
        # and with 20 bars
        #png(file='R.histogram-pct-20bars', trait, col, 'png', sep='.'), 
        #    width=1024, height=1024)
        #Hist(freqs[ ,col], groups=freqs$ab, 
	#     scale="frequency", breaks=20, col="darkgray", 
        #     main=paste("Frequencies by", trait))
        #dev.off()
        # NOTE: the genomes seem to have a camel-shaped distribution,
        # suggestive of two different populations.
        cat("        Boxplot of ", col, "by", trait, '\n')
        par(mfrow=c(2,1))         
        # plot the freqs as reference and this data plot
        Boxplot(freqs$freqs ~ freqs$ab, id=list(method="y"), 
                xlab=trait, ylab="freqs", main=paste("freqs by", trait))
        Boxplot(freqs[ , col] ~ freqs$ab, id=list(method="y"), 
                xlab=trait, ylab=col, main=paste(col, "by", trait))
        par(mfrow=c(1,1))

        cat('\n    Checking variances\n')
        cat(  '    ------------------\n\n')
        # check homogeneity of the variances
        print(Tapply(freqs[ , col] ~ freqs$ab, var, na.action=na.omit)) # variances by group
        print(bartlett.test(freqs[ , col] ~ freqs$ab))
        print(Tapply(freqs[ , col] ~ freqs$ab, var, na.action=na.omit)) # variances by group
        print(leveneTest(freqs[ , col] ~ freqs$ab, center="median"))

        cat('\n    Checking statistical support\n')
        cat(  '    ----------------------------\n\n')
        # both say that variances are different, this makes the t-test dubious,
        # but all the same, we can violate one assumption, so... let us try
        # with an independent samples t-test with Welch correction for unequal
        # variances
        cat("        Parametric\n")
        print(t.test(freqs[ , col] ~ freqs$ab, 
               alternative='two.sided', conf.level=.95, var.equal=FALSE))

        # says they are different (p = 0.0123 with fasttree, 2.2e-16 with raxml), 
        # but this could be due to the difference in frequency numbers (we have 
        # more R genomes than S genomes. And indeed a non-parametric Wilcoxon test 
        # for independent samples (Mann-Whitney U) is non-significant for the
        # fasttree data (it is for the RaXML data)
        cat("        Medians\n")
        print(Tapply(freqs[ , col] ~ freqs$ab, median, na.action=na.omit)) # medians by group
	cat("        Non-parametric\n")
        wilcox.test(freqs[ , col] ~ freqs$ab, alternative="two.sided")
        #
        # NOTE that the skews that are apparent when looking at freqs will
        # disappear when using normalized data, for then they will move
        # between 0 and 1 and not min and max. To avoid that, one might use
        # the following (but will get then the same results as above, for if
        # we normalize all together we get the same distribution shape).
        #
      } # endif (verbose)  
        cat("\n\n    SUMMARY FOR", col, "by", trait, "\n")
        cat(    '    ----------------------------------------------\n')
        cat('\n    x =', levels(freqs$ab)[1], "y =", levels(freqs$ab)[2])
        cat('\n    unilateral tests are that x is <|> y\n')
        tryCatch( {
          cat("\n    Welch T-test (2 sided) =", 
              t.test(freqs[ ,col] ~ freqs$ab, 
                     alternative='two.sided', conf.level=.95,
                     var.equal=FALSE)$p.value)
          cat("\n    Welch T-test (less) =", 
              t.test(freqs[ ,col] ~ freqs$ab, 
                     alternative='less', conf.level=.95,
                     var.equal=FALSE)$p.value)
          cat("\n    Welch T-test (greater) =", 
              t.test(freqs[ ,col] ~ freqs$ab, 
                     alternative='greater', conf.level=.95,
                     var.equal=FALSE)$p.value)
          cat("\n    Mann-Whitney U (2 sided) =", 
              wilcox.test(freqs[ ,col] ~ freqs$ab, alternative='two.sided')$p.value)
          cat("\n    Mann-Whitney U (less) =", 
              wilcox.test(freqs[ ,col] ~ freqs$ab, alternative='less')$p.value)
          cat("\n    Mann-Whitney U (greater) =", 
              wilcox.test(freqs[ ,col] ~ freqs$ab, alternative='greater')$p.value)
          },
          error=say.and.continue 
        )
    }
}




cat("\n\n=======================================\n")
cat("\n    ANALYZING DETAILED DATA\n")
cat("\n\n=======================================\n")

# Now, analyze using ponderated data
per.branch.file <- "zyg_gubbins.per_branch_statistics.csv"

per.branch.stats <- read.table(per.branch.file, header=T, sep='\t')
pbs <- per.branch.stats

for (type in c("S83L", "D87N", "S83L.D87N", "resistant") ) {
    cat("\n\nAnalyzing by", type, '\n')
    pos.names <- metadata.df$ID[ metadata.df[ , type] == 'Y' ]
    neg.names <- metadata.df$ID[ metadata.df[ , type] == 'N' ]
    
    pbs.pos <- pbs[pbs$Node %in% pos.names, ]
    pbs.pos$yn <- rep('Y', length(pbs.pos$Node))
    pbs.neg <- pbs[ pbs$Node %in% neg.names, ]
    pbs.neg$yn <- rep('N', length(pbs.neg$Node))
    data <- rbind(pbs.pos, pbs.neg)
    data$YN <- as.factor(data$yn)
    
    cat("Computing log transforms")
    #compute log-transforms
    for (i in 2:9) {
        ncols <- dim(data)[2]
        nam <- colnames(data)[i]
        data[ , ncols+1] <- log1p(data[ , i])		# Log(1+x)
        newnam <- paste("log1p", nam, sep='.')
        print(paste(nam, '->', newnam))
        colnames(data)[ ncols+1] <- newnam
    }

    # Now, proceed with the analysis of raw and log-transformed data
    for (i in c(2:9, 13:20)) {
        col <- colnames(data)[i]
        cat('\n\n=================================================================')
        cat('\nAnalyzing "', col, "by", type, '"\n')
        cat('=================================================================\n\n')
      if (verbose) {
        cat("    Summary", col, "by", type, ":\n")
        cat("    ------------------------------------\n")
        tryCatch( {
          print(numSummary(data[ , col], group=data$yn))
          },
          error=say.and.continue
        )
        cat("\n\n    Normality", col, "by", type, " :\n")
        cat(    "    --------------------------------\n")
        tryCatch( {
          print(normalityTest(data[ ,i] ~ data$YN, test="shapiro.test"))
          print(bartlett.test(data[ ,i] ~ data$YN))
          },
          error=say.and.continue
        )

        cat("\n\n    Plots", col, "by", type, ":\n")
        cat(    '    ----------------------------\n')
        tryCatch( {
          Hist(data[ ,i], groups=data$YN, 
                     col="darkgray", scale="percent", breaks=40,
                     xlab=col, ylab="percent",
                     main=paste(col, "by", type))
          Boxplot(data[ ,i] ~ data$yn,
                        id=list(method="y"), 
                        xlab=type, ylab=col,
                        main=paste(col, "by",type))
          },
          error=say.and.continue 
        )

        cat("\n\n    Statistical support for", col, "by", type, ":\n")
        cat(    '    ----------------------------------------------\n')
        tryCatch( {
          print(t.test(data[ ,i] ~ data$YN, 
                       alternative='two.sided', conf.level=.95,
                       var.equal=FALSE))
          print(t.test(data[ ,i] ~ data$YN, 
                       alternative='less', conf.level=.95,
                       var.equal=FALSE))
          print(t.test(data[ ,i] ~ data$YN, 
                       alternative='greater', conf.level=.95,
                       var.equal=FALSE))
          print(Tapply(data[ ,i] ~ data$YN, median, na.action=na.omit))
          print(wilcox.test(data[ ,i] ~ data$yn, alternative='two.sided'))
          },
          error=say.and.continue 
        )
      } # endif (verbose)
        cat("\n\n    SUMMARY FOR", col, "by", type, "\n")
        cat(    '    ----------------------------------------------\n')
        cat('\n    x =', levels(data$YN)[1], "y =", levels(data$YN)[2])
        cat('\n    unilateral tests are that x is <|> y\n')
        tryCatch( {
          cat("\n    Welch T-test (2 sided) =", 
              t.test(data[ ,i] ~ data$YN, 
                     alternative='two.sided', conf.level=.95,
                     var.equal=FALSE)$p.value)
          cat("\n    Welch T-test (less) =", 
              t.test(data[ ,i] ~ data$YN, 
                     alternative='less', conf.level=.95,
                     var.equal=FALSE)$p.value)
          cat("\n    Welch T-test (greater) =", 
              t.test(data[ ,i] ~ data$YN, 
                     alternative='greater', conf.level=.95,
                     var.equal=FALSE)$p.value)
          cat("\n    Mann-Whitney U (2 sided) =", 
              wilcox.test(data[ ,i] ~ data$YN, alternative='two.sided')$p.value)
          cat("\n    Mann-Whitney U (less) =", 
              wilcox.test(data[ ,i] ~ data$YN, alternative='less')$p.value)
          cat("\n    Mann-Whitney U (greater) =", 
              wilcox.test(data[ ,i] ~ data$YN, alternative='greater')$p.value)
          cat("\n")
          },
          error=say.and.continue 
        )
    }
    rm(pos.names, neg.names, pbs.pos,pbs.neg, data)
}


# clustering
#kcf2 <- kmeans(freqs.rs$freqs, centers=2)
#kcf2$cluster
#kcf2$betweenss
#kcf2$totss
#kcf4 <- kmeans(freqs.rs$freqs, centers=4)
#kcf4$cluster


while ( ! is.null(dev.list()) ) dev.off()	# multi-plot PDF file
sink()
