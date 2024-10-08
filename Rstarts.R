library(ape)
library(tidyr)
library(RcmdrMisc)	# for Hist() and numSummary()
library(abind)
library(e1071)
library(Deducer)
library(effectsize)


#data.dir <- '0050'
#data.dir <- '0250'
data.dir <- '0500'
outdir  <- paste(data.dir, 'Rstats', sep='.')
metadata.file <- "zyg_metadata.tsv"
gubbins.gff.file <- "zyg_gubbins.recombination_predictions.gff"
ref.genome.gff.file <- "EcoliK12.gff3"
verbose <- TRUE



# Mark a p-value with stars to reflect its significance value
tag.p.value <- function(p.value=1) {
    return(paste(
            ifelse(p.value <= 0.05, "*", ' '),
            ifelse(p.value <= 0.01, "*", ' '),
            ifelse(p.value <= 0.001, "*", ' '),
            sep='')
    )
}

# Mark Cohen's D with stars to reflect it significance value
tag.cohen.d <- function(d=0) {
    return(paste(
            ifelse(d >= 0.2, '*', ' '),
            ifelse(d >= 0.5, '*', ' '),
            ifelse(d >= 0.8, '*', ' '),
            sep='')
          )
}

tag.hedge.g <- function(g=0) {
    return(paste(
            ifelse(g >= 0.2, '*', ' '),
            ifelse(g >= 0.5, '*', ' '),
            ifelse(g >= 0.8, '*', ' '),
            sep='')
          )
}

tag.glass.delta <- function(delta) {
    return(paste(
            ifelse(delta >= 0.2, '*', ' '),
            ifelse(delta >= 0.5, '*', ' '),
            ifelse(delta >= 0.8, '*', ' '),
            sep='')
          )
}


# execute a PLOT command and save the output as a png file
as.png <- function(PLOT=NULL, 
                   file='out.png', width=1024, height=1024, 
                   overwrite=FALSE) {
        
    if (is.null(file)) {
        tryCatch(print(PLOT))
    } else if ( overwrite | ! file.exists(file)  ) {
        tryCatch( {
                png(file, width=1000, height=1000)
                print(PLOT)
            },
            finally=dev.off()
        )
    }
    return ()
}


# set output file name
outf <- function(f) {
    return( paste(outdir, f, sep="/") )
}


# normalize values between minimum and maximum
min_max_normalize <- function(x) {
    (x - min(x)) / (max(x) - min(x))
}


# a simple convenience function
printf <- function(...) cat(sprintf(...))


# guess what? This function backups a file
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


# print a condition list and continue
say.and.continue <- function(condition=list()) {
    cat('\n\n    *** ERROR ***\n')
    for (n in names(condition))
        cat(    '    ***', n, '=', condition[[n]], '***\n')
    cat('\n')
    return()
}


# Save all output
if (verbose) {
    logfile <- outf("R.Rstats.log")
} else {
    logfile <- outf("R.Rstats.sum")
}


# backup log file brefore opening a new one
make.backup(logfile)
sink(logfile, split=TRUE)

#
# Prepare the data we need
#
cat('Building statistical data bases\n')
cat('===============================\n\n')

## first we start by creating empty data.frames with appropriate column names

# resistance status
all.freqs.res <- data.frame(matrix(ncol=3, nrow=0))
names(all.freqs.res) <- c('genome', 'freqs', 'ab')

# S83L mutant status
all.freqs.s83l <- data.frame(matrix(ncol=3, nrow=0))
names(all.freqs.s83l) <- c('genome', 'freqs', 'ab')

# D87N mutant status
all.freqs.d87n <- data.frame(matrix(ncol=3, nrow=0))
names(all.freqs.d87n) <- c('genome', 'freqs', 'ab')

# double mutant status
all.freqs.s83l.d87n <- data.frame(matrix(ncol=3, nrow=0))
names(all.freqs.s83l.d87n) <- c('genome', 'freqs', 'ab')



# Apply additional transformations
transform.freqs <- function(freqs) {
    cat('\n    Expanding dataset\n')
    # Apply additional transformations
    # We can try to transform both the frequency and the normalized freq data
    #freqs <- cbind(freqs, log.freqs=log1p(freqs$freqs), sqrt.freqs=sqrt(freqs$freqs))
    #freqs <- cbind(freqs, log.norm.freqs=log1p(freqs$norm.freqs), 
    #                      sqrt.norm.freqs=sqrt(freqs$norm.freqs))
    # this also works but avoids column duplication (overwrites it if it exists)
    freqs$log.freqs=log1p(freqs$freqs)
    freqs$sqrt.freqs=sqrt(freqs$freqs)
    return (freqs)
}


# Analyze data
analyze.freqs <- function(freqs) {
    all.pvals <- data.frame(dummy='')

    # for the given columns
    for ( co in c("freqs", "log.freqs", "sqrt.freqs") ) {
        cat('\n\n    Analyzing', co, '\n')
        x <- levels(freqs$ab)[1]
        y <- levels(freqs$ab)[2]
        t.b.name <- paste('p.', co, '.t.', x, 'ne', y, sep='') # t test
        t.l.name <- paste('p.', co, '.t.', x, 'lt', y, sep='')
        t.g.name <- paste('p.', co, '.t.', x, 'gt', y, sep='')
        w.b.name <- paste('p.', co, '.U.', x, 'ne', y, sep='') # Mann-Whitney U
        w.l.name <- paste('p.', co, '.U.', x, 'lt', y, sep='')
        w.g.name <- paste('p.', co, '.U.', x, 'gt', y, sep='')
        cat('\n    x =', x, "y =", y)
        cat('\n    unilateral tests are that x is <|> y\n')
        tryCatch( {
          t.test.b <- t.test(freqs[ ,co] ~ freqs$ab, 
                     alternative='two.sided', conf.level=.95,
                     var.equal=FALSE)
          cat("\n    Welch T-test (", x, '≠', y, ") =", 
              t.test.b$p.value)
        }, error=say.and.continue )
        tryCatch( {
          t.test.l <- t.test(freqs[ ,co] ~ freqs$ab, 
                     alternative='less', conf.level=.95,
                     var.equal=FALSE)
          cat("\n    Welch T-test (", x, '<', y, ") =", 
              t.test.l$p.value)
        }, error=say.and.continue )
        tryCatch( {
          t.test.g <- t.test(freqs[ ,co] ~ freqs$ab, 
                     alternative='greater', conf.level=.95,
                     var.equal=FALSE)
          cat("\n    Welch T-test (", x, '>', y, ") =", 
              t.test.g$p.value)
        }, error=say.and.continue )
        tryCatch( {
          w.test.b <- wilcox.test(freqs[ ,co] ~ freqs$ab, alternative='two.sided')
          cat("\n    Mann-Whitney U (", x, '≠', y, ") =", 
              w.test.b$p.value,
              ifelse(w.test.b$p.value < 0.05, '\t+', '\t-'))
        }, error=say.and.continue )
        tryCatch( {
          w.test.l <- wilcox.test(freqs[ ,co] ~ freqs$ab, alternative='less')
          cat("\n    Mann-Whitney U (", x, '<', y, ") =", 
              w.test.l$p.value,
              ifelse(w.test.l$p.value < 0.05, '\t+', '\t-'))
        }, error=say.and.continue )
        tryCatch( {
          w.test.g <- wilcox.test(freqs[ ,co] ~ freqs$ab, alternative='greater')
          cat("\n    Mann-Whitney U (", x, '>', y, ") =", 
              w.test.g$p.value,
              ifelse(w.test.g$p.value < 0.05, '\t+', '\t-'))
        }, error=say.and.continue )
        cat('\n')
        pvals <- data.frame(
    		    t.test.b$p.value, 
                    t.test.l$p.value, 
                    t.test.g$p.value,
                    w.test.b$p.value,
                    w.test.l$p.value,
                    w.test.g$p.value
                    )
        names(pvals) <- c(t.b.name, t.l.name, t.g.name, w.b.name, w.l.name, w.g.name)
        #print(all.pvals)
        #print('')
        #print(pvals)
        all.pvals <- cbind(all.pvals, pvals)
    }
    return(all.pvals[,-1])
}
ga

#
# Now we are ready
# Procceed with the analysis
#

n.exp <- 0
for ( sub.data in dir(data.dir) ) { 
    n.exp <- n.exp + 1
    # for all the other data sets but the first one which is already loaded
    cat('\n\n============================\n       ', sub.data, '\n')
    cat('============================\n')
    metadata.f <- paste(data.dir, sub.data, metadata.file, sep='/')
    gubbins.gff.f <- paste(data.dir, sub.data, gubbins.gff.file, sep='/')
    ref.genome.gff.f <- paste(data.dir, sub.data, ref.genome.gff.file, sep='/')

    cat("    Reading Gubbins output\n")
    gubbins.GFF <- readRDS(paste(data.dir, sub.data, "zyg_R/R.gubbins.GFF.Rds", sep='/'))

    cat("    Reading recombination events\n")
    rec.events <- readRDS(paste(data.dir, sub.data, "zyg_R/R.rec.events.Rds", sep='/'))
    
    cat("    Reading phylogenetic tree\n")
    tree <- readRDS(paste(data.dir, sub.data, "zyg_R/R.tree.Rds", sep='/'))
    
    cat("    Reading metadata\n")
    metadata.df <- readRDS(paste(data.dir, sub.data, "zyg_R/R.metadata.df.Rds", sep='/'))
    metadata <- tidyr::as_tibble(read.table(metadata.f, header=T, sep='\t', comment.char="?"))
    
    cat("    Reading reference genome\n")
    refgenome.GFF <- readRDS(paste(data.dir, sub.data, "zyg_R/R.refgenome.GFF.Rds", sep='/'))

    cat("    Reading rec frequency per base\n")
    rec.freq <- readRDS(paste(data.dir, sub.data, "zyg_R/R.rec.freq.Rds", sep='/'))
    
    cat("    Reading rec events per genome\n")
    rec.genome <- readRDS(paste(data.dir, sub.data, "zyg_R/R.rec.genome.Rds", sep='/'))
    
    genome.recombinations <- data.frame(rec=rec.genome)
    colnames(genome.recombinations) <- c('genome', 'rec.Freq')
    head(genome.recombinations)

    
    cat('\n\n    Preparing datasets\n')
    cat(    '    ------------------\n\n')
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
    cat("    By S83L status\n")
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

    cat("    By D87N status\n")
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

    cat("    By S83L+D87N status\n")
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

    # Now add the data to the global data
    all.freqs.res <- rbind(all.freqs.res, freqs.res)
    all.freqs.s83l <- rbind(all.freqs.s83l, freqs.s83l)
    all.freqs.d87n <- rbind(all.freqs.d87n, freqs.d87n)
    all.freqs.s83l.d87n <- rbind(all.freqs.s83l.d87n, freqs.s83l.d87n)

    # At this point we may consider doing partial statistics and
    # saving them for a bootstrap
     cat("    resistance\n")
     freqs.res <- transform.freqs(freqs.res)
     r.p <- analyze.freqs(data.frame(freqs.res))
     names(r.p) <- paste("resis", names(r.p), sep='.')
     cat("    S83L\n")
     freqs.s83l <- transform.freqs(freqs.s83l)
     s.p <- analyze.freqs(freqs.s83l)
     names(s.p) <- paste("S83L", names(s.p), sep='.')
     cat("    D87N\n")
     freqs.d87n <- transform.freqs(freqs.d87n)
     d.p <- analyze.freqs(freqs.d87n)
     names(d.p) <- paste("D87N", names(d.p), sep='.')
     cat("    S83LD87N\n")
     freqs.s83l.d87n <- transform.freqs(freqs.s83l.d87n)
     s.d.p <- analyze.freqs(freqs.s83l.d87n)
     names(s.d.p) <- paste("S83L+D87N", names(s.d.p), sep='.')
     
     pvals <- data.frame(r.p, s.p, d.p, s.d.p)
     if (n.exp == 1)
         boot.pvals <- pvals
     else
         boot.pvals <- rbind(boot.pvals, pvals)

}


# now count significants and attempt a pseudo-bootstrap
# =====================================================
cat('\n\n==================================================\n')
cat(    '               PSEUDO - BOOTSTRAP                 \n')
cat(    ' ANALYZING FREQUENCIES OF SIGNIFICANT DIFFERENCES \n')
cat(    '==================================================\n\n')
i <- 0
for ( co in colnames(boot.pvals) ) {
    i <- i + 1
    p.vals <- boot.pvals[ , co]
    n.vals <- length(p.vals)
    n.sig <- sum(p.vals <= 0.05)
    n.uns <- sum(p.vals > 0.05)
    pct.sig <- (n.sig / length(p.vals)) * 100
    pct.uns <- (n.uns / length(p.vals)) * 100
    mn.t <- mcnemar.test(
            matrix(c(n.sig, n.uns, 0, n.vals), nrow=2, byrow=T),
            correct=T)
    ml.t <- likelihood.test(
            matrix(c(n.sig, n.uns, 0, n.vals), nrow=2, byrow=T),
            conservative=T)
#  if (n.sig > n.uns)
#if (i > j - 6)
    cat(co, " ",
        n.sig, " p ≤ 0.05 (", pct.sig, "%)\t",
        #n.uns, " p > 0.05 (", pct.uns, "%)\t",
        #'p(MN)=', mn.t$p.value, ' ', tag.p.value(mn.t$p.value), 
        #'\t',
        'p(ML)=', ml.t$p.value, ' ', tag.p.value(ml.t$p.value),
        '\n',
        sep='')
#if (i == j) break
    if (! i %% 6) cat('\n')
}



cat('\n\n==================================================\n')
cat(    'COMPUTING GLOBAL STATISTICS (ALL GENOMES TOGETHER)\n')
cat(    '==================================================\n\n')

# Compute stattistics using all unique genomes together
# =====================================================
# remove duplicates
uni.freqs.res <- unique(all.freqs.res, by='genome')
uni.freqs.s83l <- unique(all.freqs.s83l, by='genome')
uni.freqs.d87n <- unique(all.freqs.d87n, by='genome')
uni.freqs.s83l.d87n <- unique(all.freqs.s83l.d87n, by='genome')
cat('\nUnique sizes:\n')
cat('uni.freqs.res', dim(uni.freqs.res), '\n')
cat('uni.freqs.s83l', dim(uni.freqs.s83l), '\n')
cat('uni.freqs.d87n', dim(uni.freqs.d87n), '\n')
cat('uni.freqs.s83l.d87n', dim(uni.freqs.s83l.d87n), '\n')

freqs.list <- list(
	resistance=uni.freqs.res,
        S83L=uni.freqs.s83l,
        D87N=uni.freqs.d87n,
        S83L.D87N=uni.freqs.s83l.d87n
        )
fname <- outf(paste("data/R.freqs.Rds", sep='.'))
saveRDS(freqs.list, fname)
fname <- outf(paste("data/R.freqs.Rdata", sep='.'))
save(freqs.list, file=fname)



#
# ANALYZE THE FREQS DATASETS
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
    fname <- outf(paste("data/R.freqs.by", trait, "Rds", sep='.'))
    saveRDS(freqs, fname)
    fname <- outf(paste("data/R.freqs.by", trait, "tab", sep='.'))
    write.table(freqs, file=fname, sep='\t')

    for (col in c("freqs", "log.freqs", "sqrt.freqs"#,
                  #"norm.freqs", "log.norm.freqs", "sqrt.norm.freqs"
                  ) )
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
            for ( nt in c("shapiro.test","cvm.test","sf.test",
                          "ad.test","lillie.test","pearson.test") )
            {

                # for very large datasets, it is better to use a qqplot
                # and inspect it visually
                out.png <- outf(paste('plots/R.qqplot', col, '.png', sep=''))
                as.png(
                    qqnorm(freqs[ , col], main=paste(col, "Normal Q-Q Plot")),
                    out.png)

	        # we cannot use shapiro or sf with too large samples
                if ( (dim(freqs)[1] >= 5000) & (nt == "shapiro.test" | nt == "sf.test") ) { 
                    next
                }
                # for too large samples, this is indeed useless, as they
                # will always fail the normality test as soon as there is
                # any minimal deviation
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
            out.png <- outf( paste('plots/R.histogram-pct', trait, col, 'png', sep='.'))
            as.png(
                Hist(freqs[ ,col], groups=freqs$ab, 
	             scale="percent", breaks=40, col="darkgray", 
                     xlab=col, ylab="percent",
                     main=paste(col, "by", trait)),
                out.png)
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
            out.png <- outf(paste('plots/R.Boxplot', trait, col, 'png', sep='.'))
            as.png( {
                par(mfrow=c(2,1))         
                # plot the freqs as reference and this data plot
                Boxplot(freqs$freqs ~ freqs$ab, id=list(method="y"), 
                        xlab=trait, ylab="freqs", main=paste("freqs by", trait))
                Boxplot(freqs[ , col] ~ freqs$ab, id=list(method="y"), 
                        xlab=trait, ylab=col, main=paste(col, "by", trait))
                par(mfrow=c(1,1))
                },
                out.png)
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
            tryCatch( {
              cat("        Parametric\n")
              tst <- (t.test(freqs[ , col] ~ freqs$ab, 
                     alternative='two.sided', conf.level=.95, var.equal=FALSE))
              cat("        Welch Two Sample t-test\n")
              cat("        p.value",  
                  tst$p.value, tag.p.value(tst$p.value), '\n')
              cat("        95% confidence interval",
                           tst$conf.int[1], '..', tst$conf.int[2], '\n')
              cat("        means (", levels(freqs$ab)[1], levels(freqs$ab)[2], ")",
                           tst$estimate[1], tst$estimate[2], '\n')
            }, error=say.and.continue )
            # this could be due to the difference in frequency numbers (if we have 
            # more R genomes than S genomes.

            tryCatch( {
              cat("        Non-parametric\n")
              tst <- wilcox.test(freqs[ , col] ~ freqs$ab, alternative="two.sided")
              cat("        Wilcoxon rank sum test with continuity correction\n")
              cat("        p.value",  
                  tst$p.value, tag.p.value(tst$p.value), '\n')
              medians <- Tapply(freqs[ , col] ~ freqs$ab, median, na.action=na.omit)
              cat("        medians (",levels(freqs$ab)[1], levels(freqs$ab)[2], ")",
                  medians[1], medians[2], '\n')
            }, error=say.and.continue )
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
          tst <- t.test(freqs[ ,col] ~ freqs$ab, 
                     alternative='two.sided', conf.level=.95,
                     var.equal=FALSE)
          cat("\n    Welch T-test (2 sided) =", 
              tst$p.value, tag.p.value(tst$p.value))
        }, error=say.and.continue )
        tryCatch( {
          tst <- t.test(freqs[ ,col] ~ freqs$ab, 
                     alternative='less', conf.level=.95,
                     var.equal=FALSE)
          cat("\n    Welch T-test (less) =", 
              tst$p.value, tag.p.value(tst$p.value))
        }, error=say.and.continue )
        tryCatch( {
          tst <- t.test(freqs[ ,col] ~ freqs$ab, 
                     alternative='greater', conf.level=.95,
                     var.equal=FALSE)
          cat("\n    Welch T-test (greater) =", 
              tst$p.value, tag.p.value(tst$p.value))
        }, error=say.and.continue )
        tryCatch( {
          tst <- wilcox.test(freqs[ ,col] ~ freqs$ab, alternative='two.sided')
          cat("\n    Mann-Whitney U (2 sided) =", 
              tst$p.value, tag.p.value(tst$p.value))
        }, error=say.and.continue )
        tryCatch( {
          tst <- wilcox.test(freqs[ ,col] ~ freqs$ab, alternative='less')
          cat("\n    Mann-Whitney U (less) =", 
              tst$p.value, tag.p.value(tst$p.value))
        }, error=say.and.continue )
        tryCatch( {
          tst <- wilcox.test(freqs[ ,col] ~ freqs$ab, alternative='greater')
          cat("\n    Mann-Whitney U (greater) =", 
              tst$p.value, tag.p.value(tst$p.value))
          }, error=say.and.continue )
        # effectsize provides functions for estimating the common indices 
        # of standardized differences such as Cohen’s d (cohens_d()), 
        # Hedges’ g (hedges_g()) for both paired and independent samples 
        # (Cohen 1988; Hedges and Olkin 1985), and Glass’ Δ (glass_delta()) 
        # for independent samples with different variances (Hedges and 
        # Olkin 1985).
        tryCatch( {
          cd <- cohens_d(freqs[ , col] ~ freqs$ab)
          cat("\n    Cohen's D =", cd$Cohens_d, 
              "[", cd$CI_low, ",", cd$CI_high, "]", 
              tag.cohen.d(abs(cd$Cohens_d))
             )
        }, error=say.and.continue )
        # The only difference between Cohen’s d and Hedges’ g is that 
        # Hedges’ g takes each sample size into consideration when 
        # calculating the overall effect size.
        #
        # Thus, it’s recommended to use Hedge’s g to calculate effect 
        # size when the two sample sizes are not equal.
        tryCatch( {
          hg <- hedges_g(freqs[ , col] ~ freqs$ab)
          cat("\n    Hedge's G =", hg$Hedges_g, 
              "[", gd$CI_low, ",", gd$CI_high, "]", 
              tag.hedge.g(abs(hg$Hedges_g))
             )
        }, error=say.and.continue )
        tryCatch( {
          gd <- glass_delta(freqs[ , col] ~ freqs$ab)
          cat("\n    Glass's Delta =", gd$Glass_delta, 
              "[", gd$CI_low, ",", gd$CI_high, "]", 
              tag.glass.delta(abs(gd$Glass_delta)),
              interpret_cohens_d((abs(gd$Glass_delta)), rules="gignac2016") # or "cohen1988"
             )
        }, error=say.and.continue )
        cat('\n')
    }
}


sink()
