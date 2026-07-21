This section is dedicated to two visualization methods for population structure analyses, genomic PCAs and admixture plots. Both of these will be run primarily in R. The instructions are designed to be one or two files that are run piecemeal.

## PCA

First, we will load our packages, adegenet and vcfR are the most critical here.

    # install relevant packages
    install.packages(c("adegenet", "ade4", "vcfR", 
                       "parallel", "viridis", "tidyverse"))
    
    # load packages
    library("adegenet")
    library("ade4")
    library("vcfR")
    library("parallel")
    library("tidyverse")
    library("viridis")

Next we will set our working directory. This will be the path to your workshop directory. Once you are there, you will read in your ingroup only VCF.

    # set working directory
    setwd('C:/Documents/PostDoc/Workshop_working_directory')
    
    # read the VCF in as a vcfR object
    ingroup_vcf <- read.vcfR("output/NAME_ingroup_75_thinned.vcf")

You will then convert that VCF to a "genlight" object and assign population identifiers.

    # convert the VCF to a genlight file
    ingroup_genlight <- vcfR2genlight(ingroup_vcf)
    
    # assign populations to the genlight based population table
    pop(ingroup_genlight) <- filter(read.csv("sample_table_NAME.csv"), 
                                    Island!="Outgroup")$Island

Finally, you will run the PCA itself. This is parallelized (doesn't need to be for these small datasets) and we will retain 4 principal components. You can re-run with more if you want to explore further.

    # run the PCA using 2 cores of your PC, and keep 4 PCs
    ingroup_pca <- glPca(ingroup_genlight, n.cores=2, nf=4)
    
We will now plot the output of that PCA. This is the "naive" plot, only the first two PCs. Note that the first argument in "magma" should be approximately the number of islands in your dataset.

    # make a plot using the s.class program
    s.class(ingroup_pca$scores[,c(1,2)], 
            pop(ingroup_genlight), col=magma(10, begin=.8, end=0),
            clab=1, cstar=1, cell=2.5)

Now we will change the first line in the PCA arguments to plot different sets of principal components.

    # Let's try it again with PC 3 and 4, then 1 and 3
    s.class(ingroup_pca$scores[,c(3,4)],
            pop(ingroup_genlight), col=magma(10, begin=.8, end=0),
            clab=1, cstar=1, cell=2.5)
    s.class(ingroup_pca$scores[,c(1,3)],
            pop(ingroup_genlight), col=magma(10, begin=.8, end=0),
            clab=1, cstar=1, cell=2.5)

The utility of each of these plots can be visualized by a barplot of how much variance each principal component explains.

    # barplot of the eigenvalues to tell you the percent of variance each PC explains
    barplot(ingroup_pca$eig/sum(ingroup_pca$eig), 
            main="Variance explained per PC",
            col=heat.colors(length(ingroup_pca$eig)))

## sNMF

The sNMF algorithm is a way of computing admixture fractions, a rough way of detecting hybridization. It comes with MANY caveats, but it is important to know how to conduct as a first pass for detecting hybrids. You can either make a new script of append this to the bottom of your PCA script.

First, we will install and load our packages (including two sNMF-specific functions).

    # install LEA
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("LEA")
    
    # Load in additional tools for sNMF
    source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
    source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")
    
    # Load main library, LEA, and a color scheme
    library(LEA)
    library(RColorBrewer)

Now set your working directory and make a .geno file from your VCF.

    # set working directory
    setwd('C:/Documents/PostDoc/Workshop_working_directory')
    
    # vcf2geno function from LEA make a geno file in a specified directory
    vcf2geno("output/NAME_ingroup_75_thinned.vcf", # input path 
             output = "output/ NAME _ingroup_75_thinned.geno") # output path

Now run sNMF! Here I run if for up to 7 populations, which may be excessive for you. Note that the alpha parameter is a unique thing to sNMF's algorithm that lets you variably punish admixture.
    
    # run sNMF for 1-7 pops, 20 reps per pop, modest admixture punish (alpha)
    ingroup_snmf = snmf("output/pachycephala_ingroup_75_thinned.geno", ploidy=2, 
                      K = 1:7,  alpha = 100, project = "new", entropy = T,  repetitions = 20)

Now plot the cross-validation criterion. The "ideal" number of populations minimizes this number, but I reccomend taking that with a strong grain of salt. Biologically meaningful population structure can occur at higher K values, and sometimes higher values are just noise.

    # plot to decide optimal K, helpful but not absolute
    plot(ingroup_snmf, cex = 1.2, col = "lightblue", pch = 19)

Now we will plot the actual data! First we will make a color palette. I like to use gradient palettes for this, but vary it depending on what is most effective. Change your colors to something that fits your species!

    # make a color palette
    BlackYellow <- colorRampPalette(c("yellow", "black"))

Now make a function that will actually plot the sNMF output for a given number of populations (k value). The palette is used as an argument here, so make sure to adjust that accordingly.

    # make function for plotting admixture output
    plot_sNMF <- function(input, k_val, colors = BlackYellow(k_val+1)){
      # picks best run best on cross entropy
      best_run <- which.min(cross.entropy(input, K = k_val))
      # makes q matrix of ancestry coeffs
      q_matrix <- Q(input, K = k_val, run = best_run)
      # plots the output, space makes blank between indivs
      barplot(t(q_matrix), col = colors, border = NA, space = 0.25, xlab = "Individuals", ylab = "Admixture coefficients", horiz=FALSE)
    }

Now plot as many values of K as makes sense. Reduce this number (and the value in mfrow) based on what makes sense for your dataset.
    
    # plotting parameters, note the number in mfrow should be your k
    par(mfrow=c(7,1), mar=c(0,2,0,0), oma=c(1,2,1,0))
    
    # plot sNMF for each K
    plot_sNMF(ingroup_snmf, 1)
    plot_sNMF(ingroup_snmf, 2)
    plot_sNMF(ingroup_snmf, 3)
    plot_sNMF(ingroup_snmf, 4)
    plot_sNMF(ingroup_snmf, 5)
    plot_sNMF(ingroup_snmf, 6)
    plot_sNMF(ingroup_snmf, 7)

Note that ordering the sNMF input can be tricky. It defaults to the same order as the VCF. You can change it manually by going back to command line to make a new sorted VCF. Like this:

Make sure bcftools is installed (it should be in your env)
Compress vcf: bgzip in.vcf
Index vcf: tabix in.vcf.gz
Order vcf: bcftools view -S sorted_list in.vcf.gz > out.vcf
Decompress vcf: gunzip in.vcf.gz
