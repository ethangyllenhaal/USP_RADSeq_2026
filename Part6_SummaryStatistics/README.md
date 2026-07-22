Summary statistics are a very powerful way to describe population genetic data, which can be used for gaining a strong, quantitative sense of what is going on. It can take a while to get a good handle on what they mean for a given group, but once you understand them they can be the fastest way to understanding what is going on with your group of interest.

This section will start by opening the summary statistics file output by populations and selecting columns of interest that relate to genetic diversity. Then, add a column denoting island sizes for your focal taxa based on a quick search. If you want to do this more objectively, I STRONGLY reccomend PleistoDist by collaborator Daivd Tan: https://davidbirdtan.com/pleistodist/introduction/

## Correlation plot of genetic diversity

Once your file is obtained, we will look at how our statistics of interest are correlated with each other. We will use tidyverse for data processing (and future plotting) and corrplot for visualizing multiple regression coefficients.

    # we will only need the tidyverse and corrplot
    library(tidyverse)
    library(corrplot)
    
    # change to working directory
    setwd('C:/Documents/PostDoc/Workshop_working_directory')

We will then read in our file and make a correlation matrix for all columns except your Population column (or whatever you chose to name it).

    # read in the file
    sumstats <- read.csv("summary_stats.csv")
    
    # make a correlation matrix, removing the Population column
    corr_matrix <- cor(select(sumstats, !Population))

Now make two corrplots. First a simple one, then one that is more complex, showing precise values on the upper axis and rough visualizations of the correlations on the lower.

    # look at a correlation plot of all values
    corrplot(corr_matrix)
    
    # more complex
    corrplot.mixed(corr_matrix, order = 'AOE', tl.cex=0.8, tl.pos="lt",
             lower="ellipse", # lower corner is a visual depiction
             upper="number") # upper is the correlation coefficient

## Island size and genetic diversity

Now, chose a diversity statistic and see how strongly it is correlated with island size. Do this both visually and with a linear model

    ggplot(data=sumstats, aes(x=log(Size), y=Pi)) +
      geom_point() + 
      geom_smooth(method="lm") +
      theme_bw()

    model = lm(data=sumstats, Pi~log(Size))
    summary(model)

## Quantifying divergence with FST

We will now make a heatmap of pairwise FST values for all islands in our dataset! This will show us where divergence stands out. This method inspired by Devon DeRaad: https://devonderaad.github.io/zosterops.rad/zosterops.pairwise.fst.html. I tend to do this with a chain of commandline tools, in case alternatives are of interest: https://github.com/ethangyllenhaal/FijiPachyRad/tree/main/04_pairwise_Fst.

This is intended to be added to the same script as before. Start by loading and as needed installing these lbiraries.

    # load packages for FST
    library(vcfR)
    library(adegenet)
    library(StAMPP)
    library(reshape)

Next, read in your data

    # read in VCF
        ingroup_vcf <- read.vcfR("output/pachycephala_ingroup_75_thinned.vcf")
    
    # read in sample table
    samples <- filter(read.csv("sample_table_pachycephala.csv", sep=","), # read in the CSV
                  Island!="Outgroup") # remove outgroup

Now, just like in the PCA, convert it to a genlight and assign populations.

    # make genlight and set population like PCA
    ingroup_genlight <- vcfR2genlight(ingroup_vcf)
    ingroup_genlight@pop <- as.factor(samples$Island)

Now calculate the FST, convert it to a matrix, and finally use reshape to convert it to long format for plotting (check out the difference, this is key for making effective analyses).

    # calculate FST
    pairwise_fst <- stamppFst(ingroup_genlight)
    
    # convert to matrix, fill in both columns
    fst_matrix <- pairwise_fst$Fsts
    fst_matrix[upper.tri(fst_matrix)] <- t(fst_matrix)[upper.tri(fst_matrix)]
    
    # melt to make proper heatmap input
    heat_input <- reshape::melt(fst_matrix)

Finally, plot it with ggplot.

    # plot with ggplot
    ggplot(data = heat_input, aes(x=X1, y=X2, fill=value)) + 
      geom_tile() + # tile plot
      geom_text(data=heat_input,aes(label=round(value, 2))) + # add values
      theme_minimal() + # remove labels
      scale_fill_gradient2(low = "white", high = "red", space = "Lab", name="Fst") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), # axis labels
            axis.text.y = element_text(angle = 45, hjust = 1))
