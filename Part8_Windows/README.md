<img width="1684" height="330" alt="image" src="https://github.com/user-attachments/assets/772e93a2-499f-4727-bf1d-6485ca94601e" />We will cap off this workshop with a simple windowed FST analysis. This analysis looks at where population structure is concentrated across the genome. Unfortunately, RADSeq is quite coarse for it, but you can still learn important take-aways. At the very least, getting familiar with how to make simple manhattan plots is valuable.

First, you need to choose a two non-overlapping sets of samples to compare to eachother. This can be done in many ways. You can do it off of geography, off of where population structure is centered in your group, or off of phenotype. For the last one, Ethan may have a suggestion, so just ask. You can even use it to examine oddities in your PCA, such as groups splitting up in ways you don't expect (it could be a chromosomal inversion!).

For making your individual lists, just make files called "keep_X" and "keep_Y" for each X and Y you want to compare (e.g., keep_viti and keep_vanua for comparing two island groups, keep_white and keep_yellow for throat colors.\). Copy that list of individuals into the files. The more individuals the better.

Next, we will re-run populations again for the ingroup, this time keeping all variants.

    mkdir populations_out/windows
    populations --in-path stacks_out/ingroup/ --popmap popmap_ingroup \
                --out-path populations_out/windows --vcf
    mv populations_out/windows/populations.snps.vcf output/window_input.vcf

We will then run VCFtools on those variants for each set of comparisons you want to do (e.g., for each X and Y). If you do more than one comparison, just copy/paste the VCFtools command.

    #!/bin/bash
    
    vcftools --vcf output/window_input.vcf \
            --max-missing 0.75 --min-alleles 2 --max-alleles 2 \
            --weir-fst-pop keep_X --weir-fst-pop keep_Y \
            --fst-window-size 100000 \
            --out output/windowed_fst_X_Y_75

Unfortunately, this output is not accepted by the program we will use. Therefore, we will arbitrarily replace chromosome names with numbers. In practice, it is best to order these from largest to smallest chromosome, with unplaced scaffolds then the sex chromosome at the end. We won't worry about that for now, though.

#!/bin/bash

    grep -v CHROM output/windowed_fst_X_Y_75.windowed.weir.fst | \
        cut -f 1 | sort | uniq > chrom_list
    
    track=1
    for chrom in $(cat chrom_list); do
        sed -i "s/$chrom/$track/g" output/windowed_fst_X_Y_75.windowed.weir.fst
        track=$((track+1))
    done

Finally, we will use an R script for plotting. The only new package is qqman. It is primarily made up of a plotting function, so we can easily plot multiple different comparisons.

    install.packages("qqman")
    library(tidyverse)
    library(qqman)
    
    setwd('C:/Documents/PostDoc/Workshop_working_directory')
    
    plot_fst <- function(input, title_name="Windowed FST"){
      fst <- read.table(input, header=TRUE)
      fstsubset <- fst[complete.cases(fst),]
      SNP<-c(1:(nrow(fstsubset)))
      snp_df<-data.frame(SNP,fstsubset)
      manhattan(snp_df,chr="CHROM",bp="BIN_START",p="WEIGHTED_FST",
                snp="SNP",logp=FALSE, ylab="Weir and Cockerham Fst", 
                main=title_name, cex=1.5)
    }
    
    plot_fst("output/windowed_fst_X_Y_75.windowed.weir.fst", 
             "X vs Y Windowed FST")
    
