This section is all about making output. The slides discuss a lot of options, but this will describe some very simple outputs with populations we will use downstream.

The first and practically only step is to run populations! We will do it out of a script. This outputs a VCF and summary stats for a random SNP at 75% completeness.

    #!/bin/bash

    mkdir populations_out/ingroup
    populations --in-path stacks_out/ingroup/ --popmap popmap_ingroup \
            --out-path populations_out/ingroup --vcf --min-samples-overall 0.75 \
            --write-random-snp  --smooth --threads 2

Populations has pretty messy output, so I like to take the specific files I need and move them to a real output directory.

    mv populations_out/ingroup/populations.snps.vcf output/pachycephala_ingroup_75_thinned.vcf
    mv populations_out/ingroup/populations.sumstats_summary.tsv output/pachycephala_ingroup_75_thinned_sumstats.tsv

If we are running ahead of schedule, we will go through these output options and explore vcftools.
