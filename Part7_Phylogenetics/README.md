Now we will construct a phylogeny for our dataset. The trickiest part of this is making the input, because Stacks doesn't normally output sequence data per-individual. To "trick" it, we will add an "Individual" column to our sample table, with names at least 6 and at most 10 characters long. Then, we will make a new population map.

    cut -f 1,7 -d ',' sample_table_Pachycephala.csv | \
    	tail -n+2 | \
    	tr , \\t > popmap_individual

Next, we will make a new bash script with nano (e.g., nano stacks_phylogenetics.sh) which will have contents like:

    #!/bin/bash

    # run gstacks
    gstacks -I bams -M popmap_individual -O stacks_out/individual/

    # make a new directory
    mkdir populations_out/individual
    # run populations
    populations --in-path stacks_out/individual/ --popmap popmap_individual \
                --out-path populations_out/individual --phylip-var-all \
                --min-samples-overall 0.75
    # move and rename the output
    grep -v "#" populations_out/individual/populations.all.phylip > output/NAME_75.phylip

Finally, we will run IQTree. But first, we need to install it into our conda environment (make sure the workshop environment is activated!).

    conda install -c bioconda -c conda-forge iqtree

Now run IQTree, here with a full model selection protocol, not reccomended for very large datasets.

    iqtree -s output/NAME_75.phylip -m MFP -B 1000 -T 2

We will now add concordance factors to this tree, which tell us how well our data really support parts of the tree.

    iqtree -t output/NAME_75.phylip.treefile \
           --scf 5000 -s output/NAME_75.phylip \
           --prefix NAME_SCF_75 -T 2

We will view our trees in PearTree: https://peartree.live/. It can be a bit tricky to figure out, but try playing around with the software to learn it. Make sure to save whatever you get.
