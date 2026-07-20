<img width="513" height="80" alt="image" src="https://github.com/user-attachments/assets/f4438428-fe68-46d5-afb5-ff2cc4721963" />This section focuses on learning command line. It starts by exploring a few standard commands. I reccomend working with a cheatsheet like this as well: https://linuxstans.com/bash-cheat-sheet/

Start by navigating your computer

    # pwd lists out current directory
    pwd

    # ls lists the files in a directory
    ls

    # the -l option (sometimes coded as ll) is very helpful for getting detailed information on size, edit time, etc
    ls -l
    ll

Now we will move around and edit the directory structure

    # cd lets you change your directory
    # mkdir lets you make a new directory
    
    # make a new directory for this workshop, then change to it
    mkdir workshop_dir
    cd workshop_dir

    # make subdirectories, you can chain them together
    mkdir reads bams stacks_out populations_out output

We will use miniconda for all command line (non R) program installations. Install miniconda for linux here: https://www.anaconda.com/docs/getting-started/miniconda/install/linux-install

    curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh

We will install everything in one command!

    conda create -n rad-workshop-env -c bioconda -c conda-forge bwa stacks vcftools samtools

And activate it with

    conda activate rad-workshop-env

While in person people download reads, they will make their first script, the comp sci 101 classic, Hello World! Make a new script like.

    touch hello_world.sh

Edit that file with nano

    nano hello_world.sh

The contents will be a header that tells the computer how to run it and a line to print a phrase

    #!/bin/bash
    echo Hello world!

Write out (save) with control+O, exit with control+X. Then, run it!

    sh hello_world.sh

Open your sample list file in a spreadsheet program or a command line viewer like less. You will make a sample list of the first column, cutting out the first column, then removing the first line, and writing it out. At the end, view it with less.

    cut -f 1 -d ',' sample_table_NAME.csv |
        tail -n+2 >
        sample_list_NAME

    less sample_list_NAME

Now make a population map the same way.

    cut -f 1,4 -d ',' sample_table_NAME.csv | 
        tail -n+2 | 
        tr , \\t > 
        popmap_full

Name subset it to the ingroup using grep.

    grep Outgroup popmap_full
    grep –v Outgroup popmap_full
    grep –v Outgroup popmap_full > popmap_ingroup
