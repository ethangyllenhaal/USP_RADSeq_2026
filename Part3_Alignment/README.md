This is the computationally intensive part of the workshop. We will look at reads, align them, then call SNPs with gstacks.

Let's start with counting reads. You can get four times the number of reads like this:

    zcat reads/sample_name.*gz | wc -l

Now you can do it in a loop in a bash script (a bit basic, there are ways of doing the math and making it a single line per sample we can cover if we have time):

    #!/bin/bash

    # write a new file
    > read_counts.txt

    for sample in $(cat sample_list_pachycephala); do
         echo $sample >> read_counts.txt
         zcat reads/${sample}.*gz | wc -l >> read_counts.txt
    done

Next we will index the reference we will align reads to

    bwa index name_of_reference.fna

And we will now align all of our reads to the reference!

    for sample in $(cat sample_list_NAME); do
        bwa mem -t 2 -M name_of_reference.fna \
        reads/${sample}.*.gz \
        | samtools view -b -F 4 \ 
        | samtools sort -o bams/${sample}.bam
    done 

While you wait, read up on DeNovo stacks! It is very easy to run, a single line, but has many caveats: https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12775 

    denovo_map.pl -T 2 -o denovo -m 10 -M 5 -n 5 --popmap popmap_ingroup --samples reads/

Once your BWA run is done, run gstacks on both scales:

    gstacks -I bams -M popmap_full -O stacks_out/full/
    gstacks -I bams -M popmap_ingroup -O stacks_out/ingroup/
