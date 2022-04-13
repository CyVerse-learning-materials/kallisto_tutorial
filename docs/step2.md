\_

![Home\_Icon](./img/homeicon.png){width="25" height="25"}\_ [Learning
Center Home](http://learning.cyverse.org/)

Build Kallisto Transcriptome Index
==================================

**Description:**

As described in the , RNA-Seq reads are efficiently mapped through a
pseudoalignment process against a reference transcriptome index. We will
build the index in this step.

------------------------------------------------------------------------

\*\*Input <Data:**>

*Build Kallisto Index and Quantify Reads*
-----------------------------------------

We will now index the Arabidopsis transcriptome imported from Ensembl.
This transcriptome can be used multiple times for future Kallisto
analyses and only needs to be made once. In this tutorial, we have 36
fastq files (18 pairs), so you will need to add these to the Kallisto
analyses. Kallisto uses a 'hash-based' pseudo alignment to deliver
extremely fast matching of RNA-Seq reads against the transcriptome
index.

> 1.  If necessary, login to the CyVerse .
> 2.  In the App panel, open the **Kallisto v.0.43.1** app or click this
>     link:
> 3.  Name your analysis, and if desired enter comments and click
>     'Next.' In the App's 'Input' section under 'The transcript fasta
>     file supplied (fasta or gzipped)' browse to and select the
>     transcriptome imported in the previous section.
> 4.  Under Paired of single end choose the format used in your
>     sequencing.
>
>     <div class="admonition">
>
>     Sample data
>
>     For the sample data, choose **Paired**
>
>     </div>
>
>     <div class="admonition note">
>
>     For single-end data you will also need to choose fragment length
>     and fragment standard deviation values in the apps "Options"
>     section. You may also adjust settings for strand-specific reads.
>
>     </div>
>
> 5.  Under 'FASTQ Files (Read 1)' navigate to your data and select all
>     the left-read files (usually R1). For paired-end data also enter
>     the right-read files (usually R2) .
>
>     <div class="admonition">
>
>     Sample data
>
>     For the sample data, navigate to
>     **/iplant/home/shared/cyverse\_training/tutorials/kallisto/00\_input\_fastq\_trimmed**
>
>     -   For FASTQ Files (Read 1) choose all 18 files ending labeled R1
>         (e.g. *SRR1761506\_R1\_001.fastq.gz\_fp.trimmed.fastq.gz*)
>     -   For FASTQ Files (Read 2) choose all 18 files ending labeled R2
>         (e.g. *SRR1761506\_R2\_001.fastq.gz\_fp.trimmed.fastq.gz*)
>
>     </div>
>
> 6.  If desired adjust the bootstrap value (See for recommendations);
>     Click 'Next' to continue.
>
>     <div class="admonition">
>
>     Sample data
>
>     We will use 25.
>
>     </div>
>
> 7.  If desired adjust the resources required and/or click 'Next.'
>
>     <div class="admonition">
>
>     Sample data
>
>     For the sample data, we will not specify resources.
>
>     </div>
>
> 8.  Finally click 'Launch Analyses' to start the job. Click on the
>     Analyses menu to monitor the job and results.

**Output/Results**

Kallisto jobs will generate and index file and 3 output files per read /
read-pair:

------------------------------------------------------------------------

**Description of results and next steps**

First, this application runs the 'kallisto index' command to build the
the index of the transcriptome. Then the 'kallisto quant' command is run
to do the pesudoalignment of the RNA-Seq reads. Kallisto quantifies
RNA-Seq reads against an indexed transcriptome and generates a folder of
results for each set of RNA-Seq reads. Sleuth will be used to examine
the Kallisto results in R Studio.

------------------------------------------------------------------------

**Fix or improve this documentation**

-   Search for an answer:
-   Ask us for help: click on the lower right-hand side of the page
-   Report an issue or submit a change:
-   Send feedback: [<Tutorials@CyVerse.org>](Tutorials@CyVerse.org)

------------------------------------------------------------------------

![Home\_Icon](./img/homeicon.png){width="25" height="25"}\_ [Learning
Center Home](http://learning.cyverse.org/)
