# RNA-Seq with Kallisto and Sleuth

## Goal

**Analyze RNA-Seq data for differential expression**. [Kallisto](https://pachterlab.github.io/kallisto/manual) is a quick, highly-efficient software for quantifying transcript abundances in an
RNA-Seq experiment. Even on a typical laptop, Kallisto can quantify 30
million reads in less than 3 minutes. Integrated into CyVerse, you can
take advantage of CyVerse data management tools to process your reads,
do the Kallisto quantification, and analyze your reads with the Kallisto
companion software [Sleuth](https://pachterlab.github.io/sleuth/about) in an R-Studio environment.

------------------------------------------------------------------------

## Manual Maintainer(s)

Who to contact if this manual needs fixing. You can also email
[<Tutorials@CyVerse.org>](Tutorials@CyVerse.org)

| Maintainer | Institution | Contact |
|---|---|---|
| Jason Williams | CyVerse / Cold Spring Harbor Laboratory | [Williams@cshl.edu](Williams@cshl.edu) |

------------------------------------------------------------------------

## Prerequisites

### Downloads, access, and services

*In order to complete this tutorial you will need access to the
following services/software*

| Prerequisite | Preparation/Notes | Link/Download |
|---|---|---|
| CyVerse account | You will need a CyVerse account to complete this exercise | [CyVerse User Portal](https://user.cyverse.org/) |

### Platform(s)

*We will use the following CyVerse platform(s):*

| Platform	| Interface	| Link |
|---|---|---|
| Data Store | GUI/Command line	| [Data Store](https://data.cyverse.org/) |
| Discovery Environment	| Web/Point-and-click | [Discovery Environment](https://de.cyverse.org/) |

### Application(s) used

**Discovery Environment App(s):**

| App name | Version | Description | App link | Notes/other links |
|---|---|---|---|---|
| Kallisto-v.0.43.1 | 0.43.1 | Kallisto v.0.43.1 | [Kallisto app](https://de.cyverse.org/de/?type=quick-launch&quick-launch-id=6132e25c-6576-4c84-bd6f-9e343e5ef03a&app-id=c341ba8c-30ad-11e8-8fb4-008cfa5ae621) | [Kallisto manual](https://pachterlab.github.io/kallisto/manual) |
| RStudio | Sleuth 0.30.0 | RStudio with Sleuth (v.0.30.0) and dependencies | [Sleuth app](https://de.cyverse.org/de/?type=quick-launch&quick-launch-id=3125ee9a-9f0c-4f4c-8efd-aa6f7ea00405&app-id=8eb1291c-34ea-11eb-b90c-008cfa5ae621) | [Sleuth manual](https://pachterlab.github.io/sleuth/about) |

### Input and example data

*In order to complete this tutorial you will need to have the following
inputs prepared*

| Input File(s)	| Format | Preparation/Notes | Example Data |
|---|---|---|---|
| RNA-Seq reads	| FastQ (may also be compressed, e.g. fastq.gz) | These reads should have been cleaned by upstream tools such as [Trimmomatic](https://cyverse-trimmomatic-quickstart.readthedocs-hosted.com/en/latest/)	| [Example FastQ files](http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/00_input_fastq_trimmed) |
| Reference transcriptome | fasta | Transcriptome for your organism of interest	| [Example transcriptome](http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/01_input_transcriptome) |

## Sample Data and Working with Your Own Data

!!! Note "Sample data"
        **About the Sample Dataset** In this tutorial, we are using publicly
        available data from the SRA. This tutorial will start with cleaned and
        processed reads. The SRA experiment used data from bioproject [PRJNA272719](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA272719). The abstract from that project is reprinted here:

        'To survey transcriptome changes by the mutations of a DNA demethylase
        ROS1 responding to a phytohormone abscisic acid, we performed the
        Next-gen sequencing (NGS) associated RNA-seq analysis. Two ROS1 knockout
        lines (ros1-3, ros1-4; Penterman et al. 2007 [PMID: 17409185]) with
        the wild-type Col line (wt) were subjected. Overall design: Three
        samples (ros1-3, ros1-4 and wt), biological triplicates, ABA or mock
        treatment, using Illumina HiSeq 2500 system' [citation](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA272719).

!!! Tip
        **Working with your own data**

        If you have your own FASTQ files upload them to CyVerse using
        instructions in the CyVerse (e.g. iCommands/Cyberduck).
