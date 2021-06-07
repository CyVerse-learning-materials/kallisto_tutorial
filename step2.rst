.. include:: cyverse_rst_defined_substitutions.txt
.. include:: custom_urls.txt

|CyVerse_logo|_

|Home_Icon|_
`Learning Center Home <http://learning.cyverse.org/>`_


Build Kallisto Transcriptome Index
-----------------------------------

**Description:**

As described in the |Kallisto paper|,
RNA-Seq reads are efficiently mapped through a pseudoalignment process against a
reference transcriptome index. We will build the index in this step.

..
	#### Comment: short text description goes here ####

----

**Input Data:**

.. list-table::
    :header-rows: 1

    * - Input
      - Description
      - Example
    * - Reference transcriptome
      - fasta
      - |Example transcriptome|
    * - RNA-Seq Reads
      - Cleaned fastq files
      - |Example fastq files|


*Build Kallisto Index and Quantify Reads*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We will now index the Arabidopsis transcriptome imported from Ensembl. This
transcriptome can be used multiple times for future Kallisto analyses and only
needs to be made once. In this tutorial, we have 36 fastq files (18 pairs), so
you will need to add  these to the Kallisto analyses. Kallisto uses a 'hash-based' pseudo alignment to deliver extremely fast matching
of RNA-Seq reads against the transcriptome index.

  1. If necessary, login to the CyVerse |discovery_enviornment|.

  2. In the App panel, open the **Kallisto v.0.43.1** app or click this
     link:
     |kallisto app|

  3. Name your analysis, and if desired enter comments and click
     'Next.' In the App's 'Input' section under 'The transcript fasta
     file supplied (fasta or gzipped)' browse to and select the
     transcriptome imported in the previous section.

  4. Under Paired of single end choose the format used in your
     sequencing.

     .. admonition:: Sample data

       For the sample data, choose **Paired**

     .. note::

        For single-end data you will also need to choose fragment
        length and
        fragment standard deviation values in the apps "Options"
        section.
        You may also adjust settings for strand-specific reads.

  5. Under 'FASTQ Files (Read 1)' navigate to your data and select all
     the left-read files (usually R1). For paired-end data also enter
     the right-read files (usually R2) .

     .. admonition:: Sample data

       For the sample data, navigate to **/iplant/home/shared/cyverse_training/tutorials/kallisto/00_input_fastq_trimmed**

       - For FASTQ Files (Read 1) choose all 18 files ending labeled R1 (e.g.
         *SRR1761506_R1_001.fastq.gz_fp.trimmed.fastq.gz*)
       - For FASTQ Files (Read 2) choose all 18 files ending labeled R2 (e.g.
         *SRR1761506_R2_001.fastq.gz_fp.trimmed.fastq.gz*)

  6. If desired adjust the bootstrap value (See |Kallisto paper| for
     recommendations); Click 'Next' to continue.

     .. admonition:: Sample data

        We will use 25.

  7. If desired adjust the resources required and/or click 'Next.'

     .. admonition:: Sample data

       For the sample data, we will not specify resources.


  8. Finally click 'Launch Analyses' to start the job. Click on the
     Analyses menu to monitor the job and results.


**Output/Results**

Kallisto jobs will generate and index file and 3 output files per read
/ read-pair:

.. list-table::
    :header-rows: 1

    * - Output
      - Description
      - Example
    * - Kallisto Index
      - This is the index file Kallisto will map RNA-Seq reads to.
      - |Example Kallisto index|
    * - abundances.h5
      - HDF5 binary file containing run info, abundance estimates,
        bootstrap estimates, and transcript length information length.
        This file can be read in by Sleuth
      - |example abundance.h5|
    * - abundances.tsv
      - plaintext file of the abundance estimates. It does not contains
        bootstrap estimates. When plaintext mode is selected; output plaintext
        abundance estimates. Alternatively, kallisto h5dump will output
        an HDF5 file to plaintext. The first line contains a header for each
        column, including estimated counts, TPM, effective length.
      - |example abundance.tsv|
    * - run_info.json
      - a json file containing information about the run
      - |example json|


----

**Description of results and next steps**

First, this application runs the 'kallisto index' command to build the the
index of the transcriptome. Then the 'kallisto quant' command is run to do the
pesudoalignment of the RNA-Seq reads. Kallisto quantifies RNA-Seq reads against
an indexed transcriptome and generates a folder of results for each set of RNA-
Seq reads. Sleuth will be used to examine the Kallisto results in R Studio.

----

**Fix or improve this documentation**

- Search for an answer:
  |CyVerse Learning Center|
- Ask us for help:
  click |Intercom| on the lower right-hand side of the page
- Report an issue or submit a change:
  |Github Repo Link|
- Send feedback: `Tutorials@CyVerse.org <Tutorials@CyVerse.org>`_


----

|Home_Icon|_
`Learning Center Home <http://learning.cyverse.org/>`_

.. |CyVerse logo| image:: ./img/cyverse_rgb.png
    :width: 500
    :height: 100
.. _CyVerse logo: http://learning.cyverse.org/
.. |Home_Icon| image:: ./img/homeicon.png
    :width: 25
    :height: 25
.. _Home_Icon: http://learning.cyverse.org/
