.. include:: cyverse_rst_defined_substitutions.txt

|CyVerse_logo|_

|Home_Icon|_
`Learning Center Home <http://learning.cyverse.org/>`_


Quantify Reads with Kallisto
----------------------------

**Description:**

Kallisto uses a 'hash-based' pseudo alignment to deliver extremely fast matching
of RNA-Seq reads against the transcriptome index. Each Kallisto job in this
tutorial will take only several minutes to complete.

----

**Input Data:**

.. list-table::
    :header-rows: 1

    * - Input
      - Description
      - Example
    * - Kallisto Index
      - Indexed transcriptome
      - |Example Kallisto index|
    * - RNA-Seq Reads
      - Cleaned fastq files
      - |Example fastq files|



*Quantify RNA-Seq reads with Kallisto*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A Kallisto analyses must be run for each mapping of RNA-Seq reads to the index.
In this tutorial, we have 36 fastq files (18 pairs), so you will need to add these
to the Kallisto analyses. It is sufficient here to launch a single Kallisto job to
examine the input and then use the completed results (which are small files) for
Sleuth analyses.

  1. If necessary, login to the CyVerse |discovery_enviornment|.

  2. Open the |kallisto quant|.

  3. Name your analysis, and if desired enter comments. In the App's 'Input' step
     under 'Index file' browse to and select the Kallisto index generated in the previous
     tutorial section. In the output directory, enter the name for the output directory
     that will be created. For this tutorial, name your output directory **pair01_wt_mock_r1**
     (This is for our first pair of WT reads, mock treatment, replicate 1).

  4. Under 'Read 1 Fastq files' and 'Read 1 Fastq files' the respective right ang left sequences.
     In this tutorial click 'Add' to select the following files located in
     *Community Data > cyverse_training > tutorials > kallisto > 00_input_fastq_trimmed*:
     - |file1|
     - |file2|

  5. Click 'Launch Analyses' to launch the job and monitor its progress.


**Output/Results**

Kallisto jobs will generate 3 files per read pair:


.. list-table::
    :header-rows: 1

    * - Output
      - Description
      - Example
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

Kallisto quantifies RNA-Seq reads against an indexed transcriptome and generates
a folder of results for each set of RNA-Seq reads. Sleuth will be used to examine
the Kallisto results in R Studio.

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
.. |discovery_enviornment| raw:: html

    <a href="https://de.cyverse.org/de/" target="_blank">Discovery Environment</a>

.. |kallisto quant| raw:: html

    <a href="https://de.cyverse.org/de/?type=apps&app-id=c304d9de-66eb-11e5-83d0-b36f5d747f5c&system-id=de" target="_blank">Kallisto-0.42.3-quant-PE App</a>

.. |file1| raw:: html

    <a href="http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/00_input_fastq_trimmed/SRR1761506_R1_001.fastq.gz_fp.trimmed.fastq.gz" target="_blank">SRR1761506_R1_001.fastq.gz_fp.trimmed.fastq.gz</a>

.. |file2| raw:: html

    <a href="http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/00_input_fastq_trimmed/SRR1761506_R2_001.fastq.gz_rp.trimmed.fastq.gz" target="_blank">SRR1761506_R2_001.fastq.gz_fp.trimmed.fastq.gz</a>

.. |Example Kallisto index| raw:: html

    <a href="http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/02_output_kallisto_index/Arabidopsis_thaliana.TAIR10.36.cdna.all.fa.index" target="_blank">Example Kallisto index</a>

.. |Example fastq files| raw:: html

    <a href="http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/00_input_fastq_trimmed" target="_blank">Example fastq files</a>

.. |example abundance.h5| raw:: html

    <a href="http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/03_output_kallisto_results/pair01_wt_mock_r1/abundance.h5" target="_blank">example abundance.h5</a>

.. |example abundance.tsv| raw:: html

    <a href="http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/03_output_kallisto_results/pair01_wt_mock_r1/abundance.tsv" target="_blank">example abundance.tsv</a>

.. |example json| raw:: html

    <a href="http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/03_output_kallisto_results/pair01_wt_mock_r1/run_info.json" target="_blank">example json</a>
