|CyVerse logo|_

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
      - `Example Kallisto index <http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/02_output_kallisto_index/Arabidopsis_thaliana.TAIR10.36.cdna.all.fa.index>`_
    * - RNA-Seq Reads
      - Cleaned fastq files
      - `Example fastq files <http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/00_input_fastq_trimmed>`_



*Quantify RNA-Seq reads with Kallisto*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A Kallisto analyses must be run for each mapping of RNA-Seq reads to the index.
In this tutorial, we have 36 fastq files (18 pairs), so you will need to launch
18 Kallisto analyses. It is sufficient here to launch a single Kallisto job to
examine the input and then use the completed results (which are small files) for
Sleuth analyses.

  1. If necessary, login to the CyVerse `Discovery Environment <https://de.cyverse.org/de/>`_

  2. Open the `Kallisto-0.42.3-Quant-PE App <https://de.cyverse.org/de/?type=apps&app-id=38159000-83da-11e5-be5b-d7c855bb70b2&system-id=de>`_

  3. Name your analysis, and if desired enter comments. In the App's 'Input' step
     under 'Index file' browse to and select the Kallisto index generated in the previous
     tutorial section. In the output directory, enter the name for the output directory
     that will be created. For this tutorial, name your output directory **pair01_wt_mock_r1**
     (This is for our first pair of WT reads, mock treatment, replicate 1).

  4. Under 'Input Read1&Reaad2 fastq files' enter a single pair (foward/reverse)
     of reads. In this tutorial click 'Add' to select the following files located in
     *Community Data > cyverse_training > tutorials > kallisto > 00_input_fastq_trimmed*:
     - `SRR1761506_R1_001.fastq.gz_fp.trimmed.fastq.gz <http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/00_input_fastq_trimmed/SRR1761506_R1_001.fastq.gz_fp.trimmed.fastq.gz>`_
     - `SRR1761506_R2_001.fastq.gz_fp.trimmed.fastq.gz <http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/00_input_fastq_trimmed/SRR1761506_R2_001.fastq.gz_rp.trimmed.fastq.gz>`_

  5. Click 'Launch Analyses' to launch the job and monitor its progress.


**Output/Results**

Each Kallisto job will generate 3 files


.. list-table::
    :header-rows: 1

    * - Output
      - Description
      - Example
    * - abundances.h5
      - HDF5 binary file containing run info, abundance estimates,
        bootstrap estimates, and transcript length information length.
        This file can be read in by Sleuth
      - `example abundance.h5 <http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/03_output_kallisto_results/pair01_wt_mock_r1/abundance.h5>`_
    * - abundances.tsv
      - plaintext file of the abundance estimates. It does not contains
        bootstrap estimates. When plaintext mode is selected; output plaintext
        abundance estimates. Alternatively, kallisto h5dump will output
        an HDF5 file to plaintext. The first line contains a header for each
        column, including estimated counts, TPM, effective length.
      - `example abundance.tsv <http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/03_output_kallisto_results/pair01_wt_mock_r1/abundance.tsv>`_
    * - run_info.json
      - a json file containing information about the run
      - `example json <http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/03_output_kallisto_results/pair01_wt_mock_r1/run_info.json>`_

----

**Summary**
~~~~~~~~~~~

Kallisto quantifies RNA-Seq reads against an indexed transcriptome and generates
a folder of results for each set of RNA-Seq reads.

----

**Next Steps:**

Sleuth will be used to examine the Kallisto results in R Studio.

----

More help and additional information
`````````````````````````````````````

..
    Short description and links to any reading materials (KEEP THIS on LAST PAGE
    of Tutorial)

Search for an answer:
    `CyVerse Learning Center <http://learning.cyverse.org>`_ or
    `CyVerse Wiki <https://wiki.cyverse.org>`_

Post your question to the user forum:
    `Ask CyVerse <http://ask.iplantcollaborative.org/questions>`_

----

**Fix or improve this documentation**

- On Github: `Repo link <https://github.com/CyVerse-learning-materials/kallisto_tutorial>`_
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
