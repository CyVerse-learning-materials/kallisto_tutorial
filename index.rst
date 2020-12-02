.. include:: cyverse_rst_defined_substitutions.txt
.. include:: custom_urls.txt

|CyVerse_logo|_

|Home_Icon|_
`Learning Center Home <http://learning.cyverse.org/>`_


RNA-Seq with Kallisto and Sleuth
=================================

Goal
----

**Analyze RNA-Seq data for differential expression**. |Kallisto manual| is a
quick, highly-efficient software for quantifying transcript abundances in an
RNA-Seq experiment. Even on a typical laptop, Kallisto can quantify 30 million
reads in less than 3 minutes. Integrated into CyVerse, you can take advantage
of CyVerse data management tools to process your reads, do the Kallisto
quantification, and analyze your reads with the Kallisto companion software
|Sleuth manual| in an R-Studio environment.


----

Manual Maintainer(s)
------------------------

Who to contact if this manual needs fixing. You can also email
`Tutorials@CyVerse.org <Tutorials@CyVerse.org>`_

.. list-table::
    :header-rows: 1

    * - Maintainer
      - Institution
      - Contact
    * - Jason Williams
      - CyVerse / Cold Spring Harbor Laboratory
      - Williams@cshl.edu


----

.. toctree::
	:maxdepth: 2

	Tutorial home <self>
	Organize Kallisto Input Data <step1.rst>
	Build  Transcriptome Index and Quantify Reads with Kallisto <step2.rst>
	Analyze Kallisto Results with Sleuth <step3.rst>


Prerequisites
-------------

Downloads, access, and services
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*In order to complete this tutorial you will need access to the following services/software*

..
	#### comment: delete any row not needed in this table ####

.. list-table::
    :header-rows: 1

    * - Prerequisite
      - Preparation/Notes
      - Link/Download
    * - CyVerse account
      - You will need a CyVerse account to complete this exercise
      - |CyVerse User Portal|


Platform(s)
~~~~~~~~~~~

*We will use the following CyVerse platform(s):*

..
	#### comment: delete any row not needed in this table ####

.. list-table::
    :header-rows: 1

    * - Platform
      - Interface
      - Link
      - Platform Documentation
      - Quick Start
    * - Data Store
      - GUI/Command line
      - |Data Store|
      - |Data Store Manual|
      - |Data Store Guide|
    * - Discovery Environment
      - Web/Point-and-click
      - |Discovery Environment|
      - |DE Manual|
      - |Discovery Environment Guide|

Application(s) used
~~~~~~~~~~~~~~~~~~~
..
	#### Comment: these tables are examples, delete whatever is unnecessary ####

**Discovery Environment App(s):**

.. list-table::
    :header-rows: 1

    * - App name
      - Version
      - Description
      - App link
      - Notes/other links
    * - Kallisto-v.0.43.1
      - 0.43.1
      - Kallisto v.0.43.1
      - |Kallisto app|
      - |Kallisto manual|
    * - RStudio Sleuth
      - 0.30.0
      - RStudio with Sleuth (v.0.30.0) and dependencies
      - |Sleuth app|
      - |Sleuth manual|


Input and example data
~~~~~~~~~~~~~~~~~~~~~~

*In order to complete this tutorial you will need to have the following inputs prepared*

..
	#### comment: delete any row not needed in this table ####

.. list-table::
    :header-rows: 1

    * - Input File(s)
      - Format
      - Preparation/Notes
      - Example Data
    * - RNA-Seq reads
      - FastQ (may also be compressed, e.g. fastq.gz)
      - These reads should have been cleaned by upstream tools
        such as |Trimmomatic|
      - |Example FastQ files|
    * - Reference transcriptome
      - fasta
      - Transcriptome for your organism of interest
      - |Example transcriptome|


Sample Data and Working with Your Own Data
-----------------------------------------------

.. admonition:: Sample data

    **About the Sample Dataset**
    In this tutorial, we are using publicly available data from the SRA. This
    tutorial will start with cleaned and processed reads. The SRA experiment
    used data from bioproject |PRJNA272719|. The abstract from that project is
    reprinted here:

    'To survey transcriptome changes by the mutations of a DNA demethylase ROS1
    responding to a phytohormone abscisic acid, we performed the Next-gen
    sequencing (NGS) associated RNA-seq analysis. Two ROS1 knockout lines
    (ros1-3, ros1-4; Penterman et al. 2007 [PMID: 17409185])
    with the wild-type Col line (wt) were subjected. Overall design:
    Three samples (ros1-3, ros1-4 and wt), biological triplicates, ABA or mock
    treatment, using Illumina HiSeq 2500 system' |citation|.

.. tip::

    **Working with your own data**

    If you have your own FASTQ files upload them to CyVerse using instructions
    in the CyVerse |Data Store Guide| (e.g. iCommands/Cyberduck).


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
`Learning Center Home <http://learning.cyverse.org/>`__
