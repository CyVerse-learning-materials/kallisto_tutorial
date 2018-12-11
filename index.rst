.. include:: cyverse_rst_defined_substitutions.txt

|CyVerse logo|_

|Home_Icon|_
`Learning Center Home <http://learning.cyverse.org/>`_

RNA-Seq with Kallisto and Sleuth
=================================

Goal
----

**Analyze RNA-Seq data for differential expression**. |Kallisto manual| is a quick, highly-efficient
software for quantifying transcript abundances in an RNA-Seq experiment. Even on
a typical laptop, Kallisto can quantify 30 million reads in less than 3 minutes.
Integrated into CyVerse, you can take advantage of CyVerse data management tools
to process your reads, do the Kallisto quantification, and analyze your reads
with the Kallisto companion software |Sleuth manual| in an R-Studio environment.


----

.. toctree::
	:maxdepth: 2

	Tutorial home <self>
	Organize Kallisto Input Data <step1.rst>
	Build Kallisto Transcriptome Index <step2.rst>
	Quantify Reads with Kallisto <step3.rst>
	Analyze Kallisto Results with Sleuth <step4.rst>


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
    * - Atmosphere access (optional)
      - This tutorial will use R studio in Atmosphere; if desired you can
        complete these sections by installing the Sleuth tools on your own R
        instance
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
    * - Atmosphere
      - Command line (ssh) and/or Desktop (VNC)
      - |Atmosphere|
      - |Atmosphere Manual|
      - |Atmosphere Guide|

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
    * - Kallisto-0.42.3-index
      - 0.42.3
      - Kallisto-0.42.3-index
      - |Kallisto Index|
      - |Kallisto manual|
    * - Kallisto-0.42.3-quant-PE
      - 0.42.3
      - Kallisto Quantification
      - |Kallisto Quant App|
      - |Kallisto manual|


**Atmosphere Image(s):**

.. list-table::
    :header-rows: 1

    * - Image name
      - Version
      - Description
      - Link
      - Notes/other links
    * - CyVerse Training Workshop
      - 1.1.2
      - Image for use at CyVerse Training Workshops
      -	|CyVerse Training Workshop Image|
      - This image has Kallisto and Sleuth installed. Once started, an
        R Studio Server will be available at your image ip address, port 8787
        (e.g. image.ip.address:8787)



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
      - Fastq (may also be compressed, e.g. fastq.gz)
      - These reads should have been cleaned by upstream tools
        such as |Trimmomatic|
      - |Example fastQ files|
    * - Reference transcriptome
      - fasta
      - Transcriptome for your organism of interest
      - |Example transcriptome|


----

**Fix or improve this documentation**

Search for an answer:
|CyVerse Learning Center| or
|CyVerse Wiki|

Post your question to the user forum:
|Ask CyVerse|

----

|Home_Icon|_
`Learning Center Home <http://learning.cyverse.org/>`__

.. |Kallisto Index| raw:: html

   <a href="https://de.cyverse.org/de/?type=apps&app-id=ee7ce21e-645c-11e5-a295-d7e4d1e00ab5&system-id=de" target="blank">Kallisto-0.42.3-index</a>

.. |Kallisto manual| raw:: html

   <a href="https://pachterlab.github.io/kallisto/manual" target="blank">Kallisto manual</a>

.. |Kallisto Quant App| raw:: html

   <a href="https://de.cyverse.org/de/?type=apps&app-id=c304d9de-66eb-11e5-83d0-b36f5d747f5c&system-id=de" target="blank">Kallisto-0.42.3-quant-PEp</a>

.. |CyVerse Training Workshop Image| raw:: html

   <a href="https://atmo.cyverse.org/application/images/1479" target="blank">CyVerse Training Workshop Image</a>

.. |Sleuth manual| raw:: html

   <a href="https://pachterlab.github.io/sleuth/about" target="blank">Sleuth manual</a>

.. |Trimmomatic| raw:: html

    <a href="https://cyverse-trimmomatic-quickstart.readthedocs-hosted.com/en/latest/" target="blank">Trimmomatic</a>

.. |Example fastQ files| raw:: html

    <a href="http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/00_input_fastq_trimmed" target="blank">Example fastQ files</a>

.. |Example transcriptome| raw:: html

    <a href="http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/01_input_transcriptome" target="blank">Example transcriptome</a>

.. || raw:: html

    <a href="" target="blank"></a>
