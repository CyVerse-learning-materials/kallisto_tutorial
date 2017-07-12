|CyVerse logo|_

|Home_Icon|_
`Learning Center Home <http://learning.cyverse.org/>`_

RNA-Seq with Kallisto and Sleuth
=================================

Goal
----

**Analyze RNA-Seq data for differential expression**. `Kallisto <https://pachterlab.github.io/kallisto/about>`_ is a quick, highly-efficient
software for quantifying transcript abundances in an RNA-Seq experiment. Even on
a typical laptop, Kallisto can quantify 30 million reads in less than 3 minutes.
Integrated into CyVerse, you can take advantage of CyVerse data management tools
to process your reads, do the Kallisto quantification, and analyze your reads
with the Kallisto companion software `Sleuth <https://pachterlab.github.io/sleuth/about>`_ in an R-Studio environment.


----

.. toctree::
	:maxdepth: 2

	Tutorial home <self>
	Organize Kallisto Input Data <step1.rst>
	Build Kallisto Transcriptome Index <step2.rst>
	Quantify Reads with Kallisto <step3.rst>
	Analyze Kallisto Results with Sleuth <step4.rst>

..
	#### Comment:This tutorial can have multiple pages. The table of contents assumes
	you have an additional page called 'Step Owo' with content located in 'step1.rst'
	Edit these titles and filenames as needed ####

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
      - `Register <https://user.cyverse.org/>`_
    * - Atmosphere access (optional)
      - This tutorial will use R studio in Atmosphere; if desired you can
        complete these sections by installing the Sleuth tools on your own R
        instance
      - `Request Access <https://user.cyverse.org/>`_

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
      - `Data Store <http://www.cyverse.org/data-store>`_
      - `Data Store Manual <https://wiki.cyverse.org/wiki/display/DS/Data+Store+Table+of+Contents>`_
      - `Guide <http://learning.cyverse.org/projects/cyverse-discovery-environment-guide/>`__
    * - Discovery Environment
      - Web/Point-and-click
      - `Discovery Environment <https://de.cyverse.org/de/>`_
      - `DE Manual <https://wiki.cyverse.org/wiki/display/DEmanual/Table+of+Contents>`_
      - `Guide <http://learning.cyverse.org/projects/cyverse-discovery-environment-guide/>`__
    * - Atmosphere
      - Command line (ssh) and/or Desktop (VNC)
      - `Atmosphere <https://atmo.cyverse.org>`_
      - `Atmosphere Manual <https://wiki.cyverse.org/wiki/display/atmman/Atmosphere+Manual+Table+of+Contents>`_
      - `Guide <https://cyverse-atmosphere-guide.readthedocs-hosted.com/en/latest/>`__


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
    * - Kallisto-0.42.3-INDEX
      - 0.42.3
      - Kallisto Index Builder
      -	`DE App link <https://de.cyverse.org/de/?type=apps&app-id=ffd24602-923e-11e5-843a-e7021d2c7752&system-id=de>`__
      - `Original documentation <https://pachterlab.github.io/kallisto/manual>`_
    * - Kallisto-0.42.3-Quant-PE
      - 0.43.3
      - Kallisto Quantification
      - `DE App link <https://de.cyverse.org/de/?type=apps&app-id=38159000-83da-11e5-be5b-d7c855bb70b2&system-id=de>`__
      - `Original documentation`_


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
      -	`Image <https://atmo.cyverse.org/application/images/1479>`_
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
        such as `Trimmomatic <https://cyverse-trimmomatic-quickstart.readthedocs-hosted.com/en/latest/>`_
      - `Example fastq files <http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/00_input_fastq_trimmed>`_
    * - Refference transcriptome
      - fasta
      - Transcriptome for your organism of interest
      - `Example transcriptome <http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/01_input_transcriptome>`_


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
