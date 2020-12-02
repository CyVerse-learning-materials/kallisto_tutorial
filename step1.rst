.. include:: cyverse_rst_defined_substitutions.txt
.. include:: custom_urls.txt

|CyVerse_logo|_

|Home_Icon|_
`Learning Center Home <http://learning.cyverse.org/>`_


Organize Kallisto Input Data
-----------------------------

**Description:**

Kallisto has relatively few input requirements. You will need to have either
single or paired end reads, as well as a reference transcriptome. It is
suggested that your RNA-Seq reads are analyzed using |FastQC|, followed by any
additional trimming and filtering using and application such as |Trimmomatic|.

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

*Importing Reference  Transcriptome*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example, we will import a reference transcriptome for Arabidopsis from
Ensembl. In many cases you can find an appropriate transcriptome from Ensembl
for your organism of interest, or provide your own fasta-formmatted transcriptome.

  1. Go to the Ensembl Plants Arabidopsis page: |ensembl|.

  2. In the 'Gene annotation' section, click on the 'Download genes, cDNAs,
     ncRNA, proteins' 'FASTA' link: |fasta link|.

  3. The transcriptome files will be located in the 'cdna' folder: |ftp link|.

  4. Right-click/command-click on the
     `Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz` link and copy the link location/URL (usually right-click) to the clipboard

  5. Login to the CyVerse |discovery_enviornment|.

  6. Open the Data window. In your home directory, create a folder to organize
     your Kallisto project

  7. In the created folder, go to the Data window's 'Upload' menu, and select,
     'Import from URL...'; paste in the Ensembl link and click 'Import from URL' to begin the import

**Output/Results**

.. list-table::
    :header-rows: 1

    * - Output
      - Description
      - Example
    * - Reference transcriptome
      - fasta
      - |Example transcriptome|


----

**Description of results and next steps**

This example transcriptome will be indexed in the next step. RNA-Seq reads will
be mapped against this set of transcripts. Once you have the transcriptome and
your RNA-Seq reads, you can proceed with the next step. We suggest organizing
your RNA-Seq reads in the folder created in step 6 above.

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
