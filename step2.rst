|CyVerse logo|_

|Home_Icon|_
`Learning Center Home <http://learning.cyverse.org/>`_


Build Kallisto Transcriptome Index
-----------------------------------

**Description:**

As described in the `Kallisto paper <https://www.nature.com/nbt/journal/v34/n5/full/nbt.3519.html>`_,
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
      - `Example transcriptome <http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/01_input_transcriptome>`_

*Build Kallisto Index*
~~~~~~~~~~~~~~~~~~~~~~~

We will now index the Arabidopsis transcriptome imported from Ensembl. This
transcriptome can be used multiple times for future Kallisto analyses and only
needs to be made once.

  1.  If necessary, login to the CyVerse `Discovery Environment <https://de.cyverse.org/de/>`_

  2. Open the `Kallisto-0.42.3-INDEX App <https://de.cyverse.org/de/?type=apps&app-id=ffd24602-923e-11e5-843a-e7021d2c7752&system-id=de>`_

  3. Name your analysis, and if desired enter comments. In the App's 'Input' step
     under 'Index name' enter a name for your index. For this tutorial, name
     your index **Arabidopsis_thaliana.TAIR10.36.cdna.all.fa.index**.

  4. If desired adjust the k-mer length (See `Kallisto paper`_ for recommendations);
     we will use the default.

  5. Click 'Launch Analyses' to start the job. Click on the Analyses button
     to monitor the job and results.


**Output/Results**

.. list-table::
    :header-rows: 1

    * - Output
      - Description
      - Example
    * - Kallisto Index
      - This is the index file Kallisto will map RNA-Seq reads to.
      - `Example Kallisto index <http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/02_output_kallisto_index/Arabidopsis_thaliana.TAIR10.36.cdna.all.fa.index>`_

----

**Next Steps:**

With the index, we will now use the Kallisto Quant App to do a pesudoalignment
of the RNA-Seq reads.

----

More help and additional information
------------------------------------

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
