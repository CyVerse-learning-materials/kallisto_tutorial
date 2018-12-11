.. include:: cyverse_rst_defined_substitutions.txt

|CyVerse logo|_

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

*Build Kallisto Index*
~~~~~~~~~~~~~~~~~~~~~~~

We will now index the Arabidopsis transcriptome imported from Ensembl. This
transcriptome can be used multiple times for future Kallisto analyses and only
needs to be made once.

  1.  If necessary, login to the CyVerse |discovery_enviornment|.

  2. Open the |kallisto index|.

  3. Name your analysis, and if desired enter comments. In the App's 'Input' step
     under 'Index name' enter a name for your index. For this tutorial, name
     your index **Arabidopsis_thaliana.TAIR10.36.cdna.all.fa.index**.
  4. For 'Fasta file' browse to the transcriptome imported in the previous section.

  5. If desired adjust the k-mer length (See |Kallisto paper| for recommendations);
     we will use the default.

  6. Click 'Launch Analyses' to start the job. Click on the Analyses button
     to monitor the job and results.


**Output/Results**

.. list-table::
    :header-rows: 1

    * - Output
      - Description
      - Example
    * - Kallisto Index
      - This is the index file Kallisto will map RNA-Seq reads to.
      - |Example Kallisto index|

----

**Description of results and next steps**

With the index, we will now use the Kallisto Quant App to do a pesudoalignment
of the RNA-Seq reads.

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
.. |discovery_enviornment| raw:: html

    <a href="https://de.cyverse.org/de/" target="_blank">Discovery Environment</a>

.. |kallisto index| raw:: html

    <a href="https://de.cyverse.org/de/?type=apps&app-id=ee7ce21e-645c-11e5-a295-d7e4d1e00ab5&system-id=de" target="_blank">Kallisto-0.42.3-index App</a>

.. |Kallisto paper| raw:: html

    <a href="https://www.nature.com/nbt/journal/v34/n5/full/nbt.3519.html" target="_blank">Kallisto paper</a>

.. |Example transcriptome| raw:: html

    <a href="http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/01_input_transcriptome" target="blank">Example transcriptome</a>

.. |Example Kallisto index| raw:: html

    <a href="http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/02_output_kallisto_index/Arabidopsis_thaliana.TAIR10.36.cdna.all.fa.index" target="_blank">Example Kallisto index</a>
