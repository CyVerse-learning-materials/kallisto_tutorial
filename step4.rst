|CyVerse logo|_

|Home_Icon|_
`Learning Center Home <http://learning.cyverse.org/>`_


Analyze Kallisto Results with Sleuth
--------------------------------------

**Description:**

`Sleuth <https://pachterlab.github.io/sleuth/about>`_ is a program for analysis
of RNA-Seq experiments for which transcript abundances have been quantified with
kallisto. In this tutorial, we will use R Studio being served from an Atmosphere
instance.

----

**Input Data:**

.. list-table::
    :header-rows: 1

    * - Input
      - Description
      - Example
    * - Kallisto results folder(s)
      - Outputs from a Kallisto quantification
      - `Example directory of Kallisto results <http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/03_output_kallisto_results>`_


*Import Kallisto results into R and analyze with Sleuth*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. tip::

    If you are not familiar with Atmosphere and transfering data into an instance
    see the `Atmosphere guide <https://cyverse-atmosphere-guide.readthedocs-hosted.com/en/latest/>`_

We will import the Kallisto results into an RStudio session being run from
an Atmosphere image. Then we will follow a R script based on the  `Sleuth Walkthoughs <https://pachterlab.github.io/sleuth/walkthroughs>`_

  1. If necessary, launch the `CyVerse Workshop Training Image <https://atmo.cyverse.org/application/images/1479>`_;
     a 'Small1' instance size should be sufficient.

  2. When the instance becomes available, open a web browser and navigate to your
     R studio session. This will be located at URL based on your image ip:
     "http ://your.image.ip:8787";

  3. Login to Atmosphere via Webshell or SSH; initialize iCommands, and use the
     following commands to import the test data and R scripts to your home directory:

     .. code:: bash

        iget -rPV /iplant/home/shared/cyverse_training/tutorials/kallisto/03_output_kallisto_results
        iget -rPV /iplant/home/shared/cyverse_training/tutorials/kallisto/04_sleuth_R

  4. Load the `sample_kallisto_script.R <http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/04_sleuth_R/sample_kallisto_script.R>`_
     into your R Studio session. Use the `Kallisto_demo_tsv <http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/04_sleuth_R/kallisto_demo.tsv>`_ file for the step
     where you describe experimental design metadata.


.. code:: R

     # This tutorial is based on the one by the Pachterlab here:
     # https://pachterlab.github.io/sleuth_walkthroughs/boj/analysis.html

     # Load sleuth
     suppressMessages({library("sleuth")})

     # modify the following line so the location of your kallisto results (each
     # in their own folder) are specified. Hint: don't use slashes to
     # separate your directory names, use a comma-separated list of
     # quoted directory names.
     sample_id <- dir(file.path("./03_output_kallisto_results/"))

     # Running the next line should give you the names of your samples
     # assuming your kallisto result folders have been appropriately named
     sample_id

     #specify the full path relative to this script
     kal_dirs <- file.path("./03_output_kallisto_results",sample_id)

     # should show full path to results containing the kallisto outputs
     kal_dirs

     # create a table for the experimental design
     s2c <- read.table(file.path("kallisto_demo.tsv"),
                       header = TRUE,
                       stringsAsFactors = FALSE,
                       sep = "\t")

     # Check the table
     s2c

     # Add the directories as a final column called path
     s2c <- dplyr::mutate(s2c, path = kal_dirs)

     # check the table

     print(s2c)


     # install cowplot and load for better plots
     install.packages("cowplot")
     library(cowplot)


     # prepare to associate transcripts to genes - you will have to
     # do a bit of checking to ensure you have the proper names
     # and attributes for your genome of interest
     # this step is not essential
     # see some instructions on getting your gene names here:
     # http://www.ensembl.info/blog/2015/06/01/biomart-or-how-to-access-the-ensembl-data-from-r/

     #load biomaRt
     library(biomaRt)

     #load arabidopsis genes
     mart <- biomaRt::useMart(biomart = "plants_mart",
                              dataset = "athaliana_eg_gene",
                              host = "plants.ensembl.org")
     #get target transcripts

     ttg <- biomaRt::getBM(
       attributes = c("ensembl_transcript_id", "transcript_version",
                      "ensembl_gene_id", "external_gene_name", "description",
                      "transcript_biotype"),
                       mart = mart)

     # do some renaming
     ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,
                          ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

     # Check your listing of transcripts
     head(ttg)

     # Create the "Sleuth Object" a data structure that contains our results

     so <- sleuth_prep(s2c,
                       target_mapping = ttg,
                       aggregation_column = 'ens_gene',
                       extra_bootstrap_summary = TRUE)

     # Check the structure of the data with a PCA plot
     # PCA by treatment shows clear difference between ABA and mock
     # good! this was SRA data - so glad this worked!
     # You will get some expected warnings

     plot_pca(so, color_by = 'treatment_s')

     # you can also add labels to the plot
     plot_pca(so, color_by = 'treatment_s',
              text_labels = TRUE)

     # We can also see genes involved in the the 1st PC by looking
     # at the loadings (primary genes whose linear combinations define
     # the principal components)

     plot_loadings(so, pc_input = 1)

     # See more about the gene most influential in this dataset
     # https://www.arabidopsis.org/servlets/TairObject?type=gene&name=At2g34420.1


     # Testing for differential genes
     # Create the full model (i.e. a model that contains all covariates)
     # then create an additional model with respect to one of the covariates
     # finally compare both models to identify the differences

     so <- sleuth_fit(so, ~ genotype_variation_s + treatment_s, 'full')
     so <- sleuth_fit(so, ~genotype_variation_s, 'treatment_s')
     # likelihood ratio test
     so <- sleuth_lrt(so, 'treatment_s', 'full')

     #get the full results
     full_results <- sleuth_results(so, 'treatment_s:full', 'lrt',
                                    show_all = FALSE)

     # filter out the significant genes according to a set qvalue

     sleuth_significant <- dplyr::filter(full_results, qval <= 0.05)

     #view the first 20 genes in this list
     head(sleuth_significant, 20)


     # we can write the entire gene table out
     write.csv(full_results,
               file = "sleuth_results.csv" )

     # we can also use Rshiny to do some interactive visualizations

     sleuth_live(so)


..
	#### Comment: Suggested style guide:
	1. Steps begin with a verb or preposition: Click on... OR Under the "Results Menu"
	2. Locations of files listed parenthetically, separated by carets, ultimate object in bold
	(Username > analyses > *output*)
	3. Buttons and/or keywords in bold: Click on **Apps** OR select **Arabidopsis**
	4. Primary menu titles in double quotes: Under "Input" choose...
	5. Secondary menu titles or headers in single quotes: For the 'Select Input' option choose...
	####
----


**Summary**
~~~~~~~~~~~

Together, Kallisto and Sleuth are quick, powerful ways to analyze RNA-Seq data.

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
