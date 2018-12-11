.. include:: cyverse_rst_defined_substitutions.txt

|CyVerse logo|_

|Home_Icon|_
`Learning Center Home <http://learning.cyverse.org/>`_


Analyze Kallisto Results with Sleuth
--------------------------------------

**Description:**

|Sleuth manual| is a program for analysis
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
      - |Example directory of Kallisto results|


*Import Kallisto results into R and analyze with Sleuth*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. tip::

    If you are not familiar with Atmosphere and transfering data into an instance
    see the |Atmosphere guide|

We will import the Kallisto results into an RStudio session being run from
an Atmosphere image. Then we will follow a R script based on the  |Sleuth Walkthoughs|

  1. If necessary, launch the |CyVerse Training Workshop Image|;
     a 'Small1' instance size should be sufficient.

  2. When the instance becomes available, open a web browser and navigate to your
     R studio session. This will be located at URL based on your image ip:
     "http ://your.image.ip:8787";

  3. Login to Atmosphere via Webshell or SSH; initialize iCommands, and use the
     following commands to import the test data and R scripts to your home directory:

     .. code:: bash

        iget -rPV /iplant/home/shared/cyverse_training/tutorials/kallisto/03_output_kallisto_results
        iget -rPV /iplant/home/shared/cyverse_training/tutorials/kallisto/04_sleuth_R

  .. tip::

     If you have problem writing to RStduio libraries, use the Atmosphere web-webshell
     or connect to Atmosphere via ssh and run the following command. See the
     |Atmosphere Guide| for more info on connecting to your instance.

        .. code:: bash

         sudo chown -R $USER /usr/local/lib/R/site-library

    Your CyVerse username is your sudo password.

  4. Load the |sample_kallisto_script.R|
     into your R Studio session. Use the |Kallisto_demo_tsv| file for the step
     where you describe experimental design metadata.


.. code:: R

    ---
    title: "Sleuth RNA-Seq Tutorial - Arabidopsis"
    output:
    html_document:
    df_print: paged
    ---

    This tutorial will take you through a sample RNA-Seq analysis using performed by [kallisto](https://pachterlab.github.io/kallisto/about), using an RNA-Seq R package [Sleuth](https://pachterlab.github.io/sleuth/about). This tutorial is based on the one by the [Pachter lab](https://pachterlab.github.io/sleuth_walkthroughs/boj/analysis.html).

    ## Install Kallisto

    Installation instructions for Kallisto are [available here](https://pachterlab.github.io/sleuth/download). This make take a few minutes so we have done most of these items for you. If you were running this on your own installation of R, you would need to uncomment and run all these tools.

    ```{r message=FALSE}

    #source("http://bioconductor.org/biocLite.R")
    #biocLite("rhdf5")
    #install.packages("devtools")
    #install.packages("rlang")
    #devtools::install_github("pachterlab/sleuth")
    ```


    Next, we need to load the Sleuth library to begin. We will also check the version:
    ```{r message=FALSE}
    require("sleuth")
    packageVersion("sleuth")
    ```

    ## Locate sample names and describe our experimental design

    We need to provide Sleuth with our sample names:

    ```{r}
    sample_id <- dir(file.path("~/03_output_kallisto_results/"))
    sample_id
    ```

    We also need to get the file paths to our results files.
    ```{r}
    kal_dirs <- file.path("~/03_output_kallisto_results", sample_id)
    ```

    We also need a table that provides more meaningful names for describing our experiment...

    ```{r}
    s2c <- read.table(file.path("~/04_sleuth_R/kallisto_demo.tsv"),
                 header = TRUE,
                 stringsAsFactors = FALSE,
                 sep = "\t")
    ```



    We will add our filepaths to the table

    ```{r}
    s2c <- dplyr::mutate(s2c, path = kal_dirs)
    ```

    Let's view the table we have created:
    ```{r}
    s2c
    ```

    We will also install a plotting library and some other functions we will need...

    ```{r}
    require("gridExtra")
    install.packages("cowplot")
    library("cowplot")
    ```

    ## Get gene names from biomaRt

    Before we look at our data, lets use [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) tools will allow us to pull in recognizable gene names from a database.

    First let's install and load BiomRart. **(this will take a few minutes!)**.

    ```{r echo=FALSE, message=FALSE, warning=FALSE}
    source("https://bioconductor.org/biocLite.R")
    biocLite("biomaRt", suppressUpdates=TRUE)
    ```

    Next we need to determine which biomaRt to use. This can be a little complex so
    be sure to read their [documentation](https://www.bioconductor.org/packages/devel/bioc/vignettes/biomaRt/inst/doc/biomaRt.html) and this [blog post](https://nsaunders.wordpress.com/2015/04/28/some-basics-of-biomart/) is also helpful.

    ```{r}
    library(biomaRt)
    marts <- listMarts()
    marts
    ```

    If you are not working with these Ensembl data bases you may want to check out documentation on [using BiomaRts other than Ensembl](https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html#using-a-biomart-other-than-ensembl). We are using plants, so

    ```{r}
    marts <- listMarts(host = "plants.ensembl.org")
    marts
    ```

    For now, remember that we will want to use `plants_mart`.

    Next, we need to choose a specific dataset.

    ```{r}
    plants_mart <- useMart("plants_mart", host = "plants.ensembl.org" )
    listDatasets(plants_mart)
    ```

    After a little looking, its the `athaliana_eg_gene` dataset that we need. Finally, we need to update our `plants_mart` to be more specific.

    ```{r}
    plants_mart <- useMart("plants_mart", dataset = "athaliana_eg_gene", host="plants.ensembl.org" )
    ```

    Now we want to get specific attributes from the list of genes we can import from biomart

    ```{r}
    listAttributes(plants_mart)
    ```

    We can choose whichever of these we'd like to use. Let's get transcript ids, gene ids, a description, and gene names. Notice, there are many things you may
    want to come back for. We must get the transcript id because these are the names of the transcripts that were used in our Kallisto quantification.

    ```{r}
    t2g <- getBM(attributes = c("ensembl_transcript_id",
                           "ensembl_gene_id",
                           "description",
                           "external_gene_name"),
            mart = plants_mart)
    ```

    We need to make sure the `ensembl_transcript_id` column is named `target_id`

    ```{r}
    ttg <- dplyr::rename(t2g, target_id= ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
    ```


    ### Prepare data for Sleuth

    first we need to alter our experimental design so that we consider the full transcriptome sample to be the "control" to compare to...

    ```{r}
    s2c$genotype_variation_s <- as.factor(s2c$genotype_variation_s)
    s2c$genotype_variation_s <- relevel(s2c$genotype_variation_s, ref = "wild type")
    ```


    Now we need to tell Sleuth both about the Kallisto results and the gene names (and gene descriptions/metadata) we obtained from biomaRt. The `sleuth_prep` function does this.

    ```{r}
    so <- sleuth_prep(s2c,
                 full_model = ~genotype_variation_s,
                 target_mapping = ttg,
                 read_bootstrap_tpm=TRUE,
                 extra_bootstrap_summary = TRUE)
    ```


    ### Examine Sleuth PCA

    Next, we should check to see if our samples (and replicates) cluster on a PCA (as should expect) or if there are outliers:

    ```{r}
    plot_pca(so, color_by = 'genotype_variation_s', text_labels = TRUE)
    ```

    We can also see genes involved in the the 1st PC by looking at the loadings (primary genes whose linear combinations define the principal components)

    ```{r}
    plot_loadings(so, pc_input = 1)
    ```

    ## Differential expression testing with Sleuth

    Now we need to run a few functions that will test for differential expression (abundance).

    First we will create a model

    ```{r}
    so <- sleuth_fit(so, ~genotype_variation_s, 'full')
    so <- sleuth_fit(so, ~1, 'reduced')
    so <- sleuth_lrt(so, 'reduced', 'full')
    ```

    Now we can get the results of this analysis

    ```{r}
    full_results <- sleuth_results(so, 'reduced:full', 'lrt',
                              show_all = FALSE)
    head(full_results)
    ```



    Let's add  Wald test
    ```{r}
    wald_test <- colnames(design_matrix(so))[2]
    so <- sleuth_wt(so, wald_test)
    ```

    And start a Shiny Browser

    ```{r}
    sleuth_live(so)
    ```

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
Together, Kallisto and Sleuth are quick, powerful ways to analyze RNA-Seq data.

----

Additional information, help
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

.. |Sleuth manual| raw:: html

   <a href="https://pachterlab.github.io/sleuth/about" target="blank">Sleuth</a>

.. |Example directory of Kallisto results| raw:: html

    <a href="http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/03_output_kallisto_results" target="_blank">Example directory of Kallisto results</a>

.. |Sleuth Walkthoughs| raw:: html

    <a href="https://pachterlab.github.io/sleuth/walkthroughs" target="_blank">Sleuth Walkthoughs</a>

.. |CyVerse Training Workshop Image| raw:: html

   <a href="https://atmo.cyverse.org/application/images/1479" target="blank">CyVerse Training Workshop Image</a>

.. |sample_kallisto_script.R| raw:: html

    <a href="http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/04_sleuth_R/sample_kallisto_script.R" target="_blank">sample_kallisto_script.R</a>

.. |Kallisto_demo_tsv| raw:: html

    <a href="http://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/tutorials/kallisto/04_sleuth_R/kallisto_demo.tsv" target="_blank">Kallisto_demo_tsv</a>
