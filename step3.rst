.. include:: cyverse_rst_defined_substitutions.txt
.. include:: custom_urls.txt

|CyVerse_logo|_

|Home_Icon|_
`Learning Center Home <http://learning.cyverse.org/>`_


Analyze Kallisto Results with Sleuth
--------------------------------------

**Description:**

|Sleuth manual| is a program for analysis of RNA-Seq experiments for which
transcript abundances have been quantified with Kallisto. In this tutorial, we
will use R Studio being served from an VICE instance.

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



We will import the Kallisto results into an RStudio session being run from
an Atmosphere image. Then we will follow a R script based on the  |Sleuth Walkthoughs|


  1. If necessary, login to the CyVerse |discovery_enviornment|.

  2. In the App panel, open the **Sleuth RStudio** app or click this link:
     |Sleuth app|

  3. Name your analysis, and if desired enter comments.

  4. (Optional) In the 'Notebooks' section, under 'Select an RMarkdown
     notebook to run' select a notebook.

     .. admonition:: Sample data

       For the sample data, navigate to and select
       **/iplant/home/shared/cyverse_training/tutorials/kallisto/04_sleuth_R/sleuth_tutorial.Rmd**

  5. In the 'Datasets' section, under 'Data for analysis (outputs of Kallisto
     quantification)' choose the folders containing quantification information for all sets of reads.

     .. admonition:: Sample data

       For the sample data, navigate to and select
       **/iplant/home/shared/cyverse_training/tutorials/kallisto/03_output_kallisto_results**

  6. In the 'Datasets' section, under 'Study design file' choose a TSV file
     describing the samples and study design (see |Sleuth manual|).

     .. admonition:: Sample data

       For the sample data, navigate to and select
       **/iplant/home/shared/cyverse_training/tutorials/kallisto/04_sleuth_R/kallisto_demo.tsv**

    .. tip::

       See the Example study design (|Kallisto_demo_tsv|) TSV file. You can
       create and edit your own in a spreadsheet editing program. The |Sleuth
       manual| explains this file and more is described in this tutorial's RMarkdown notebook.

  7. Click 'Launch Analyses' to start the job. Click on the Analyses button
     to monitor the job and results. In your notifications, you will find a
     link to your VICE session ("Access your running analyses here"); this may
     take a few minutes to become active.


  8. In your RStudio session, double click on the **sleuth_tutorial.Rmd** file
     and follow the tutorial by pressing the green "play" triangles in each
     section of code. The code for the notebook is replicated below:



   .. code:: R

              ---
        title: "Sleuth RNA-Seq Tutorial - Arabidopsis"
        output:
          html_document:
            df_print: paged
        ---

        This tutorial will take you through a sample RNA-Seq analysis using [kallisto](https://pachterlab.github.io/kallisto/about), using an RNA-Seq R package [Sleuth](https://pachterlab.github.io/sleuth/about). This tutorial is based on the one by the [Pachter lab](https://pachterlab.github.io/sleuth_walkthroughs/boj/analysis.html)

        ## Step 1: Load Sleuth and accessory libraries

        Next, we need to load the Sleuth library to begin. We will also check the version:


        ```{r message=FALSE}
        require("sleuth")
        packageVersion("sleuth")
        ```

        We will also install a plotting library and some other functions we will need...

        ```{r echo=FALSE, message=FALSE, warning=FALSE}
        library("gridExtra")
        library("cowplot")
        ```

        We will also use [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) tools will allow us to pull in recognizable gene names from a database.


        ```{r echo=FALSE, message=FALSE, warning=FALSE}
        library("biomaRt")
        ```

        ## Step 2: Load experimental design and label kallisto outputs with metadata

        ### Locate sample names and describe our experimental design

        We need to provide Sleuth with our sample names:

        ```{r}
        sample_id <- dir(file.path("~/kallisto_qaunt_output/"))
        sample_id
        ```

        We also need to get the file paths to our results files.
        ```{r}
        kal_dirs <- file.path("~/kallisto_qaunt_output", sample_id)
        ```

        We also need a table that provides more meaningful names for describing our experiment...

        ```{r}
        s2c <- read.table(file.path("~/kallisto_demo.tsv"),
                          header = TRUE,
                          stringsAsFactors = FALSE,
                          sep = "\t")
        ```

        We will add our file paths to the table

        ```{r}
        s2c <- dplyr::mutate(s2c, path = kal_dirs)
        ```

        Let's view the table we have created:
        ```{r}
        s2c
        ```

        ## Step 3: Load gene names from Ensembl

        Next we need to determine which biomaRt to use. This can be a little complex so be sure to read their [documentation](https://www.bioconductor.org/packages/devel/bioc/vignettes/biomaRt/inst/doc/biomaRt.html). This [blog post](https://nsaunders.wordpress.com/2015/04/28/some-basics-of-biomart/) is also helpful.

        ```{r}
        marts <- listMarts()
        marts
        ```

        If you are not working with these Ensembl data bases you may want to check out documentation on [using BiomaRts other than Ensembl](https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html#using-a-biomart-other-than-ensembl). We are using plants, so...

        ```{r}
        marts <- listMarts(host = "plants.ensembl.org")
        marts
        ```

        For now, remember that we will want to use `plants_mart`.

        Next, we need to choose a specific data set.

        ```{r}
        plants_mart <- useMart("plants_mart", host = "plants.ensembl.org" )
        listDatasets(plants_mart)
        ```

        After a little looking, its the `athaliana_eg_gene` data set that we need. Finally, we need to update our `plants_mart` to be more specific.

        ```{r}
        plants_mart <- useMart("plants_mart", dataset = "athaliana_eg_gene", host="plants.ensembl.org" )
        ```

        Now we want to get specific attributes from the list of genes we can import from biomart

        ```{r}
        listAttributes(plants_mart)
        ```

        We can choose whichever of these we'd like to use. Let's get transcript ids, gene ids, a description, and gene names. Notice, there are many things you may
        want to come back for. We must get the transcript id because these are the names of the transcripts that were used in our Kallisto quantification.

        ```{r echo=FALSE, message=FALSE, warning=FALSE}
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


        ##Step 4: Prepare data for Sleuth

        first we need to alter our experimental design so that we consider the full transcriptome sample to be the "control" to compare to...

        ```{r}
        s2c$genotype_variation_s <- as.factor(s2c$genotype_variation_s)
        s2c$genotype_variation_s <- relevel(s2c$genotype_variation_s, ref = "wild type")

        ```


        Now we need to tell Sleuth both about the Kallisto results and the gene names (and gene descriptions/metadata) we obtained from biomaRt. The `sleuth_prep` function does this.

        ```{r warning=FALSE}
        so <- sleuth_prep(s2c,
                     full_model = ~genotype_variation_s,
                     target_mapping = ttg,
                     read_bootstrap_tpm=TRUE,
                     extra_bootstrap_summary = TRUE)
        ```

        ##Step 5: Initial data exploration

        ### Examine Sleuth PCA

        Next, we should check to see if our samples (and replicates) cluster on a PCA (as should expect) or if there are outliers. When we plot by condition, we'd expect that similar colors group together.

        ```{r}
        library(cowplot)
        ggplot2::theme_set(theme_cowplot())
        plot_pca(so, color_by = 'genotype_variation_s', text_labels = TRUE)
        ```

        Let's try plotting by treatment


        ```{r}
        plot_pca(so, color_by = 'treatment_s', text_labels = TRUE)
        ```


        We can also see genes involved in the the 1st PC by looking at the loadings (primary genes whose linear combinations define the principal components)


        ```{r}
        plot_loadings(so, pc_input = 1)
        ```

        Let's see how this "influential" gene (at least as far as PCA tells us) looks by condition
        ```{r}
        plot_bootstrap(so, 'AT2G34420.1', color_by = 'genotype_variation_s')
        ```

        Let's see how this "influential" gene (at least as far as PCA tells us) looks by treatment

        ```{r}
        plot_bootstrap(so, 'AT2G34420.1', color_by = 'treatment_s')
        ```


        ##Step 6: Modeling, testing, and results exploration

        ### Differential expression testing with Sleuth

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














----


**Summary**
Together, Kallisto and Sleuth are quick, powerful ways to analyze RNA-Seq data.


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
