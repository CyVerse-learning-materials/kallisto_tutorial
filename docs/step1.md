\_

![Home\_Icon](./img/homeicon.png){width="25" height="25"}\_ [Learning
Center Home](http://learning.cyverse.org/)

Organize Kallisto Input Data
============================

**Description:**

Kallisto has relatively few input requirements. You will need to have
either single or paired end reads, as well as a reference transcriptome.
It is suggested that your RNA-Seq reads are analyzed using , followed by
any additional trimming and filtering using and application such as .

------------------------------------------------------------------------

\*\*Input <Data:**>

*Importing Reference Transcriptome*
-----------------------------------

In this example, we will import a reference transcriptome for
Arabidopsis from Ensembl. In many cases you can find an appropriate
transcriptome from Ensembl for your organism of interest, or provide
your own fasta-formmatted transcriptome.

> 1.  Go to the Ensembl Plants Arabidopsis page: .
> 2.  In the 'Gene annotation' section, click on the 'Download genes,
>     cDNAs, ncRNA, proteins' 'FASTA' link: .
> 3.  The transcriptome files will be located in the 'cdna' folder: .
> 4.  Right-click/command-click on the
>     Arabidopsis\_thaliana.TAIR10.cdna.all.fa.gz link and copy the link
>     location/URL (usually right-click) to the clipboard
> 5.  Login to the CyVerse .
> 6.  Open the Data window. In your home directory, create a folder to
>     organize your Kallisto project
> 7.  In the created folder, go to the Data window's 'Upload' menu, and
>     select, 'Import from URL...'; paste in the Ensembl link and click
>     'Import from URL' to begin the import

**Output/Results**

------------------------------------------------------------------------

**Description of results and next steps**

This example transcriptome will be indexed in the next step. RNA-Seq
reads will be mapped against this set of transcripts. Once you have the
transcriptome and your RNA-Seq reads, you can proceed with the next
step. We suggest organizing your RNA-Seq reads in the folder created in
step 6 above.

------------------------------------------------------------------------

**Fix or improve this documentation**

-   Search for an answer:
-   Ask us for help: click on the lower right-hand side of the page
-   Report an issue or submit a change:
-   Send feedback: [<Tutorials@CyVerse.org>](Tutorials@CyVerse.org)

------------------------------------------------------------------------

![Home\_Icon](./img/homeicon.png){width="25" height="25"}\_ [Learning
Center Home](http://learning.cyverse.org/)
