# About what's here

This directory contains data files and source code for analyzing gene models from:

* Clark, et al. (2019) [Expanding Alternative Splicing Identification by Integrating Multiple Sources of Transcription Data in Tomato](https://www.ncbi.nlm.nih.gov/pubmed/31191588)

Note that the annotations produced in this article used RNA-Seq data aligned to the SL3.0 Heinz genome assembly. The spliced aligner `tophat` was used with default parameters.
This was problematic because the default intron size parameter for `tophat` (as with most spliced aligners) is too large. 
As a result, many alignments spanned neighboring genes, which introduced false positives. 

To see some examples, view the image files in the `results` directory.

To see the Clark et al. gene models aligned to the SL3.0 genome:

* Start [Integrated Genome Browser](https://bioviz.org).
* Open the SL3.0 (S_lycopersicum_Feb_2017) genome assembly (use the **Current Genome** tab).
* Select **Clark et al. 2019** in the **Available Data** section of the **Data Access** panel.
* Click **Load Data** button to load the data.





