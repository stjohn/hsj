# HSJ #

![Image](examples/comp_4taxaCROPPED.jpg)

<sup><sup>Silhouettes modified from [phylopic.org](phylopic.org) (Scutigerella immaculata by  Janssen, Prpic, Damen, & Keesey; Lithobius forficatus by B. Lang).</sup></sup>

Includes code, examples and some data files for computing a scalable disparity metric for datasets with hierarchical characters and scoring phylogenetic trees under the optimality criteria for phylogenetic trees.

## What's here: ##

The repo is organized into 3 folders, corresponding to [code](code), [illustrative examples](examples), and the [data files](data) used in [Hopkins and St. John, 2021]( https://doi.org/10.1093/sysbio/syab005).

## To cite: ##

If you use the disparity metric code, please reference:

+ [Hopkins, M. J. and K. St. John. 2018.](https://doi.org/10.1098/rspb.2018.1784) A new family of dissimilarity metrics for discrete character matrices that include inapplicable characters and its importance for disparity studies. Proceedings of the Royal Society B 285:20181784.

and the Claddis R package:

+ [Lloyd, G. 2016.](https://doi.org/10.1111/bij.12746) Estimating morphological diversity and tempo with discrete character-taxon matrices: implementation, challenges, progress, and future directions. Biological Journal of the Linnean Society 118:131–151.
+ [Lloyd, G. T. 2018.](https://doi.org/10.1111/pala.12380) Journeys through discrete-character morphospace: synthesizing phylogeny, tempo, and disparity. Palaeontology 61:637–645.

If you use the tree scoring code, please reference:

+ [Hopkins, M.J. and K. St. John.  2021.]( https://doi.org/10.1093/sysbio/syab005) Incorporating Hierarchical Characters into Phylogenetic Analysis. Accepted to Systematic Biology.

It is built on the extendible framework of Brazeau et al. 2019 and uses the R packages Claddis, TreeTools, Phangorn, TreeSearch and R scripts for the HSJ approach:

+ [Brazeau, M. D., T. Guillerme, and M. R. Smith. 2019.](https://doi.org/10.1093/sysbio/syy083) An algorithm for Morphological Phylogenetic Analysis with Inapplicable Data. Systematic Biology, 68:619–631.
+ [Brazeau, M. D., M. R. Smith, and T. Guillerme. 2017.](https://zenodo.org/record/815372#.XvEpXS3MyL4) MorphyLib: a library for phylogenetic analysis of categorical trait data with inapplicability. Version 0.0.1-alpha (http://www.morphyproject.org/)
+ [Hopkins, M. J. and K. St. John. 2018.](https://doi.org/10.1098/rspb.2018.1784) A new family of dissimilarity metrics for discrete character matrices that include inapplicable characters and its importance for disparity studies. Proceedings of the Royal Society B, 285:20181784.
+ [Lloyd, G. 2016.](https://doi.org/10.1111/bij.12746) Estimating morphological diversity and tempo with discrete character-taxon matrices: implementation, challenges, progress, and future directions. Biological Journal of the Linnean Society, 118:131–151.
+ [Lloyd, G. T. 2018.](https://doi.org/10.1111/pala.12380) Journeys through discrete-character morphospace: synthesizing phylogeny, tempo, and disparity. Palaeontology, 61:637–645.
+ [Schliep K.P. 2011.](https://doi.org/10.1093/bioinformatics/btq706) phangorn: phylogenetic analysis in R. Bioinformatics, 27(4) 592-593. R package version 2.5.5.
+ [Smith, M. R. 2018.](https://rdrr.io/cran/TreeSearch/) TreeSearch: phylogenetic tree search using custom optimality criteria. R package version 0.4.0.
