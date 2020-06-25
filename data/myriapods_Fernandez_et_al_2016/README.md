# HSJ: Myriapods dataset #

 A character matrix from a large study of living and fossil myriapods [Fernandez et al 2016](https://doi.org/10.1093/sysbio/syw041). The original matrix was coded using a reductive coding strategy [sensu Strong and Lipscomb 1999](https://doi.org/10.1006/clad.1999.0114). We modified it slightly so that the characters were coded hierarchically as expected by the HSJ approach [Appendix 3; Hopkins and St. John 2018](https://doi.org/10.1098/rspb.2018.1784). For this study, we also removed tertiary characters and modified a few characters to meet the assumption that all controlling primary characters are binary with the "0" token representing the character state for which contingent secondary characters are inapplicable. The dataset has 47 taxa and 201 characters, 130 of which are primary characters.

 Annotated list of files:
 + [fernandez_et_al_2016_newcoding.nex](fernandez_et_al_2016_newcoding.nex): Nexus file of 47 taxa and 130 characters.  Missing data coded as '?'; inapplicable characters coded as '-'.
 + [fernandez_et_al_2016_prepped_for_fitch.nex](fernandez_et_al_2016_prepped_for_fitch.nex):  The dataset with all inapplicable characters coded as '?', prepared for the Fitch approach where inapplicable characters are treated as missing characters (instead of an additional symbol).
 + [fernandez_et_al_2016_primaries_only.nex](fernandez_et_al_2016_primaries_only.nex):  The dataset restricted to the primary characters only.
 + [fernandez_et_al_2016_type.txt](fernandez_et_al_2016_type.txt):  The text file that specifies primary characters, secondary characters and which primaries they depend.
 + [fernandez_paup_search_primaries_only.tre](fernandez_paup_search_primaries_only.tre): The 34,560 best scoring trees from a PAUP search, using primary characters only.
