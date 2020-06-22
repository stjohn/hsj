# HSJ: Myriapods dataset #

 A character matrix from a large study of living and fossil myriapods [Fernandez et al 2016](https://doi.org/10.1093/sysbio/syw041), modified slightly to use hierarchical coding and removing tertiary characters.

 Annotated list of files:
 + [fernandez_et_al_2016_newcoding.nex](fernandez_et_al_2016_newcoding.nex): Nexus file of 23 taxa and 49 unordered characters.  Missing data coded as '?'; inapplicable characters coded as '-'.
 + [fernandez_et_al_2016_prepped_for_fitch.nex](fernandez_et_al_2016_prepped_for_fitch.nex):  The dataset with all inapplicable characters coded as '?', prepared for the Fitch approach where inapplicable characters are treated as missing characters (instead of an additional symbol).
 + [fernandez_et_al_2016_primaries_only.nex](fernandez_et_al_2016_primaries_only.nex):  The dataset restricted to the primary characters only.
 + [fernandez_et_al_2016_type.txt](fernandez_et_al_2016_type.txt):  The text file that specifies primary characters, secondary characters and which primaries they depend.
 + [fernandez_paup_search_primaries_only.tre](fernandez_paup_search_primaries_only.tre): The 34,560 best scoring trees from a PAUP search, using primary characters only.
 +
 [fernandez_paup_search.tre](fernandez_paup_search.tre): The 2940 best scoring trees from a PAUP search with the full matrix.
