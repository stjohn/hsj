# HSJ: Brachiopods dataset #

Dataset from [Cusack et al 1999](https://doi.org/10.1111/1475-4983.00098):  A moderately-sized study of fossil brachiopods.  Originally coded with non-additive binary coding, it has been recoded using redundant/hierarhical coding for this study.In the original matrix, characters were coded using non-additive binary coding, and we recoded them using redundant coding; the new matrix has 23 taxa and 49 characters, 39 percent of which are secondary characters.


Annotated list of files:
+ [cusack_et_al_1999_newcoding.nex](cusack_et_al_1999_newcoding.nex): Nexus file of 23 taxa and 49 unordered characters.  Missing data coded as '?'; inapplicable characters coded as '-'.
+ [cusack_et_al_1999_prepped_for_fitch.nex](cusack_et_al_1999_prepped_for_fitch.nex):  The dataset with all inapplicable characters coded as '?', prepared for the Fitch approach where inapplicable characters are treated as missing characters (instead of an additional symbol).
+ [cusack_et_al_1999_primaries_only.nex](cusack_et_al_1999_primaries_only.nex):  The dataset restricted to the primary characters only.
+ [cusack_et_al_1999_type.txt](cusack_et_al_1999_type.txt):  The text file that specifies primary characters, secondary characters and which primaries they depend.
+ [cusack_paup_search_primaries_only.tre](cusack_paup_search_primaries_only.tre): The 54,524 best scoring trees from a PAUP search, using primary characters only.
+
[cusack_paup_search.tre](cusack_paup_search.tre): The 373 best scoring trees from a PAUP search with the full matrix.
