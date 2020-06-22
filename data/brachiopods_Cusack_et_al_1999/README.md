# HSJ: Brachiopods dataset #


Dataset from [Cusack et al 1999](https://doi.org/10.1111/1475-4983.00098):  A moderately-sized study of fossil brachiopods.  Originally coded with a mix of non-additive binary and composite coding, it has been recoded using reductive coding [sensu Strong and Lipscomb 1999](https://doi.org/10.1006/clad.1999.0114) for this study. For example, the first characters in the original matrix coded for the presence of a pitted surface in the larval shell (as opposed to smooth). The subsequent five characters described variation in the expression of the pits, but coded any taxon with a smooth shell as "absent" for these characters (which overweights the absence of pits because it is effectively coded six times). We removed the "absent" character state from these characters and considered all taxa with smooth shells as "inapplicable" for this character.  Thus the first character was designated a controlling primary character with five secondary characters contigent on it. Similarly, other multistate traits that included an "absent" state and more than one states describing variation in the trait when present were broken into two or more characters, one defining the absence or presence of the trait and the rest describing the expression of the trait when present.  The final dataset has 23 taxa and 49 characters, 29 of which are primary characters.


Annotated list of files:
+ [cusack_et_al_1999_newcoding.nex](cusack_et_al_1999_newcoding.nex): Nexus file of 23 taxa and 49 unordered characters.  Missing data coded as '?'; inapplicable characters coded as '-'.
+ [cusack_et_al_1999_prepped_for_fitch.nex](cusack_et_al_1999_prepped_for_fitch.nex):  The dataset with all inapplicable characters coded as '?', prepared for the Fitch approach where inapplicable characters are treated as missing characters (instead of an additional symbol).
+ [cusack_et_al_1999_primaries_only.nex](cusack_et_al_1999_primaries_only.nex):  The dataset restricted to the primary characters only.
+ [cusack_et_al_1999_type.txt](cusack_et_al_1999_type.txt):  The text file that specifies primary characters, secondary characters and which primaries they depend.
+ [cusack_paup_search_primaries_only.tre](cusack_paup_search_primaries_only.tre): The 54,524 best scoring trees from a PAUP search, using primary characters only.
+ [cusack_paup_search.tre](cusack_paup_search.tre): The 373 best scoring trees from a PAUP search with the full matrix.
