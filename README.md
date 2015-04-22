README
======

Perturbations to the gut microbiota result in a loss of colonization resistance against gastrointestinal pathogens such as *Clostridium difficile*. Although *C. difficile* infection is commonly associated with antibiotic use, the precise alterations to the microbiota associated with this loss in function are unknown. We used a variety of antibiotic perturbations to generate a diverse array of gut microbiota structures, which were then challenged with *C. difficile* spores. Across these treatments we observed that *C. difficile* resistance was never attributable to a single organism, but rather it was the result of multiple microbiota members interacting in a context-dependent manner. Using relative abundance data, we built a machine learning regression model to predict the levels of *C. difficile* that were found 24 hours after challenging the perturbed communities. This model was able to explain 77.2% of the variation in the observed number of *C. difficile* per gram of feces. This model revealed important bacterial populations within the microbiota, which correlation analysis alone did not detect. Specifically, we observed that populations associated with the *Porphyromonadaceae*, *Lachnospiraceae*, *Lactobacillus*, and *Alistipes* were protective and populations associated with *Escherichia* and *Streptococcus* were associated with high levels of colonization. In addition, a population affiliated with *Akkermansia* indicated a strong context dependency on other members of the microbiota. Together, these results indicate that individual bacterial populations do not drive colonization resistance to *C. difficile*. Rather, multiple diverse assemblages act in concert to mediate colonization resistance.



Overview
--------

    project
    |- README          # the top level description of content
    |
    |- doc/            # documentation for the study
    |  |- notebook/    # preliminary analyses (dead branches of analysis)
    |  +- paper/       # manuscript(s), whether generated or not
    |
    |- data            # raw and primary data, are not changed once created
    |  |- references/  # reference files to be used in analysis
    |  |- raw/         # raw data, will not be altered
    |  +- process/     # cleaned data, will not be altered once created
    |
    |- code/           # any programmatic code
    |- results         # all output from workflows and analyses
    |  |- tables/      # text version of tables to be rendered with kable in R
    |  |- figures/     # graphs, likely designated for manuscript figures
    |  +- pictures/    # diagrams, images, and other non-graph graphics
    |
    |- scratch/        # temporary files that can be safely deleted or lost
    |
    |- study.Rmd       # executable Rmarkdown for this study, if applicable
    |- study.md        # Markdown (GitHub) version of the *Rmd file
    |- study.html      # HTML version of *.Rmd file
    |
    +- Makefile        # executable Makefile for this study, if applicable
