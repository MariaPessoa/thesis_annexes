# thesis_annexes
## Software developed and main results obtained in the Master dissertation in Bioinformatics: 
<br>A model validation pipeline for healthy tissue genome-scale metabolic models. 

## Abstract:
<br>In the past few years, high-throughput experimental methods have made omics data available for several
layers of biological organization, enabling the integration of knowledge from individual components into
complex modelssuch as genome-scale metabolic models (GSMMs). These can be analysed by constraintbased modelling (CBM) methods, 
which facilitate in silico predictive approaches.

Human metabolic models have been used to study healthy human tissues and their associated
metabolic diseases, such as obesity, diabetes, and cancer. Generic human models can be integrated with
contextual data through reconstruction algorithms to produce context-specific models (CSMs), which are
typically better at capturing the variation between different tissues and cell types. As the human body
contains a multitude of tissues and cell types, CSMs are frequently adopted as a means to obtain accurate
metabolic models of healthy human tissues.

However, unlike microorganisms’ or cancer models, which allow several methods of validation
such as the comparison of in silico fluxes or gene essentiality predictions to experimental data, the
validation methods easily applicable to CSMs of healthy human tissue are more limited. Consequently,
despite continued efforts to update generic human models and reconstruction algorithms to extract high
quality CSMs, their validation remains a concern.

This work presents a pipeline for the extraction and basic validation of CSMs of normal human
tissues derived from the integration of transcriptomics data with a generic human model. All CSMs were
extracted from the Human1 generic model recently published by Robinson et al. (2020), relied on the
open-source Troppo Python package and in the fastCORE and tINIT reconstruction algorithms
implemented therein. CSMs were extracted for 11 healthy tissues available in the GTEx v8 dataset.

Prior to extraction, machine learning methods were applied to threshold selection for gene scores
conversion. The highest quality models were obtained with a global threshold applied to the omics data
directly. The CSM validation strategy focused on the total number of metabolic tasks passed as a
performance indicator. Lastly, this work is accompanied by Jupyter Notebooks, which include a beginner
friendly model extraction guide.
