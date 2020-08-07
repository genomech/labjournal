# Report BGRS-2020

Mechanisms of chromatin organization are complicated and not completely understood.​



For prediction EP and other genomic interactions the following epigenetic data were used: 
* ChIP-seq profiles: describing chromatin binding of architectural proteins;
* RNA-seq profiles, which show gene expression;
* E1 values, which show compartments A and B;
* genomic distance. The farther regions apart, the less chance of contact. The algorithm was restricted to mid-range contacts (less than 1.5Mb), because almost all EP contacts lie within this distance.

...

Prediction can be significantly improved by training algorithm on higher resolution data.




## Abstract

Alterations of spatial contacts of chromatin may lead to developmental disorders and cancer.
We have recently shown that 3-dimensional genomic contacts could be predicted in normal and rearranged genomes using epigenetic information.
Here we describe an extension of this algorithm, characterized by more accurate and higher-resolution predictions, and its web implementation which allows researchers and clinical doctors to generate patient-specific predictions of genome architecture.

## 3D organization

Важность 3D организации в целом

## Изменчивость 3D организации

## Motivation and Aim

### Motivation

It was recently shown that chromatin contacts could be predicted using widely-available epigenetic information and machine-learning methods.
These predictions could be used to dissect molecular mechanisms which connect chromosomal rearrangements with human diseases, and thus could provide useful resource for clinical genetics.
However, the majority of the developed predictive algorithms are not available as a stand-alone or online software;
therefore, their usage in clinics is limited.
Moreover, the resolution of these algorithms is often low, thus estimation of the clinical significance of changes in 3D-genome architecture is challenging.

### Aim

We aimed to extend our recently developed computational tool, 3DPredictor, to predict chromatin contacts at high resolution (1kb).
Moreover, we aimed to develop a web application that allows to predict genome organization in normal and rearranged genomes, for researchers and clinical doctors who are not experienced in bioinformatics.

## Methods

We used ML-based algorithm to train models and predict chromatin 3D organization.
Training datasets which high-resolution prediction required were taken from *Krietenstein et al.*

For web service development we used the following tools:

1. Front-end: PureCSS framework, PHP
2. Back-end:
	* Pipeline: Bash script
	* Preprocessing CTCF data: gimmemotifs
	* Preprocessing RNA-seq data: Python3 script
	* Hi-C map: Juicer Tools

## Results

We showed that prediction can be significantly improved by training algorithm on higher resolution data.
Also, we made a web service which performs prediction of chromatin 3D-structure for a region of interest.
Trained models for *H. sapiens* and *M. musculus* at 1- and 5-kb resolutions are provided.
Web-3DPredictor is hosted on the ICG server, available [here](https://genedev.bionet.nsc.ru/Web\_3DPredictor).
Both source code and trained models could be obtained on GitHub to deploy local server elsewhere.

16:9
