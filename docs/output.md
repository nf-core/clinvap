# nf-core/clinvap: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [VEP](#vep) - Variant effect prediction via Ensembl VEP
* [Report Generation](#report) - Aggregate clinical variant annottion report

## VEP

[VEP](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html) is used to annotate raw VCF files. The annotation includes the predicted effects of the mutations on the protein function.

**Output directory: `results`**

* `<VCF>.vcf`
  * Annotated VCF file

## Report Generation

[Report Generation](#report) processes the information found in the annotated VCF file and conducts further annotations including driver gene type and therapeutic suggestions.
The resulting json report includes 6 main categories:

1. "mskdg": Somatic Mutations in Known Driver genes  
    List of cancer driver genes along with the observed mutations in the patient.
2. "ptp_da": Summary of Cancer Drugs Targeting the Affected genes  
    List of cancer drugs targeting the mutated gene.
3. "ptp_ia": CIViC Summary of Drugs Targeting the Affected genes  
    Therapies that have evidence of targeting the affected gene
4. "mskpe": Somatic Mutations with Known Pharmacogenetic Effect
    List of drugs that directly targets the observed variant of the gene  
5. "ref": References  
    The publications as an evidence of the found annotation/information
6. "appendix": Appendix  
    All the somatic variants of the patient with their dbSNP and COSMIC IDs.

**Output directory: `results`**

* `<REPORT>.json`
  * Clinical report in json format.
