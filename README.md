# Genome for B96 maize inbred line assembly using PacBio technology

This is a repository for B96 genome assembly, scaffolding, and annotation processes. B96, a maize inbred line, is highly resistant to multiple herbivores, including European corn borer (<i>Ostrinia nubilalis</i>), two spotted spider mites (<i>Tetranychus urticae</i>) and complex species of thrips. Transcriptome study for B96 compared to the more susceptible inbred line B73 has elucidated a large number defense genes being highly overexpressed under constitutive and/or induced conditions ([Ji et al.]()). By sequencing the B96 genome using long-sequencing technology, we provided a high-quality genome resouce for further herbivore-resistance study in maize and to facilitate the development of a more sustainable agriculture.

## Table of content

- Overview of PacBio-sequencing data[Overview-of-PacBio-sequencing-data]
- <i>De novo</i> assembly using long-sequencing data[#De-novo-assembly-using-long-sequencing-data]
- Genome scaffolding based on reference genome B73 v5.0[#Genome-scaffolding-based-on-reference-genome-B73-v5.0]
- Gene annotation pipeline using transcriptome evidence[#Gene-annotation-pipeline-using-transcriptome-evidence]
- 

## Overview of PacBio-sequencing data


## <i>De novo</i> assembly using long-sequencing data

## Genome scaffolding based on reference genome B73 v5.0

## Gene annotation pipeline using transcriptome evidence
We used repeat-masked genome for gene annotation. Given the fact that some gene space (with RNA-evidence) might be masked from gene model prediction, we developed a script to correct genome spaces which have RNA support by taking RNA-seq alignment file. 
```bash
# use the unmask_RNA.py (parallel process)
mpiexec -n 10 python unmask_RNA.py -fasta <fasta> -bam <RNA_bam> -O <output>
```
##
