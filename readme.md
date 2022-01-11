# greporeseq: The GREPore-seq Analysis Package

## Table of contents

* [Features](#Features)
* [Dependencies](#Dependencies)
* [Getting Set Up](#Getting-Set-Up)
  
  * [Install Miniconda2](#Install-Miniconda2)
  * [Install Dependencies](#Install-Dependencies)
* [Running the Full Analysis Pipeline](#Running-the-Full-Analysis-Pipeline)
  
  * [Quickstart](#Quickstart)
  * [Writing A Manifest Excel](#Writing-A-Manifest-file)
  * [A Full Manifest File Example](#A-Full-Manifest-File-Example)
  * [Pipeline Output](#Pipeline-Output)
  
  

## Features

The package implements a pipeline consisting of a read preprocessing module followed by a visualization module. The preprocessing module takes raw reads(FASTQ) from a pooled multi-sample Oxford Nanopore sequencing run as input. Reads are demultiplexed into sample-specific FASTQs using Grepseq information.

[outline picture]

The individual pipeline steps are:

1. **Make_reference**: Based on the information entered in the RefInfo.yaml file, produce the corresponding reference sequence FASTA files.
2. **Demultiplex**: A multi-sample Oxford Nanopore sequencing run is demultiplexed into sample-specific read FASTQ files based on the information in DemultiplexInfo.yaml.
3. **Visualization**: The demultiplexed read files are aligned to reference using Minimap2-ax map-ont algorithm with default parameters ([Li. H, 2018](https://doi.org/10.1093/bioinformatics/bty191)), then sorted and indexed by using samtools([Petr Danecek, 2021](https://doi.org/10.1093/gigascience/giab008)).

After all the above steps have been completed, Results can be viewed using IGV.

[![igv_panel.png](https://github.com/lisiang/GREPore-seq/blob/master/igv_panel.png)]

## Dependencies

* Python(3.6+)
* porechop
* samtools==1.10
* minimap2
* pigz
* yaml

## Getting Set Up

### Install Miniconda2

We recommend you use Anaconda to install all dependencies:

Install conda, take miniconda2 as an example:

1. Download the Miniconda installer to your Home directory.
   `wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh -O ~/miniconda.sh` 
2. Install Miniconda quietly, accepting defaults, to your Home directory.
   `bash ~/miniconda.sh -b -p`
3. Remove the Miniconda installer from your Home directory.
   `rm ~/miniconda.sh`
4. Test Miniconda install
   `source ~/.bashrc`
   `conda --help`

### Install Dependencies

Use conda install all dependencies:

1. Make a conda virtual environment for greporeseq.

   `conda create -n greporeseq`

2. Activate the conda greporeseq environment.
   `conda activate greporeseq`

3. Install greporeseq dependencies by entering the greporeseq directory and running
   `conda install --file requirements.txt -y`

## Running the Full Analysis Pipeline

### Quickstart

To run the full greporeseq analysis pipeline, you must first have create 2 manifest YAML file that describes all pipeline inputs. Once you have done so, you can simply run

`python /path/to/greporeseq.py all -n ***.fastq.gz -r ***RefInfo.yaml -d ***DemultiplexInfo.yaml`

to run the entire pipeline. Below are specific instructions detailing how to write the manifest file.

If you wish to run an example on our abridged test data, You can do it by the following command

`python ./greporeseq/greporeseq.py all -n ./test/N26chop.fastq.gz -r ./test/N26RefInfo.yaml -d ./test/N26DemultiplexInfo.yaml`

### Writing A Manifest file

When running all steps of greporeseq, it is necessary to describe each sample that needs to be demultiplex. YAML file can easily store all the input information and be read by greporeseq.

the fields contained in the `DemultiplexInfo.yaml` are:

* **reference_id**: Used to fill in the reference id. Note: the reference_id must be from the RefInfo.
* **left_150bp**: Used to fill in the 150bp sequence at the left end of the sequenced sample.
* **right_150bp**: Used to fill in the 150bp sequence at the right end of the sequencing sample.
* **BCprimer_F**: Used to fill in the sequence of the forward Barcode + forward primer. Required if the sequencing sample is a PCR product. Note: If there is no Barcode, only the primer sequence is required.
* **BClen_F**: Used to fill in the number that refers to the length of the Barcode. Note: If there is no Barcode, enter the number 0.
* *optional*
* **BCprimer_R**: Used to fill in the sequence of the reverse Barcode + reverse primer. Note: If there is no Barcode, only the primer sequence is required.
* **BClen_R**: Used to fill in the number that refers to the length of the Barcode. Note: If there is no Barcode, enter the number 0.
* **unique_sequence**: For filling in sample-specific sequences. Improves the accuracy of demultiplexing.

An example `DemultiplexInfo.yaml`:

```yaml
[SAMPLE_ID]: 
  reference_id: 4kAAVS1
  left_150bp: tgcaaacaggaagtgaacggggaagggagggggcttctcatctgggtgcgggaaccccacatggtacctgttagacacggcaaaacccccgtcaccacccacaggtggcgcttccagtgctcagactagggaagaggttccagcccctcc
  right_150bp: aaggagacaaagtccaggaccggctggaggggctcaacatcggaagaggggaagtcgagggagggatggtaaggaggactgcatgggtcagcacaggctgccaaagccagggccagttaaagcgactccaatgcggaagagagtaggtcg
  BCprimer_F: atgataactaggTGCAAACAGGAAGTGAACGG
  BClen_F: 12
#optional
  BCprimer_R: 
  BClen_R: 
  unique_sequence: 
```

the field contained in the `RefInfo.yaml` :

**sequence**: Used to fill in the reference sequence. Note: BC sequences are usually not included in the reference sequence.

An example `RefInfo.yaml`:

```yaml
[REFERENCE_ID]: 
  sequence: [REFERENCE_SEQUENCE]
```

When too many samples need to be analysed at once, it can be very time consuming to enter the manifest YAML file manually, so we have provided an EXCEL sheet, `GREPore-seq_Template.us.xlsx`, to **quickly produce manifest manifest YAML file**. Fill in the appropriate columns with the input information. **Last column of the sheet, the YAML column**, contains formulas that will automatically convert the information into YAML format. **Copy and paste into a .yaml file, removing the double quotes before and after,** to obtain a correctly formatted manifest YAML description file.

### **A Full Manifest File Example**

Below is an example of a full manifest file.

`DemultiplexInfo.yaml`

```yaml
N26_WW1590_sg25_4kBC26_2020_06_19_853T4d_wt_4d: 
  reference_id: 4kAAVS1
  left_150bp: tgcaaacaggaagtgaacggggaagggagggggcttctcatctgggtgcgggaaccccacatggtacctgttagacacggcaaaacccccgtcaccacccacaggtggcgcttccagtgctcagactagggaagaggttccagcccctcc
  right_150bp: aaggagacaaagtccaggaccggctggaggggctcaacatcggaagaggggaagtcgagggagggatggtaaggaggactgcatgggtcagcacaggctgccaaagccagggccagttaaagcgactccaatgcggaagagagtaggtcg
  BCprimer_F: atgataactaggTGCAAACAGGAAGTGAACGG
  BClen_F: 12
#optional
  BCprimer_R: 
  BClen_R: 
  unique_sequence: 
N26_WW1591_sg25_4kBC27_2020_06_19_853T4d_RNPsyn25_4d: 
  reference_id: 4kAAVS1
  left_150bp: tgcaaacaggaagtgaacggggaagggagggggcttctcatctgggtgcgggaaccccacatggtacctgttagacacggcaaaacccccgtcaccacccacaggtggcgcttccagtgctcagactagggaagaggttccagcccctcc
  right_150bp: aaggagacaaagtccaggaccggctggaggggctcaacatcggaagaggggaagtcgagggagggatggtaaggaggactgcatgggtcagcacaggctgccaaagccagggccagttaaagcgactccaatgcggaagagagtaggtcg
  BCprimer_F: catctcatctcgTGCAAACAGGAAGTGAACGG
  BClen_F: 12
#optional
  BCprimer_R: 
  BClen_R: 
  unique_sequence: 
```

`RefInfo.yaml`:

```yaml
4kAAVS1: 
  sequence: tgcaaacaggaagtgaacggggaagggagggggcttctcatctgggtgcgggaaccccacatggtacctgttagacacggcaaaacccccgtcaccacccacaggtggcgcttccagtgctcagactagggaagaggttccagcccctcctccttcagagccaggagtcctggcccccagcccctcctgccttaaacccagccaggtccttccaagggtcaagctcggaaaccaccccagcagatactctgcaggaacgaagccgtgggcccagggctatgcagggtggaggaaggccaccctgtgctgggacagactcaggggcctgggcgggactcccagaggggtgagacagctgcacacctgtgtgcctgggccccaggctgtcacactccagttcactgaggccccctctgcacggggccctgcagccaggggctgacacgggccaccgtttctcattcttcccttaggggtccaaaacttggggggacaaaagccgaagtccagggggtcggaggagggacttgccccaggccttgtggacactgggtgggctccgggacctgaactggagctgaggaaggagtgaagctaaactcctagatccacgggataaattaccccccaagtccctcacctctccaaagctgcccatctggaggaggcgggagggagctacgagggccaagagcatgaggtcatggaaactcgggctgtgaaggggccgcacgtgccctgggaacgggatgaactcggctcgtttatttccacccagttgtcatggcgataggggaggggggcaaggagagcaatgggcctttccctttcaaggacctgcccagtacaggcatccctgtgaaagatgcctgaggcctgggcaccagggactccagagtccaggcccaacccctccccattcaacccaggaggccaggccccagcccttccgccctcagatgaaggagtccaggcccccagcctctccccattcagacccaggggtccaggcccagccccgcctccctaagacccagaagtccaggcccccagcccctcctccctcagacccacgagtccaggccccagcccctcctccctcggacccaggagtccaggcccccagtccctccaccctcagacccaggagtccaggccccagcccctcctccctcggacccaggagtccaggccccagcccctcctctctcaaacccaggagcccaggcccccagctcttctctgttcagccctaagaatcctggctccagcccctcctactctagcccccaaccccctagccactaaggcaattggggtgcaggaatgggggcagggtaccagcctcaccaagtggttgataaacccacgtggggtaccctaagaacttgggaacagccacagcaggggggcgatgcttggggacctgcctggagaaggatgcaggacgagaaacacagccccaggtggagaaactggccgggaatcaagagtcacccagagacagtgaccaaccatccctgttttcctaggactgagggtttcagtgctaaaactaggctgtcctgggcaaacagcataagctggtcaccccacacccagacctgacccaaacccagctcccctgcttcttggccacgtaacctgagaagggaatccctcctctctgaaccccagcccaccccaatgctccaggcctcctgggataccccgaagagtgagtttgccaagcagtcaccccacagttggaggagaatccacccaaaaggcagcctggtagacagggctggggtggcctctcgtggggtccaggccaagtaggtggcctggggcctctgggggatgcaggggaagggggatgcaggggaacggggatgcaggggaacggggctcagtctgaagagcagagccaggaacccctgtagggaaggggcaggagagccaggggcatgagatggtggacgaggaagggggacagggaagcctgagcgcctctcctgggcttgccaaggactcaaacccagaagcccagagcagggccttagggaagcgggaccctgctctgggcggaggaatatgtcccagatagcactggggactctttaaggaaagaaggatggagaaagagaaagggagtagaggcggccacgacctggtgaacacctaggacgcaccattctcacaaagggagttttccacacggacacccccctcctcaccacagccctgccaggacggggctggctactggccttatctcacaggtaaaactgacgcacggaggaacaatataaattggggactagaaaggtgaagagccaaagttagaactcaggaccaacttattctgattttgtttttccaaactgcttctcctcttgggaagtgtaaggaagctgcagcaccaggatcagtgaaacgcaccagacggccgcgtcagagcagctcaggttctgggagagggtagcgcagggtggccactgagaaccgggcaggtcacgcatcccccccttccctcccaccccctgccaagctctccctcccaggatcctctctggctccatcgtaagcaaaccttagaggttctggcaaggagagagatggctccaggaaatgggggtgtgtcaccagataaggaatctgcctaacaggaggtgggggttagacccaatatcaggagactaggaaggaggaggcctaaggatggggcttttctgtcaccaatcctgtccctagtggccccactgtggggtggaggggacagataaaagtacccagaaccagagccacattaaccggccctgggaatataaggtggtcccagctcggggacacaggatccctggaggcagcaaacatgctgtcctgaagtggacataggggcccgggttggaggaagaagactagctgagctctcggacccctggaagatgccatgacagggggctggaagagctagcacagactagagaggtaaggggggtaggggagctgcccaaatgaaaggagtgagaggtgacccgaatccacaggagaacggggtgtccaggcaaagaaagcaagaggatggagaggtggctaaagccagggagacggggtactttggggttgtccagaaaaacggtgatgatgcaggcctacaagaaggggaggcgggacgcaagggagacatccgtcggagaaggccatcctaagaaacgagagatggcacaggccccagaaggagaaggaaaagggaacccagcgagtgaagacggcatggggttgggtgagggaggagagatgcccggagaggacccagacacggggaggatccgctcagaggacatcacgtggtgcagcgccgagaaggaagtgctccggaaagagcatccttgggcagcaacacagcagagagcaaggggaagagggagtggaggaagacggaacctgaaggaggcggcagggaaggatctgggccagccgtagaggtgacccaggccacaagctgcagacagaaagcggcacaggcccaggggagagaatgcaggtcagagaaagcaggacctgcctgggaaggggaaacagtgggccagaggcggcgcagaagccagtagagctcaaagtggtccggactcaggagagagacggcagcgttagagggcagagttccggcggcacagcaagggcactcgggggcgagaggagggcagcgcaaagtgacaatggccagggccaggcagatagaccagactgagctatgggagctggctcaggttcaggagagggcagggcagggaaggagacaaagtccaggaccggctggaggggctcaacatcggaagaggggaagtcgagggagggatggtaaggaggactgcatgggtcagcacaggctgccaaagccagggccagttaaagcgactccaatgcggaagagagtaggtcg
4kBCL11A37: 
  sequence: gtgtggtgttcggagtcctaagagcccccactagctcagaaatggacttagttgacctcccccattagcagcatggagagtcaaggagatgacttctaccttgccaaaggccttgggaagaaagacagcatcaaggtctcacacaacactccagggaggcagctgctgcccagtgctgtggacagcaaagcttcagtgcaggaaattaagattccccctgcctccccctcccccatcctcatcagcttggccatggcagggctgggggatcagaggtgaacaggaagcagaaggacccctgggggagacagggcctccagtgggaccagagctgagtggcctcaggcagtggcggaagctgattaaaggaaggtacggggagtggaggggaagtggacaaaagacaggacagccatcttagacaacaatgcaagggggagaaactgaagaaaacagaacagagaccactactggcaataaacagagagaaagtgaagccccatgggtgaggcacacctacattacttaagaaacctgagcacattcttacgcctagggcaataaatacatccttgagctacacaggctaagcaagagtgagagagggtgatgctgacaggccacatgggagagtgggaagacgtgggctgggagctgggagtttggcttctcatctgtgcatggcctctaaactgggcagtgaccatggcctggtcacctccccactctggacctgggttgcccctctgtaaacaaggaggttgtaataaattatctccaataccctaatgtcttataaatcttatgcaatttttgccaagatgggagtatggggagagaagagtggaaacggcccagagctcagtgagatgagatatcaaaggggacgaaaagtgttcattccatctccctaatctccaattggcaaagccagacttggggcaatacagactggttctgtgatgacaaataactcctagctcattcctaatgatttatcaccaaatgttctttcttcagctggaatttaaaatatggactcatccgtaaaataggaataataatagtatatgcttcatagggtttgtatgaaaataaaatgagtgcgtatttgtaaagttcctagagcagagtaagtgctccgagcttgtgaactaaaatgctgcctcctggtatttattagttacacctcagcagaaacaaagttatcaggccctttccccaattcctagtttgggtcagaagaaaagggaaaagggagaggaaaaaggaaaagaatatgacgtcagggggaggcaagtcagttgggaacacagatcctaacacagtagctggtacctgataggtgcctatatgtgatggatgggtggacagcccgacagatgaaaaatggacaattatgaggaggggagagtgcagacaggggaagcttcacctcctttacaattttgggagtccacacggcatggcatacaaattatttcattcccattgagaaataaaatccaattctccatcaccaagagagccttccgaaagaggcccccctgggcaaacggccaccgatggagaggtctgccagtcctcttctaccccacccacgcccccaccctaatcagaggccaaacccttcctggagcctGTGATAAAAGCAACTGTTAGcttgcactagactagcttcaaagttgtattgaccctggtgtgttatgtctaagagtagatgccatatctcttttctggcctatgttattacctgtatggactttgcactggaatcagctatctgctcttacttatgcacacctggggcatagagccagccctgtatcgcttttcagccatctcactacagataactcccaagtcctgtctagctgccttccttatcacaggaatagcacccaaggtccatcagtacctcagagtagaaccccctataaactagtctggtttgcccatggggcacagtcaggctgttttccagggtggggtgcagacattctctgcctgttgtgatgcttacatataacgtcataacagacacacgtatgtgttgtgatccctgtggtttgagagtttggagcttccctaaaagtcaaaatattctcaatgggccctcaatcagcacatacacacaaaaggtacctggaaaactgtaattcttttcctgctcaaagacaggcaattcaataccccttcccccaaccaaaaacccttgccaccatgggagcctggggcagagaaggcacagtgaagtcaaactgtaattccaggctctaaatggtgctgtcatttttctgagagtctctaaattacaagggtgttttcactattcttagctattttttaaaacacctaagaaacatactgcagctctggaaaagagaacaaacaaaccaaagagaagggatccagaggtcaccctcatatgtgaaaagtcaattgataatgaaggctttaggataaccggaggggagatgattgaaagcaatgcacctgtgcaggaaatggattacggaaacagggaattgttcatgaaatcccagaaaaccagaaccgggaaagttctggaagtcggaaaaacaaatcatgacttaagcaatggaagtccaatacacgtttacagaatgccttgtcccacgaggcaacacaggctaccacagatgggggacagggtgggagtggaccatcccagtggtgttactgaggggcaaagggatagccctatgaggcaagtgtccagggcagaactggagctttgtgaaaccatttcccaggcagagacagagcactaggctggtgctgccagtctgacaataagtctgccattgtcctctggtcagctctggacacacagcaaaagtgagttcagagtagcctgaagcaggaaagagggaagagaggaggataacacctatcttccactttgctgcaggttcaaggcaaggatttgagacagttaccccttctggaagagcctggtgagtacatctctcctgccttgtacaaccctctctcctcaccgactttctctcccagcagccagcaggggcctgggccatttatggaatgcaagccctgaccacacagacttacttacatgccaggacagccaccaggtagcctttcccactctaggttccactgtgagtgctctctctctctctctctctcactatgctcccaagaggagtcttacatcaaccccttcctcaaatctccctcactggatgtcacagtcataggcctgaaaagcagcatgcaaactgaatttttgtaaagcaggacccatttccccatggacagtcataagagatgagtgaacacaatgtagcacttaatttctgtcttcacgattacttcacgataaatctggattccaaagggactataagctctcacatggaaggaagcaagatctctactcctcccccagtgttgagtggacagggagtacaccgcagacacctgttggccaaccaattctaattccctttagctagcatcccctaagctagagctagagctagagctatttccttgcagccttccttttctctagcaaagtccttccatgcagtagctaatgacctgtaaacacttaatgagctagagaaacattccattgaaaggaataccactgtgcatccttttgtaaagaggggggaaaatcttttgtaaaacgaagcatcgcctttaactgctctgtttgatcaagtcagatttttcagaatatgaatagctagtattcaagcatatatgaactgtctttaagttaatcaatccctagaaactagccctcaggttagcaggccaaggatatatgagagtgctttgaagtctagacttaaactgccgctcct

```

### Pipeline Output

When running the full pipeline, the results of each step are outputted in a separate folder for each step. The output folders and their respective contents are as follows:

* **Reference**: Contains the FASTA file corresponding to the reference sequence input in RefInfo.
* **demultiplexed**: Contains the demultiplexed FASTQ files, and a brief statistics  file with the number of files/reads.
* **visualization**: Contains the `.sorted.bam` and `.bai` for each demultiplexed FASTQ file.











