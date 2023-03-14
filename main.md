Analysis Workflow for PAC-seq data
===================================

Introduction
------------

Poly(A) tail sequencing (PAC-seq) is a powerful method for genome-wide profiling of mRNA 3' ends. In this guide, we will illustrate the workflow for analyzing PAC-seq data in Homo sapian, which includes three replicates with reads containing a poly(A) tail.

1.Data Processing
---------------

The first step is to identify reads with valid poly(A) tail and trim the tail. We will use the `MAP_findTailAT.pl` script, which is a Perl script developed in-house. The script can be executed with the following command:

bash

```Bash
perl MAP_findTailAT.pl -in ./input.fastq -poly T or A -ml 25 -mp 6 -mg 5 -mm 2 -mr 2  -mtail 6 -debug F -bar 8 -odir "/output-path/" -suf "sample_name"
```

Here is a brief description of the parameters:

*   `-in`: input FASTA or FASTQ file
*   `-poly`: sequence type of poly(A) tail. Use `A` if the sequence has As tail, `T` if it has Ts tail
*   `-ml`: minimum length after poly(A) trimming
*   `-mp`: minimum length of successive poly(A) (default=8)
*   `-mg`: margin from the start (poly=T) or to the end (poly=A) (default=5)
*   `-mm`: number of mismatches between TNNTTT or ANNAAA (default=2)
*   `-mr`: minimum number of Ts in the tail (default=3)
*   `-mtail`: minimum length of trimmed tail (default=8)
*   `-bar`: length of barcode if there is one
*   `-odir`: output directory (default is the same as input)
*   `-suf`: output filename suffix (default: xx.suf.T/A.fq)
*   `-debug`: debug mode (default: T)

After poly(A) tail trimming, we will use Trimmomatic to remove adapter sequences and low-quality reads. Here is the command:

bash

```bash
java –jar  /path/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 16  -phred 33  ./sample_name.A.fq  ./sample_name.fq /path/ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:25
```

2.Sequence Alignment with STAR
----------------------------

We will use STAR, one of the most popular tools for sequence alignment, to map the clean PAC-seq data to the reference genome. First, we need to generate genome indices with the following command:

bash

```bash
STAR --runMode genomeGenerate \
   --runThreadN 16  \
   --genomeFastaFiles ./ref.fa \
   --sjdbGTFfile ./ref_annotation.gff3 \
   --sjdbGTFtagExonParentTranscript Parent \
   --genomeDir    ./index
```

Then we can map the clean PAC-seq data to the reference genome using the following command:

```bash
STAR --runThreadN 20 \
     --readFilesIn ./sample_name-rep1.fq
     --genomeDir ./index \
     --outFileNamePrefix   ./sample_name-rep1\
     --outMultimapperOrder Random \
     --outFilterMultimapNmax 1
```

3.Identification of Polyadenylation Sites
---------------------------------------

`PAT2PA2PAC.sh` is a shell script that identifies poly(A) sites using in-house Perl scripts, bedtools, and bedmap. The script has two main steps:

1.  Identifying poly(A) sites by filtering out internal-priming reads.
2.  Grouping poly(A) tags into poly(A) site clusters (PAC).

Usage:



```
bash  /PAT2PA2PAC.sh  genome   /input-path/ /output-path/   distance  input-file
```

For example,



```
bash  /PAT2PA2PAC.sh  ref.fa /input-sam-file-path/  /output-path/  24  "*.out.sam"
```

or



```
bash  /PAT2PA2PAC.sh  ref.fa /input-sam-file-path/  /output-path/  24  "sample_name1.out.sam  sample_name2.out.sam   sample_name.3out.sam "
```

### Output

The above workflow produces multiple output files:

a) `all.PAC.header` – the column name of the output file `all.PAC.PATcount`.

*   `chr`: the name of the chromosome or scaffold
*   `UPA_start`: the start position of the PAC, with sequence numbering starting at 1
*   `UPA_end`: the end position of the PAC, with sequence numbering starting at 1
*   `strand`: the direction of PAC, defined as + (forward) or – (reverse)
*   `PAnum`: the number of poly(A) sites in the PAC region
*   `tot_tagum`: the total expression level of all samples in the PAC region (count)
*   `coord`: represents reference poly(A) sites in the PAC region
*   `refPAnum`: the number of poly(A) site in this "coord"
*   `sampleName`: the PAC expression level of each sample (count)

Example:

```
cat all.PAC.header
chr	UPA_start	UPA_end	strand	PAnum	tot_tagnum	coord	refPAnum sample_name1.out.sam  sample_name2.out.sam   sample_name.3out.sam 
```

b) `all.PAC.PATcount` – the PAC result, and the file "all.PAC.header" provides the corresponding column names.

Example:



```
head all.PAC.PATcount
1	227	 227	+	1	1	227	1	0	0	1
1	5743	 5750	+	4	7	5744	3	1	6	0
1	5850	 5852	+	3	4	5850	2	1	1	2
1	5888	 5916	+	5	14	5895	8	1	9	4
```

c) `all.PAC.info` – the part of "all.PAC.PATcount" with the first eight columns, including "chr", "UPA\_start", "UPA\_end", "strand", "PAnum", "tot\_tagnum", "coord", "refPAnum".

Example:

```
head all.PAC.info
1	227	227	+	1	1	227	1
1	5743	5750	+	4	7	5744	3
1	5850	5852	+
```
 
 
