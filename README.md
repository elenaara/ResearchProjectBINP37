## README - Metatranscriptomic and metagenomic binning to explore host:symbiont co-evolution in mutualistic symbioses
- ***Elena Aramendía Cotillas***
- github repository: https://github.com/elenaara/ResearchProjectBINP37

-------------------
In this project I worked both with metatranscriptomic RNA-seq samples from cultured microbial consortia including breviate and prokaryotic organisms, and with publicly available SSU rRNA sequences from Breviates and Arcobacter to study their occurrence in different environments.

**1) Metatranscriptomic data from cultured breviate-prokaryote consortia.**
Starting data were raw paired-end reads containing metatranscriptomic RNA-seq data from eight sequencing projects of a microbial consortia consisting of a breviate and prokaryotic species.
Eight different strains were analysed: Breviata anathema, BLO, FB10N2, LRM1B, LRM2N6, SaaBrev, Lenisia limosa and Pygsuia biforma. The first six species derive from in-house projects, while the L. limosa (Hamann et al., 2016) and P.biforma (Brown et al., 2013) data was previously published and reads are publicly available in NCBI under BioProject PRJNA277740 (L.limosa) and PRJNA185780 (P. biforma). For L.limosa, we had both metatranscriptomic and metagenomic data.
Raws were trimmed using FastQC and Trimmomatic. Then, on one hand, I used the Phyloflash pipeline to assemble the SSU rRNA (16/18S) sequences from the reads. This allows us to assign a taxonomy to the different SSU sequences found and thus identify the different species found in the sample.
In parallel, I used the Trinity software to get assembled transcripts from the reads.

**2) ENVIRONMENTAL DATA.**
For the public SSU rRNA sequences from Breviates and Arcobacter, I searched for matches of these sequences in the SRA archive using the IMNGS platform.This platform outputs matches to your query sequences found in publicly available SRA experiments. The SRA sequences obtained from IMNGS were further filtered using SILVAngs, which asigns a taxonomy to the sequences according to the SILVA database. Then I used the information on NCBI for each SRA number to get the environment from which the data was collected, and studied the environment were Breviates are found, and if Breviate and Arcobacter sequences co-occurred in each experiment.

-----------
### SOFTWARE USED
- FastQC v0.11.9
- Trimmomatic v0.39
- Phyloflash v3.4
- Trinity v2.14.0
- TransDecoder v5.5.0
- DIAMOND v2.0.13
- MEGAN Community Edition v6.22.2
- SILVAngs v1.4
- cd-hit v4.8.1
- MAFFT v7.310
- JalView v2.11.2.2
- iqtree v2.0.3

Most steps were performed in the Kebnekaise cluster at HPC2N.

----------------
### 1 - METATRANSCRIPTOMIC CULTURED DATA
#### 1.1 - CLEANING THE READS - FastQC and Trimmomatic
Each file containing metatranscriptomic reads in FASTQ format was analyzed using FastQC, and then filtered using Trimmomatic with the following settings:

- **ILLUMINACLIP:** Remove Illumina adapters provided in the TruSeq2-PE.fa file (provided with the software). Initially Trimmomatic will look for seed matches allowing maximum 2 mismatches. These seeds will be extended and clipped if in the case of paired end reads a score of 30 is reached.
- **HEADCROP:** first 8 positions will be removed.
- **SLIDINGWINDOW:** scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 20.
- **MINLEN:** drop reads which are less than 36 bases long after these steps.

```bash
java -jar trimmomatic-0.39.jar PE -threads 4 R1.fastq R2.fastq \
1_trim2.fastq 1_unpaired_trim2.fastq 2_trim2.fastq 2_unpaired_trim2.fastq \
ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 \
HEADCROP:8 SLIDINGWINDOW:4:20 MINLEN:36
```

#### 1.2 - SSU rRNA ASSEMBLY - PhyloFlash
PhyloFlash pipeline assembled 16/18S sequences from the reads by mapping them to the SILVA 138.1 database. The software was run using the clean reads from the previous step, with the default settings except for cluster id, set to 90%. The assembled 16/18S sequences from the different samples were all collected into one file for posterior analysis. The software also assigns taxonomy to each assembled sequence according to its closest match in the SILVA database.

```bash
phyloFlash.pl -poscov -CPUs 8 -read1 $reads_F -read2 $reads_R -lib $libname -dbhome $dbhome -clusterid 90
```

#### 1.3 - TRANSCRIPT ASSEMBLY - Trinity
Trinity is a de novo transcript assembler designed for RNA-seq reads. It was run for each of the samples with the default settings.

```bash
singularity exec -e Trinity.simg  Trinity \
         --seqType fq \
         --left $reads_F  \
         --right $reads_R \
         --max_memory 50G --CPU 8 \
         --output $outfolder
```

#### 1.4 - PREDICTION OF CODING SEQUENCES - TransDecoder
TransDecoder identifies candidate coding regions within transcript sequences. Transcripts assembled by Trinity were used as input for this software.

Two steps were run;
- **Step 1:** identify  long open reading frames (ORFs). By default, identifies ORFs that are at least 100 amino acids long.
- **Step 2:** identify likely coding regions.

It outputs both nucleotide and peptide sequences for the candidate ORFs, a gff3 file with positions within the target transcripts of the final selected ORFs and a bed file describing ORF positions.

```bash
TransDecoder.LongOrfs -t $fasta_file
TransDecoder.Predict -t $fasta_file
```

#### 1.5 - TAXONOMIC ASSIGNMENT OF THE TRANSCRIPTS - DIAMOND and MEGAN code
Taxonomic assignment of the predicted transcripts was performed using DIAMOND aligner v2.0.13 and the MEGAN Community Edition (CE) software (v6.22.2). First, sequences are aligned to the nr database using DIAMOND and then the alignment file is “meganized”; the MEGAN command performs taxonomic and functional classification of the sequences and appends the resulting information to the file so that it can be directly opened in MEGAN CE. In the MEGAN GUI taxonomic assignments can be visualized and sequences assigned to a particular taxonomy can be extracted.

```bash

# DIAMOND
diamond blastx -d $diamond_nr \
      -q $inputfile \
      -o $prefix_for_outputs.daa \
      -f 100 \
      -F 15 --range-culling --top 10

# daa-meganizer ( Prepares ('meganizes') a DIAMOND .daa file for use with MEGAN )
daa-meganizer -i 2_daa_trinity/${prefix_for_outputs}.daa -mdb $megan_mapping

```
### 2 - ENVIRONMENTAL DATA
#### 2.1 - RETRIEVAL OF SEQUENCES - IMNGS.
The IMNGS platforms allows to input query sequences and search for matches of these sequences in the SRA database. SSU sequences from both Breviates and Arcobacter were used as query to retrieve all sequences from SRA experiments matching these organisms, and thus explore their occurrence in different environments, if they appear together or not, and in what proportion.

- **For the Breviates**, 17 different sequences were used as query: 10 sequences obtained from NCBI (Accession numbers: KT023596.1, KC433554.1, HQ342676.1, GU001164.1, AF153206.1, HF568854.1, HM103503.1, HM103491.1, HM103465.1) + 7 in-house sequences (PCE strain clone 1 and 2, FB10N2 strain clone 1 and 2, LRM1b strain, LRM2N6 strain). The search was conducted using a 90% identity threshold.

- **For Arcobacter**, we used 9 sequences obtained from NCBI (Accession Numbers: NR_174232.1, NR_136421.1, NR_117919.1, NR_117570.1, NR_117569.1, NR_117105.1, NR_116729.1, NR_102873.1, NR_025906.1). The search was conducted using a 97% identity threshold.

#### 2.3 - TAXONOMIC ASSIGNMENT OF SRA SEQUENCES - SILVAngs
To confirm the results of the IMNGS platform, we used the retrieved sequences as input for the SILVAngs platform. This software uses the SILVA database to perform taxonomic assignment to the query sequences. We performed this step mainly to filter the outputof the Breviate search, although the filtering was done on both outputs (Breviate and Arcobacter). Since we used a lower threshold of 90% identity in the Breviate search to increase sensitivity we expect to get many sequences belonging to different organisms. SILVAngs was run using the default settings.

Once the sequences were assigned a taxonomy according to the SILVA database, we kept the following sequences:

- **For the Breviates:** since *Breviata anthema* is classified as *Eukaryota;Amorphea;Amoebozoa;Incertae Sedis;Breviatea;Breviata;* in the SILVA database, we kept all Eukaryotic sequences annotated
as "*Incertae Sedis*" in at least one of the 5 first taxonomic levels.
- **For the Arcobacter:** we kept sequences annotated as family Arcobacteraceae.

The SILVAngs software clusters the sequences in the process, and outputs only the reference to each OTU. A custom python script was used to get all sequences asigned to the taxa of interest from the output fasta and the cluster mapping files. The script (get_clustered.py) is available on the github repository of the project.

#### 2.4 - PRESENCE IN DIFFERENT ENVIRONMENTS.
To study the different environments in which these organisms are found, the information of each SRA experiment was extracted from NCBI, and the environments of each metagenome were extracted selecting only those annotated as "_ metagenome".
Percentage of Breviate sequences found in each environment was plotted using Python package matplotlib. The script (imngs_compare.py) is available on the github repository of the project.

### 3 - PHYLOGENETIC DATASET CONTRUCTION AND INFERENCE
To build the species tree of the Arcobacteraceae, I constructed a query dataset composed of reference Arcobacter sequences (16S sequence used in the IMNGS search, accession numbers: NR_174232.1, NR_136421.1, NR_117919.1, NR_117570.1, NR_117569.1, NR_117105.1, NR_116729.1, NR_102873.1 and NR_025906.1 and 13 additional Arcobacter sequences from NCBI, accession numbers: FJ717164, FJ717168, FN773287, JQ862032, KF799201, JX570609, DQ357825, NR_116730, FJ717100, FJ203140.1, DQ917897, KF179644, FJ202783.1). This dataset was used as a query for a BLASTN search against the 16S rRNA database on NCBI with default settings. The top 100 hits from each search were collected and all these hits (removing duplicates) were used as reference (E-values were all 0). Sequences from organisms of genus Campylobacter were used as outgroup.  This reference dataset included 156 sequences.

The 16S sequences reconstructed from the metatranscriptomes (with Phyloflash) and SRA sequencing projects (retrieved with IMNGS) were totalled 3204 sequences. This dataset was clustered at 97% identity using the cd-hit software (v4.8.1) to reduce the total dataset to 798 sequences.
```bash
cd-hit -i imngs_associated_arco.fasta -o imngs_associated_arco_cd97.fasta -c 0.97
```

For the Breviate tree, similarly to Arcobacter both Phyloflash sequences and SRA sequences retrieved with IMNGS were included in the tree, and 7 closely related Apusomonadiae sequences (the closet phylum to Breviatea) were included as an outgroup.

Alignments and tree construction were performed in the same way for both groups.

Alignments were constructed using MAFFT L-INS-i.
```bash
linsi all_arcobacter_taxa.fasta > all_arcobacter_taxa.maL.fasta
```

Trim alignment using TrimAL:
```bash
trimal -in all_arcobacter_taxa.maL.fasta -out all_arcobacter_taxa.maL.g9cons60.fasta  -gt 0.9 -cons 60
```
This settings remove all positions in the alignment with gaps in 10% or more of the sequences, unless this leaves less than 60% of original alignment. In such case, print the 60% best (with less gaps) positions.

Trees were constructed using iqtree, with 1000 ultra fast boostraps and ModelFinder to perform model selection.

```bash
iqtree -s  all_arcobacter_taxa.maL.g9cons60.fasta -m MFP -B 1000
```

Trees were visualized the trees using figtree.
