# PEB_Phylogenomics

This hands-on course will introduce you to the awesomeness of phylogenomics.

You will learn the most fundamental steps in the phylogenomics pipeline, allowing you to go from a bunch of sequences to a phylogenetic tree that represent how the species are related to each other. 

The phylogenomics pipeline can become very complex, adding many steps (particularly when assembling the datasets!) and some analyses can take weeks to complete. Pipelines are modular, meaning they can (and should) be improved, as well as modified to the particular question at hand.

## Objective and dataset

We will use a dataset [from this paper](https://academic.oup.com/sysbio/article/65/6/1057/2281640). The starting point are a subset of proteins obtained from genomes/transcriptomes from 23 species of vertebrates and our aim is to reconstruct the phylogeny of these species using concatenated and coalescent approaches. In practice, we are using a subset of the full genomes/transcriptomes of these species, only to speed up computations.

Let's start by downloading the data.

Connect to [our server](https://datasciencehub.ifca.es/). Download the dataset from [this respository](https://github.com/iirisarri/UIMP-phylo_pipeline/conidae_mito_nuclear.zip) and decompress it into your preferred location.
```
wget https://github.com/iirisarri/UIMP-phylo_pipeline/vertebrate_proteomes.tar.gz
tar zxvf vertebrate_proteomes.tar.gz
```

## Inferring ortholog groups

We will infer ortholog groups with [Orthofinder] ().

```
orthofinder -f vertebrate_proteomes

```

The single-copy orthologs will be in a file called Orthogroups.csv. We can use mirlo to parse this file and create new fasta files, each containing orthologous sequences for every species. Do something like this:

It requires [JDK] (https://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html) to be installed.

```

mirlo.py -c Results_Feb28/Orthogroups.csv -i vertebrate_proteomes -o MIRLO_OUT

```

But if you look into Orthogroups.csv, you will see it is a simple tab-separated file. Thus, we can also parse it with a bit of bash. Tip: separate sequences for each orthogroup and place them into a individual files, one ID per line: then extract the sequences for each orthogroup using [seqtk](https://github.com/lh3/seqtk) or a [perl script](link).

```
# each line contains the sequences belonging to one orthogroup
split -l 1 Orthofinder_Results_Feb28/Orthogroups.csv

# except the first line, which contains the column headers and can be ignored
rm xaa 

# you can create files .taxa containing the sequences for each orthgroup
for f in x*; do name=`cat $f | cut -f1` ; tr '\t' '\n' < $f | tail -n+2 > $name.taxa; done 

# create a file containing all proteins from all species
cat vertebrate_proteomes/*faa > vertebrate_proteomes_all.fasta

# for each orthogroup, extract the sequences from the big fasta file using seqtk
for f in OG00000*taxa; do /Applications/Phylogeny/seqtk/seqtk subseq vertebrate_proteomes.fasta $f > $f.fas; done

# aternatively, use a perl oneliner
for f in OG00000*taxa; do perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' $f vertebrate_proteomes.fasta; done

```
At this point, we will have one file per gene, containing one ortholog per species.

Let's make taxon names homogeneous across ortholog groups; this is necessary for the concatenation step. You will see that the difference between headers is just a gene number, which we must remove.

```
for f in *taxa.fas; do sed -e '/>/ s/_GENE_.*//g' $f > out; mv out $f ; done

```

**NOTE**: Paralogy is a difficult issue to solve and using Orthofinder alone might not be enough. Paralogy is tricky business! [Research] (link) has shown that including paralogs into a phylogenomic dataset can bias the results, particularly when phylogenetic signal is weak. Paralogs should always be removed prior to phylogenetic inference, but identifying them is not easy. One could build single-gene trees and look for sequences producing extremely long branches or clustering outside of the remaining sequences. But this is not a trivial task and can become very time consuming for large datasets, but still not a good enough reason to not do a proper job!


## Pre-alignment quality filtering

Often, transcriptomes and genomes have stretches of erroneous, non-homologous amino acids or nucleotides, produced by sequencing errors, assembly errors, or errors in genome annotation. 

Until recently, these errors had been mostly ignored because [no automatic tool could deal with them](https://natureecoevocommunity.nature.com/users/54859-iker-irisarri/posts/37479-automated-removal-of-non-homologous-sequence-stretches-in-phylogenomic-datasets). We will use [PREQUAL] (https://doi.org/10.1093/bioinformatics/bty448), which takes sets of unaligned sequences and identifies sequence stretches sharing no evidence of homology, which are masked in the output.

```
for f in *fas; do /Applications/Phylogeny/prequal/prequal $f ; done

```

The filtered (masked) alignments are in .filtered whereas .prequal contains relevant information such as the number of residues filtered.


## Multiple sequence alignment

The next step is to infer multiple sequence alignments. They allow us to *know* which amino acids/ nucleotides are homologous. A tool that works very well and it is easy to use is [MAFFT](https://mafft.cbrc.jp/alignment/server/).

The strategy is to align each gene file separately, using a for loop.

```
for f in *filtered; do mafft $f > $f.mafft; done

```

## Alignment trimming

Some gene regions (e.g., fast-evolving regions) are difficult to align and thus positional homology is unceratin. It is unclear (probably problem-specific) whether trimming badly-aligned regions [improves](https://academic.oup.com/sysbio/article/56/4/564/1682121) or [worsens](https://academic.oup.com/sysbio/article/64/5/778/1685763) tree inferece. However, gently trimming very incomplete positions (e.g. with >80% gaps) will speeds up computation in the next steps without significant information loss.

To trim alignment positions we can use [BMGE](https://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-10-210) but several other software are also available.

Removing alignment positions with > 20% gaps.

```
for f in *mafft; do java -jar /Applications/Phylogeny/BMGE-1.12/BMGE.jar -i $f -t AA -g 0.2 -h 1 -w 1 -of $f.gt02; done

# By default, BMGE will trim high-entropy (likely fast-evolving) and incomplete positions with >20% gaps

for f in *mafft; do java -jar /Applications/Phylogeny/BMGE-1.12/BMGE.jar -i $f -t AA -of $f.bmge; done

```

It is always useful to check the intermediate results. Multiple sequence alignments can be visualized in [SeaView] (link) or [AliView] (link). Also, one could have a quick look using command line tools (`less -S`). In this case it is better to use the phylip format. We can easily reformat alignments with a simple script:


```
for f in *fa; do fasta2phylip.pl $f > $f.phy; done

```

## Concatenate alignment

Create a super-alignment by concatenating all gene files. We will use a custom script. Another good option is [FASconCAT](link), which will read in all \*.fas \*.phy or \*.nex files in the working directory and concatenate them (in a random order) into a super-alignment.

```
perl concat_fasta_partitions.pl *filtered.mafft.gt02 > vert_56g_filtered_g02.fa
mv partitionfile.part vert_56g_filtered_g02.part
```

Congratulations!! Our concatenated dataset is ready to rock!!

