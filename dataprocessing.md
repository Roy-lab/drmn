# Preparing input data for DRMN

The regulatory features could be context-specific (e.g., chromatin marks, accessibility) or context independent (e.g., sequence-specific motifs). For context-independent features, we use instances of motifs in promoter regions of genes.

![alt text](example_input/motif_small.png "Motif instances in gene promoter.")

For context-specific features, we aggregate regulatory signals (different histone modification, or chromatin accessibility) in promoter regions of genes.

![alt text](example_input/signal_small.png "Aggregated signals in gene promoter.")

Furthermore, if accessibility signals (e.g. DNase-seq/ATAC-seq) are available, we can aggregate chromatin accessibility signals in motif instances.

![alt text](example_input/qmotif_small.png "Q-Motif, aggregated signal in motif instances in gene promoter.")

We used our in-house softwares to aggregate signal in promoter regions or motif instances. 

### Generating Histone and Accessibility features
We use the program [aggregateSignal](https://github.com/Roy-lab/drmn_utils/tree/master/feature_generation_tools/aggregateSignalRegion_nonLog) to generate the Histone and Accessibility features. This tool will aggregate the signal (ATAC-seq, DNase-seq, or histone modification) in promoter region of genes. Suppose we have bam files for a particular histone mark or accessibility. To generate this feature the steps are: 

1. Use bedtools to convert bam files to count files:
```
bedtools genomecov -ibam input.bam -bg -pc > output.counts
```
2. Apply aggregateSignal program to calculate coverage in input regions. 

The program is used as:
```
./aggregateSignal promoters.txt mm10.fa.fai output.counts output.txt
```
where promoters.txt is a gff like file:
```
chr1	CONVERT	TSS_Capn8_ENSMUST00000192671	182562537	182567537	.	+	.	Capn8
chr2	CONVERT	TSS_Casc4_ENSMUST00000110586	121864469	121869469	.	+	.	Casc4
```
that contains the coordinates of each promoter (defined as +/- Xkb of TSS) and can be generated from ensembl gtf files. In our experiments, we used +/- 2500kb of the TSS. 

The output.txt reports the signal per promoter of each gene:
```
Capn8	13.9519
Casc4	4.59291
```

See feature_generation_tools in [DRMN utils](https://github.com/Roy-lab/drmn_utils/) for example of how this is used.


### Generating Q-Motif feature (accessibility signal of motifs)

We use the program [aggregateSignalMotifNet](https://github.com/Roy-lab/drmn_utils/tree/master/feature_generation_tools/aggregateSignalMotifNet) to aggregate signal (e.g., from ATAC-seq or DNA-seq) in motif instances mapped to promoter region of genes. To generate these features the steps are:

1. Use bedtools to convert bam files to count files:
```
bedtools genomecov -ibam input.bam -bg -pc > output.counts
```
2. Map motif instances to promoters by using either bedtools or our in house script, (matchMotifToGenePerTF2.py)[https://github.com/Roy-lab/drmn_utils/blob/master/feature_generation_tools/aggregateSignalMotifNet/matchMotifToGenePerTF2.py] to produce a file listing the set of motif instances mapped to promoters. Bedtools can also be used by intersecting genome-wide instances with the promoter coordinate file. Suppose this file is called motif_promoters.txt. The input file to aggregateSignalMotifNetis expected to have a gff format:
<motif regions> <tss list> <upstream window> <downstream window> <output>

```
chr21	CONVERT	TSS_Gene2_TSS1	35198809	35198821	0	+	0	Motif2;Gene1
chr21	CONVERT	TSS_Gene2_TSS1	35144166	35144177	0	+	0	Motif2;Gene2
```
that contains the coordinates of motif instances in promoters. This can be prepared with 

```
python matchMotifToGenePerTF2.py <motif regions> <tss list> <upstream window> <downstream window> <output>

*Ali, can you update with example unput motif instance and promoter list input files for your case?*

python matchMotifToGenePerTF2.py <motif regions> <tss list> 2500 2500 motif_promoters.txt
```
**_Note:_** The last column has both motif name and target gene name, separated by ";".

**_Note:_** In both programs, in the input gff like file, set the strand sign to +. 


3. Apply aggregateSignalMotifNet to quantify accessibilty of each motif:
```
./aggregateSignal motif_promoters.txt mm10.fa.fai output.counts output.txt
```
This tool will then be simialrly applied to the planned data set for each mark/measurement for each time-point/condition to generate the full set of aggregated Q-motif feature data to describe the data set.

See feature_generation_tools in [DRMN utils](https://github.com/Roy-lab/drmn_utils/) for example of how this is used.

### Feature merging

4. To merge the Q-motif feature data across time-points/conditions the data for each such condition should be merged columns wise. In practice this merging can of course be done flexibly in different context. Locally our mergeData tool was used here in practice https://github.com/Roy-lab/drmn_utils/tree/master/mergedata.

Here the input is a list of (tab-delimited) .txt data files to merge, columnwise. The output is a single .txt file containing the combined data values from all of the input data. This is as indicated in the readme for this tool at https://github.com/Roy-lab/drmn_utils/tree/master/mergedata. Each row of the merged output data matrix will represent one motif-gene pairing, and the columns will provide the data values for each condition from the input files.

Cases of motif sites with no coverage in one of the condition data sets can appear with the "<nodata>" symbol in one or more columns. In practice those motif-promoter pairings with incomplete data (instances of <nodata >) can be removed (assuming they represent a small minority of the overall feature set), or substituted with zero values. 

Note: We also concatenate the results for different motifs row-wise at this point, since the Q-motif feature names are unqiue for each gene promoter and motif pair, and this simplifies things for the next (normalization) step. 

### Feature normalization

For the merged Q-motif (or chromatin mark) data set we then log transform and quantile normalize the obtained values across cell lines/time points. Again, such transformations can flexibly be applied in varying ways. Locally this was done in matlab using a script, https://github.com/Roy-lab/drmn_utils/blob/master/quantile_normalize.m, with the following usage:

```
 quantile_normalize(<input file>,<normalized file>,logOption)
 
 or 
 
 quantile_normalize(merged_example.txt,normalized_example.txt,'doLog')
 
```

The output normalized data matrix will then represent the log transformed and quantile normalized Q-motif feature data values for use in DRMN.

For berevity, if starting with a Matlab data array <data> the output normalized matrix result can liekwise be obtained with the following two lines, and subsequently written out in any specifically preferred way:

``` 
log(<data>+1)
quantilenorm(<data>)
```

### Formatting

Finally, prepare the final per-condition feature data files for input to the DRMN algorithm. 

Here we add cell line/time point specific suffixes to gene names in order (so gene names will be specific to a cell line while regulators will be the same for all cell lines). 

In the feature file, the first line is number of regulators and number of genes in the file (tab-delim).
The rest of the file will be in 3 columns format (tab-delim):
```
feature g_cell  value
```
