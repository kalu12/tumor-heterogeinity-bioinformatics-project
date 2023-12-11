# tumor_heterogeinity_bioinformatics_project

This task involves working with two files:

cnv.txt: Columns of interest are chromosome, start, end, cn, cn1, and cn2. cn1 and cn2 are not always available; cn = cn1 + cn2. Extract regions where cn=2 but cn1=1 and cn2=1 (or not available).

.vcf: Find all variants within the VCF file that fall into the regions extracted from the CNV file. Calculate Variant Allele Frequency (VAF) for each variant. Remove variants with VAF values exceeding 60.

Use one of the 1D cluster algorithms (HDBScan, KDE, etc.) to group the remaining variants into clusters.

The output of the analysis should include:

Graph of clustered variants.
Positions of cluster centers.

This graph shows centered clusters as final output.
![image](https://github.com/Jovan53/tumor_heterogeinity_bioinformatics_project/assets/152201867/677a5ae0-9611-4aef-8acd-c6b1d3c04107)

