# Perform fast read backed phasing of heterozygous SNVs for 100k of indiviuals
This repository contains C++ scripts tailored for analyses with phased data, specifically involving variant call format (VCF) files.

### Installation
1. Ensure you have the necessary cpp libraries. Check the provided `Dockerfile` for the complete list.
2. Install [BCFtools](https://samtools.github.io/bcftools/howtos/install.html).
3. To compile the scripts, navigate to the repository's directory and run:
```
make
```

### VTBD

**Step 1: Create a file of heterozygous SNVs**. 
Generate a file containing phased sites per gene for each sample. Filter to allele frequency less than 1% using BCFtools.

```
bcftools view trio.vcf --max-af 0.01 -Ou | bcftools query -i'GT="alt"' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT\n]' | gzip > trio.phased_sites.txt.gz
```





