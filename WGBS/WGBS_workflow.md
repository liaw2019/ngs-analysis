https://pmc.ncbi.nlm.nih.gov/articles/PMC8196797/#sec4-ijms-22-05426

# 1. WGBS data processing
Step 1 - TrimGalore (trim 18 bp from 3' and 5' to eliminate adaptase tail, reads <20bp were discarded)  
Step 2 - Mapping to Arabidopsis genome (TAIR10) with BSMAP  
Step 3 - Sambamba (sort, index divided by PCR duplicates)  
Step 4 - Deduplicated read were passed for methylation call    
Step 5 - Qualimap to evaluate alignment  
Step 6 - BSMAP methratio.py to extract methylation call (cutoff of 1 count)  
Step 7 - Annotate cytosine location with table browser function  
Step 8 - Methylation analysis with Bismark  
Step 9 - Deeptools to analyze and plot methylation patterns across gene bodies and 2kb upstream/ downstream (reference bed files from UCSC genome browser and ReMap Regulatory Atlas)

# 2. Differential Methylation Analysis
