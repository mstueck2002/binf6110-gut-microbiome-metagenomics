# BINF6110 Assignment 3: Using Bioinformatic Research Tools to Conduct Taxonomic Profiling of the Human Gut Microbiome Using Shotgun Metagenomics
## Introduction

The human gut microbiome comprises diverse communities of microorganisms, including bacteria, archaea, protozoa, and fungi, that together form a complex and dynamic ecosystem. [1] These microbial communities play essential roles in host immune responses, metabolic function, nutrient dynamics, gut physiology, and the prevention and management of disease. [2] Gut microbiome organization, heterogeneity and function are strongly influenced by diet. [3] Investigating the effects of dietary variation on the microbiome is critical for guiding evidence-based nutritional decisions that enhance metabolic and intestinal health and mitigate the development of diet-related diseases. [3]  The advent of whole-genome shotgun (WGS) sequencing technologies has facilitated high-throughput metagenomic analysis of gut microbiota. [4] WGS enables effective characterization of community structures, phylogenetic diversity, intraspecies variation, functional pathway profiling, and near-complete genome reconstruction. [4] Here, we present a computational framework for shotgun metagenomic analysis of gut microbiome samples from healthy Italian adults, used to investigate how habitual diet shapes microbial abundance, diversity, and composition.

Well-founded metagenomic analysis relies on accurate assignment of taxonomic labels to sequencing reads. [5] Kraken and its optimized pipeline, Kraken2, use k-mer-based classification strategies to map genomic sequences with the lowest common ancestor taxa. [5] Kraken2 overcomes the memory-intensive limitations of Kraken by leveraging intermediate probabilistic hash tables to expedite taxonomic classification without necessitating significant RAM requirements. [5] Compared with other classification tools, including Centrifuge, CLARK, Kraken1, KrakenUniq, and Kaiju, Kraken2 achieves comparable or superior accuracy and genus-level agreement metrics while offering significantly faster processing speeds [5]. However, because Kraken2 reports classifications at the lowest common ancestor (LCA), reads from highly populated taxonomic clades with low genomic diversity are often assigned to higher taxonomic levels (e.g., genus or family), leading to underestimation of true species-level abundance. [6] To improve estimation of species- and strain-level abundance, Lu et al. (2025) developed Bracken (Bayesian Reestimation of Abundance after Classification with KrakEN), which uses a probabilistic framework to redistribute reads based on their taxonomic assignments [6]. By reassigning reads classified by Kraken, Bracken generates more accurate species-level abundance estimates with minimal false positives. [6] 

Following the estimation of species-level abundances, microbiome analyses assess differences in microbial composition within and between samples using comparative metrics known as alpha and beta diversity. [7] Alpha diversity characterizes within-sample community complexity through measures of richness, evenness (or dominance), phylogenetic diversity, and information-based indices, collectively reflecting the diversity and distribution of taxa [8]. In contrast, beta diversity quantifies differences in microbial composition across microbiomes, with metrics such as Bray-Curtis dissimilarity and the Jaccard index capturing variation in taxa abundance and presence/absence [8]. *Phyloseq* is an R-based software package for microbiome representation and analysis, providing an object-oriented framework for the efficient organization, management, and preprocessing of phylogenetic sequencing data, thereby facilitating downstream microbiome analyses. [9] The phyloseq estimate_richness() function supports various alpha diversity metrics, while the distance() and ordinate() functions can be used to compute and evaluate pairwise sample distances for beta diversity analysis. [9] *Phyloseq* also integrates with R packages including vegan, DESeq2, and ggplot2, offering a unified and reproducible framework for microbiome data analysis, and representing one of the few Bioconductor tools that supports phylogenetic trees and taxonomic clustering outputs [9].

Identifying differentially abundant taxa across microbiome samples can reveal microbial patterns associated with habitual diet, and their potential links to human health and gut physiology. Differential abundance analysis of microbiome data faces unique challenges, namely that metagenomic count data often contains a significant portion of zeros, and that microbiome datasets are compositional, meaning they do not represent absolute abundances of microbes, but rather their relative proportions within a sample. [10, 11] ANCOM-BC (Analysis of Compositions of Microbiomes with Bias Correction) is a statistical method designed to estimate absolute taxon abundances by accounting for unknown sampling fractions and correcting biases across samples. Compared with other differential abundance approaches, ANCOM-BC demonstrates strong control of false discovery rates, provides taxon-specific p-values and confidence intervals, and is well-suited for analyses with moderate to large sample sizes (n > 10). [10] 

This study evaluates a computational pipeline for the analysis of shotgun metagenomic gut microbiome samples to taxonomically classify and compare microbial diversity, abundance, and variation between individuals with differing habitual diets. 

## Methods 

### Sequencing Data Acquisition 

Metagenomic analysis was conducted on publicly available gut microbiome samples from 74 healthy Italian adults generated by the University of Naples Federico II (BioProject: PRJNA421881; SRA study: SRP126540), and sequencing was performed using the Illumina NextSeq 500 platform. Raw SRR accession files were retrieved from the NCBI Sequence Read Archive using the SRA Toolkit and converted to FASTQ format. Sequencing quality was assessed for each sample using FastQC, which reports metrics such as per-base sequence quality, GC content, sequence duplication levels, and adapter contamination. FastQC outputs were aggregated using MultiQC (BioConda), facilitating the interpretation and visualization of overall sample quality. [12]


### Taxonomic Classification and Species-Level Abundance Estimation

Assignment of taxonomic labels to sequencing reads was performed using Kraken2 in conjunction with the standard Kraken2/Bracken RefSeq database. Paired-end FASTQ files for each sample were processed to generate both classification output files, which record read-level taxonomic assignments, and report files summarizing taxonomic abundance across hierarchical levels. A confidence threshold of 0.15 was applied to balance classification sensitivity and accuracy. The --quick flag was enabled to improve computational efficiency by terminating k-mer searches after the first database match, thereby reducing runtime across all 74 samples. 

The resulting Kraken2 report files were subsequently processed using Bracken to refine species-level abundance estimates. Bracken was run using the same reference database, with parameters -l S to specify species-level classification and -r 150 to match the read length of the sequencing data. The 74 species-level report files generated by Bracken were then converted into a BIOM-format table using kraken-biom, facilitating downstream diversity and compositional analyses. 

For downstream metagenomic analysis, the BIOM-format file was imported into R using the biomformat package and converted into a *phyloseq* object using the import_biom() function. The resulting *phyloseq* object was then integrated with sample metadata, enabling the analysis of taxonomic abundance patterns stratified by diet. Rarefaction curves were generated to assess whether sequencing depth was sufficient to capture species richness. To evaluate overall taxonomic composition, sample counts were transformed to relative abundances with the transform_sample_counts() function, and the ten most abundant phyla were visualized. Samples were subsequently grouped by diet to compare relative phylum-level abundance across habitual diets.

### Alpha and Beta Diversity Analysis 

Alpha diversity metrics were selected to comprehensively capture microbial richness, diversity, entropy, and community dominance [13]. Species richness and diversity were quantified using the Chao1, ACE, Shannon, and Simpson indices, calculated from the *Phyloseq* object using the estimate_richness() function. Dominance was assessed using the Berger–Parker index (DBP), computed with the *microbiome* R package. Mean values for each alpha diversity metric were calculated across dietary groups. Phylogenetic diversity metrics were not calculated due to the absence of sequence-based phylogenetic information in the Kraken2/Bracken-derived dataset. External reference trees were not incorporated to avoid potential mismatches and the introduction of associated biases.

Beta diversity was assessed in R using the *Phyloseq* package to evaluate differences in microbial composition across habitual diets. To reduce noise from rare taxa, features with a total abundance ≤10 across all samples were removed before analysis. Community dissimilarity was quantified using both Bray-Curtis and Jaccard distance metrics. Bray-Curtis distances were calculated from abundance data to capture differences in taxon relative abundances, whereas Jaccard distances were computed from presence-absence transformed data to assess differences in community membership. Principal coordinates analysis (PCoA) was performed on each distance matrix, and ordinations were visualized with samples coloured by diet group. Statistical significance of differences in microbial community composition across dietary groups was evaluated using permutational multivariate analysis of variance (PERMANOVA), implemented in R using the adonis2() function. Separate PERMANOVA tests were conducted for Bray-Curtis and Jaccard distance matrices.

### Differential Abundance Analysis

Differential abundance analysis was performed to identify taxa associated with habitual diet using ANCOM-BC2 (Analysis of Compositions of Microbiomes with Bias Correction), which accounts for the compositional structure of microbiome datasets and corrects for biases in sequencing data while controlling false discovery rates [10]. Analyses were conducted in R using the *ANCOMBC* package on the phyloseq object. Although abundance estimates were initially generated at the species level, taxa were aggregated at the genus level (Rank6) for differential analysis to improve robustness and reduce sparsity. The model included diet as a fixed effect (fix_formula = "Diet"), and no random effects were specified. Multiple testing correction was applied using the Holm method to attempt to tightly control against false positives common in DA methods. [14] Structural zeros were identified and incorporated into the model, and sensitivity analysis with pseudo-count addition was enabled. Low-prevalence and low-depth features were filtered using a prevalence cutoff of 0.1 and a library size cutoff of 1000 reads. Additional parameters included s0_perc = 0.05 to stabilize variance estimates and neg_lb = TRUE to ensure conservative inference based on lower-bound estimates. Differentially abundant taxa were identified for each dietary contrast by filtering results for Holm-adjusted p-values < 0.05.

## Results 

### Sequencing Data Quality and Characteristics

Analysis of publicly available gut microbiome samples from 74 healthy Italian adults, sequenced using an Illumina NextSeq 500, yielded 74 paired-end datasets totalling 304.14 Gb of data. Across all samples, FastQC analysis indicated a mean sequencing depth of 34.65 million reads per sample, with an average read length of 149.35 base pairs and a mean GC content of 50.15%. Per-base sequence quality was consistently high, with mean Phred scores ≥ 27.11 across all positions, and no detectable adapter contamination was observed. Expected Illumina WGS-specific biases were observed in analysis results, including elevated per-base sequence content and per-sequence GC content, consistent with Illumina library prep coverage biases. [15]

### Global Taxonomic Abundance and Diet-Specific Taxa 

Rarefaction curves approached asymptotic plateaus across all samples, suggesting that sequencing depth was sufficient to capture the majority of species richness. (Figure 1) 

<p align="center">
  <img src="figures/F1.rcurve.png" width="600">
</p>

<p align="center">
  <small>
    <b>Figure 1. Rarefaction Curves. </b> The Vegan package function rarecurve() was used to produce a rarefaction curve (black) for each of the 74 samples. All sample curves approached asymptotic plateaus, indicating that sequencing depth was sufficient to capture the majority of species richness. 
  </small>
</p>

Taxonomic profiles across all samples revealed that the gut microbiome was predominantly composed of ten major phyla, including *Actinomycota*, *Bacillota*, *Bacteroidota*, *Chordata*, *Fusobacteriota*, *Methanobacteriota*, *Pseudomonadota*, *Thermodesulfobacteriota*, *Uroviricota*, and *Verrucomicrobiota*.

Across all samples, *Bacteroidota* and *Bacillota* were the most abundant phyla, collectively accounting for the majority of the microbial community (Figure 2). Considerable inter-individual variability was observed in the relative proportions of these dominant taxa, while lower-abundance phyla were present at consistently smaller levels.

<p align="center">
  <img src="figures/F2.phylum_abundance.png" width="600">
</p>

<p align="center">
  <small>
    <b>Figure 2. Relative Abundance of the Ten Most Abundant Phyla Across Samples. </b> Sample counts were transformed to relative abundances with the transform_sample_counts() function, and the ten most abundant phyla were visualized. Phyla are coloured according to the legend on the right. <i>Actinomycota</i>, <i>Bacillota</i>, <i>Bacteroidota</i>, <i>Chordata</i>, <i>Fusobacteriota</i>, <i>Methanobacteriota</i>, <i>Pseudomonadota</i>, <i>Thermodesulphobacteriota</i>, <i>Uroviricoda</i>, and <i>Verrucomicrobiota</i> are the 10 most abundant phyla, with <i>Bacteroidota</i> and <i>Bacillota</i> dominating relative abundances across samples. 
  </small>
</p>

When stratified by diet, *Bacteroidota* and *Bacillota* remained the predominant phyla across all dietary groups (Figure 3). Mean relative abundance profiles indicated that individuals with an omnivorous habitual diet exhibited lower proportions of *Actinomycota*, *Pseudomonadota*, and *Verrucomicrobiota* compared to both vegan and vegetarian groups. The relative contributions of minor phyla remained low and broadly consistent across diets.

<p align="center">
  <img src="figures/F3.phylum_diet_abundance.png" width="600">
</p>

<p align="center">
  <small>
    <b>Figure 3. Mean Relative Abundance of the Ten Most Abundant Phyla by Diet. </b> Sample counts were transformed to relative abundances with the transform_sample_counts() function, stratified by diet using sample metadata, and the ten most abundant phyla were visualized. Phyla are coloured according to the legend on the right. <i>Actinomycota</i>, <i>Bacillota</i>, <i>Bacteroidota</i>, <i>Chordata</i>, <i>Fusobacteriota</i>, <i>Methanobacteriota</i>, <i>Pseudomonadota</i>, <i>Thermodesulphobacteriota</i>, <i>Uroviricoda</i>, and <i>Verrucomicrobiota</i> are the 10 most abundant phyla, with <i>Bacteroidota</i> and <i>Bacillota</i> dominating relative abundances across diet groups. Mean relative abundance profiles demonstrate omnivorous individuals exhibited lower proportions of <i>Actinomycota</i>, <i>Pseudomonadota</i>, and <i>Verrucomicrobiota</i> compared to both vegan and vegetarian groups.
  </small>
</p>

### Microbial Diversity and Composition

Alpha diversity metrics were used to assess within-sample richness and diversity across dietary groups (Figure 4). Richness estimates, as measured by Chao1 and ACE, were comparable across diets, with omnivores exhibiting the highest mean richness (Chao1 mean: 4299.6), followed by vegetarians (4168.8) and vegans (3700.9). Substantial variability was observed within each group, with overlapping distributions across all diets. (Figure 4) Diversity metrics incorporating both richness and evenness (Shannon and Simpson Information indices) showed similar patterns. Omnivorous samples exhibited the highest mean Shannon diversity (3.42), while vegan (3.21) and vegetarian (3.22) groups showed slightly lower values. Simpson diversity followed a similar trend, with omnivores displaying the highest mean (0.878), compared to vegans (0.847) and vegetarians (0.837). Dominance, assessed using the Berger-Parker index, was slightly higher in vegan and vegetarian groups (0.30 and 0.30, respectively) compared to omnivores (0.27), indicating moderately lower taxonomic diversity in these groups. Overall, while minor differences in alpha diversity were observed between dietary groups, substantial overlap in distributions suggests that within-sample diversity is broadly comparable across diets, with vegans exhibiting slightly lower species richness on average.

<p align="center">
  <img src="figures/F4.alpha_diversity.png" width="600">
</p>

<p align="center">
  <small>
    <b>Figure 4. Alpha Diverity Measures Across Diets. </b> The Phyloseq package plot_richness() function was used to generate plots for Chao1, ACE, Shannon and Simpson alpha diversity measures across diet groups. Boxplots are used to demonstrate the distribution of alpha diversity measures for the samples in each dietary group, with average measures indicated by the bold horizontal lines. Samples are coloured by diet: omnivores in red, vegans in green, and vegetarians in blue. Minor differences in alpha diversity are exhibited between dietary groups; substantial distribution overlaps suggest that within-sample diversity is broadly comparable across diets. 
  </small>
</p>

Beta diversity analysis using Bray-Curtis and Jaccard distance metrics revealed no significant differences in microbial community composition between dietary groups. Principal coordinates analysis (PCoA) showed no clear clustering by diet (Figure 5). PERMANOVA results indicated that diet explained only a small proportion of the observed variation (Bray-Curtis: $R^2$ = 0.022, p = 0.743; Jaccard: $R^2$ = 0.028, p = 0.305). Together, these results indicate that neither taxonomic composition nor community membership differed substantially across dietary groups.

<p align="center">
  <img src="figures/F5.PCoA.png" width="800">
</p>

<p align="center">
  <small>
    <b>Figure 5. Beta Diversity Principal coordinates analysis. </b> Principal coordinates analysis of ordinated Bray-Curtis and Jaccard distance matrices. Samples are coloured by diet: omnivores in red, vegans in green, and vegetarians in blue. Results show no clear clustering by diet, indicating no significant differences in microbial community composition between dietary groups
  </small>
</p>

### Dietary Differential Microbial Abundance 

Genus-level differential abundance analysis using ANCOM-BC2, with diet modelled as a fixed effect and Holm correction applied, identified a limited number of significant associations. Endlipuvirus was differentially abundant between omnivorous and vegan groups (q = 0.0155), while both Endlipuvirus and Birpovirus differed between omnivorous and vegetarian groups (q = 0.0003 and q = 0.0009). (Figure 8)

<p align="center">
  <img src="figures/F8.volcanos.png" width="800">
</p>

<p align="center">
  <small>
    <b>Figure 6. Beta Diversity Principal coordinates analysis. </b> Differential abundance was assessed using ANCOM-BC2 with diet modelled as a fixed effect and Holm-adjusted p-values. The left panel shows vegans compared to omnivores, while the right panel shows vegetarians compared to omnivores. The x-axis represents log₂ fold change, and the y-axis represents −log₁₀(p-value). Genera are coloured according to significance (q ≤ 0.05), with significant taxa shown in green (vegan comparison) or blue (vegetarian comparison), and non-significant taxa in gray. <i>Endlipuvirus</i> was the only genus differentially abundant between omnivorous and vegan groups, while both <i>Endlipuvirus</i> and <i>Birpovirus</i> were differentially abundant between omnivorous and vegetarian groups.
    </small>
</p>

