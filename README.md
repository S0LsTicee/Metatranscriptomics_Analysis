# Metatranscriptomics Analysis
## Getting Started
To begin, install Linux on your system. The following tutorial will use Ubuntu 24.04 version. You can find out how to install [here](https://ubuntu.com/download/desktop). The two key software tools, KneadData and MetaPhlAn, are Python-based tools, so it's imperative to install the latest version of Python and pip to facilitate future installations.  
```ruby
sudo apt install python3 python3-pip
sudo apt install python3-pip
```

**The pipeline for metageonomics data processing is as follows:  
_SRA Download -> Convert into fastq files -> KneadData (QC + host removal) -> MetaPhlAn (taxanomic profiling)_**

KneadData is software that coordinates several external tools and is designed for quality control, creating a controlled environment that correctly links all wrappers and executables without breaking other necessary dependencies. The common practice is to use a Conda environment by installing Conda and creating a virtual environment for the analysis. Install the SRA Toolkit as well to set up for file conversion and KneadData and MetaPhlAn inside the conda environment. 
```ruby
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda create --name env_name

conda activate env_name

conda install -c bioconda sra-tools -y
pip install kneaddata
conda install -c bioconda metaphlan -y
```

KneadData was installed with Trimmomatics and BowTie2, which will be used later to filter low-quality bases, trim adaptors, and match cleaned reads to a reference genome, respectively. Make a directory for the project and separate folders for the raw SRA files, databases, fastq files, and the KneadData and MetaPhlAn outputs before downloading the corresponding reference databases for the specimen you're running an analysis on. Replace ```$DIR``` with the path to your database directory ```(e.g., ~/metagenomics_project/)```. Below is a list of the most common databases and their installation instructions:  
* Human Genome (default)
    - kneaddata_database --download human_genome bowtie2 $DIR/databases/kneaddata_db/
* Human Transcriptome
    - kneaddata_database --download human_transcriptome bowtie2 $DIR/databases/kneaddata_db/
* SILVA RNA
    - kneaddata_database --download ribosomal_RNA bowtie2 $DIR/databases/kneaddata_db/
* Specific Host Genome, such as mouse (mouse_C57BL) or cat (cat_genome)
    - kneaddata_database --download mouse_C57BL bowtie2 $DIR/databases/kneaddata_db/
    - kneaddata_database --download cat_genome bowtie2 $DIR/databases/kneaddata_db/

Download the raw SRA files and their corresponding metadata (this will be useful later when using R for further analysis) from [NCBI](https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/) and put them in the same folder designated for raw SRA files. SRA reads are either paired-end or single-end, with paired-end reads sequencing both ends of a DNA fragment, providing higher alignment confidence and enabling detection of structural rearrangements, while single-end reads sequence only one end, enabling a faster, simpler analysis process. Metagenomics is usually paired-end, but to double-check, read the metadata for the SRA file using 
```ruby
sra-stat --meta SRRXXXXXXX.sra
```

Merge the metadata file into one master file for later processing in R: 
```ruby
head -1 first_metadata.txt > metadata_merged.txt
tail -n +2 -q metadata_*.txt >> metadata_merged.txt
```

## Convert SRA to Fastq Files
Move to the folder designated to store fastq files using ```cd``` and run either of the following codes depending on whether the SRA file is paired-end or single-end. 
**Paired-end**
```ruby
fastq-dump --split-files --gzip $DIR/raw_sra/SRRXXXXXXXX/SRRXXXXXXXX.sra --outdir $DIR/fastq/
```
**Single-end**
```ruby
fastq-dump --gzip $DIR/raw_sra/SRR12345678/SRR12345678.sra --outdir $DIR/fastq/
```
**Batch Conversion for Paired-end**
```ruby
for sra_file in $DIR/raw_sra/*/*.sra; do
    fastq-dump --split-files --gzip "$sra_file" --outdir $DIR/fastq/
done
```
**Check whether the fastq files exist**
```ruby
ls $DIR/fastq
```

## Run KneadData
Move to the folder designated to store KneadData outputs using ```cd```. For paired-end reads, which is likely the case, run the following code (```--log``` is optional): 
```ruby
kneaddata \
    -i1 $DIR/fastq/SRR12345678_1.fastq.gz \
    -i2 $DIR/fastq/SRR12345678_2.fastq.gz \
    --reference-db $DIR/databases/kneaddata_db/ \
    --output $DIR/kneaddata_output/SRRXXXXXXXX/ \
    --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" \
    --threads 8 \
    --log $DIR/logs/SRRXXXXXXXX_kneaddata.log \
    --remove-intermediate-output
```

For single-end reads, run the following:
```ruby
kneaddata \
    -i $DIR/fastq/SRR12345678.fastq.gz \
    --reference-db $DIR/databases/kneaddata_db/ \
    --output $DIR/kneaddata_output/SRRXXXXXXXX/ \
    --threads 8 \
    --log $DIR/logs/SRRXXXXXXXX_kneaddata.log
```

### Some Key Concepts
**Trimmomatic Parameters:** 
- ILLUMINACLIP: cut specific Illumina adapter sequences to remove contaminants from the read.
- LEADING: cut low-quality bases below a certain threshold off the start of a read to reduce noise.
- TRAILING: cut low-quality bases below a certain threshold off the end of a read to reduce noise.
- SLIDINGWINDOW: perform a sliding window trimming that cuts when the average quality within the window falls below a specific threshold. 
- MINLEN: drop the read if it is below a specific length.
- MAXLEN: drop the read if it is longer than a specific length.
- AVGQUAL: drop the read if its average quality is below a specific threshold.

_If your KneadData crashes due to trimmomatic error, lower the MILEN gradually (e.g., 50 -> 36 -> 20) and/or the SLIDINGWINDOW (e.g., 4:20 -> 4:15; 4 is the number of RNA bases to average across and 20 is the average quality required, both numbers can be lowered to avoid error)._

**Threads** - determines the number of CPU cores used for parallel processing, affecting runtime and performance. Increasing the thread count decreases runtime and allows tools like BowTie2 and Trimmomatic to run simultaneously. The number of threads should not be arbitrary, but based on your device's CPU. To check your device's CPU, do the following and use a thread count less than the total CPU to ensure functionality of other services:
```
nproc #displays the total number of processing units available to WSL. 
lscpu #displays a detailed breakdown of CPU architecture, including the total number of CPUs, cores per socket, and threads per core. 
```

_If your KneadData crashes or runs on indefinitely, reduce the number of threads and try again._

### Run MetaPhlAn
Download the latest MetaPhlAn database. 
```ruby
metaphlan --install --bowtie2db $DIR/databases/metaphlan_db/
```
Then move to the folder designated to store MetaPhlAn outputs using ```cd``` and run the following code for paired-end read (```2>``` is optional as it has the function as log in KneadData): 
```ruby
metaphlan \   $DIR/kneaddata_output/SRRXXXXXXXX/SRRXXXXXXXX_paired_1.fastq,$DIR/kneaddata_output/SRRXXXXXXXX/SRRXXXXXXXX_paired_2.fastq \
    --bowtie2db $DIR/databases/metaphlan_db/ \
    --input_type fastq \
    --output_file $DIR/metaphlan_output/SRRXXXXXXXX_profile.txt \
    --bowtie2out $DIR/metaphlan_output/SRRXXXXXXXX.bowtie2.bz2 \
    --nproc 8 \
    2> $DIR/logs/SRRXXXXXXXX_metaphlan.log
```
For single-paired reads:
```ruby
metaphlan \
    $DIR/kneaddata_output/SRRXXXXXXXX/SRRXXXXXXXX_kneaddata_paired_1.fastq \
    --bowtie2db $DIR/databases/metaphlan_db/ \
    --input_type fastq \
    --output_file $DIR/metaphlan_output/SRRXXXXXXXX_profile.txt \
    --nproc 8
    2> $DIR/logs/SRRXXXXXXXX_metaphlan.log
```

Merge all MetaPhlAn output files into a single file in a table format to facilitate further analysis with the metadata using R. 
```ruby
merge_metaphlan_tables.py \
    $DIR/metaphlan_output/*_profile.txt \
    -o $DIR/metaphlan_output/merged_profiles.txt
```

### Some Key Concepts:
**nproc** - argument that defines the number of processor cores (threads) to use for parallel processing. --nproc 4 or higher is recommended for faster profiling of large metagenomic datasets. However, always check your device's CPU core and lower the number if MetaPhlAn crashes or runs on indefinitely. 

### Bonus Filters, Apply As Needed:
Extract species-level data:
```
grep -E "s__|clade" $DIR/metaphlan_output/merged_profiles.txt \
    > $DIR/metaphlan_output/species_profiles.txt
```
Extract genus-level data:
```
grep -E "g__|clade" $DIR/metaphlan_output/merged_profiles.txt \
    > $DIR/metaphlan_output/genus_profiles.txt
```
Extract phylum-level data:
```
grep -E "p__|clade" $DIR/metaphlan_output/merged_profiles.txt \
    > $DIR/metaphlan_output/phylum_profiles.txt
```

### Troubleshooting:
1. Check device disk space as the databases and files are large: ```df -h```
2. Monitor job process in another command prompt tab: ```tail -f $DIR/logs/SRRXXXXXXXX_kneaddata.log``` (This works for MetaPhlAn log files as well.)

## Use Google Drive and Google Colab if MetaPhlAn Keeps Crashing on the Device
Use the free 100GB 1-month plan with student credentials and upload the KneadData output files and MetaPhlAn BowTie2 databases to your Google Drive. Then, open Google Colab, activate Colab Pro using student credentials to get a 1-year free trial, and mount your drive to the colab page with the following code:
```ruby
from google.colab import drive
drive.mount('/PATH/TO/DRIVE_FOLDER')
```

Install MetaPhlAn and BowTie2 locally in Colab:
```ruby
!pip install metaphlan
!apt-get install -y bowtie2
```

Install BowTie2 Database to a designated folder:
```ruby
db_drive_path = "/PATH/TO/metatranscriptomics/metaphlan_db"
!mkdir -p {db_drive_path}
!metaphlan --install --bowtie2db /PATH/TO/metatranscriptomics/metaphlan_db
```

Run MetaPhlAn with ```nohup```, which keeps the process running even if the session disconnects or times out, and ```verbose```, which prints detailed progress messages while MetaPhlAn runs:
```ruby
!nohup metaphlan /PATH/TO/KNEADDATA_FOLDER/SRRXXXXXXXX_1_kneaddata_paired_1.fastq,/PATH/TO/KNEADDATA_FOLDER/SRRXXXXXXXX_1_kneaddata_paired_2.fastq \
    --input_type fastq \
    --nproc 4 \
    --bowtie2db /PATH/TO/metatranscriptomics/metaphlan_db/ \
    --bowtie2out /PATH/TO/METAPHLAN_OUTPUT/SRRXXXXXXXX_bowtie2.bz2 \
    -o /PATH/TO/METAPHLAN_OUTPUT/SRRXXXXXXXX_profile.txt \
    --verbose > /PATH/TO/METAPHLAN_OUTPUT/metaphlan_log.txt 2>&1 &

print("MetaPhlAn running in background. Check log file for progress.")
```

Check the progress in the log file and if MetaPhlAn is still running:
```ruby
!tail -f /PATH/TO/METAPHLAN_OUTPUT/metaphlan_log.txt
!ps aux | grep metaphlan
```

## Download and Set Up RStudio Web Server
While RStudio can be installed as a desktop software, the web version is lighter and easier to access with less lag. You can find instructions for installing RStudio on Linux or Windows [here](https://posit.co/download/rstudio-desktop). The following are the steps to install RStudio Server on Linux:
```ruby
sudo apt update #Update package list
sudo apt install -y r-base r-base-dev #Install R
R --version #Check version

wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2023.12.1-402-amd64.deb #Downloads R Server
sudo dpkg -i rstudio-server-2023.12.1-402-amd64.deb
sudo apt --fix-broken install -y #Install any missing dependencies
```

Once the downloads are complete, start the RStudio server:
```ruby
sudo systemctl start rstudio-server #Start the service
sudo systemctl enable rstudio-server #Enable it to start on boot when WSL activates
sudo systemctl status rstudio-server #Check server status
```

Access RStudio on the web from your browser and log in with your Linux username and password:```http://localhost:YOUR_SERVER_IP``` (usually, ```http://localhost:8787``` would work as the default for your machine). 

Manage and troubleshoot the RStudio server, especially when it lags due to an overwhelming input history:
```ruby
sudo systemctl stop rstudio-server #Stop service
rm -rf /home/yourusername/.local/rstudio/ #Clear RStudio session data for specific user
rm /home/yourusername/.RData #Clear RStudio workspace data
sudo systemctl restart rstudio-server #Restart service

sudo journalctl -u rstudio-server #Check logs if an error occurs
```

## R Code for Metatranscriptomic Analysis 
See the R code for various analyses in the other files of this repository. The following would offer simple explanations of the functionality of each analysis and how to interpret the output data. 

### Alpha Diversity Analysis 
Alpha diversity determines what and how many species are present in a particular environment within one sample. It measures the richness and evenness using metrics such as OTUs, Chao1 (estimated richness), and Shannon (diversity & evenness). The result is usually visualized with boxplots that compare groups utilizing non-parametric tests like Wilcoxon. 

**OTU (Operational Taxonomic Unit)** - groups closely related genetic sequences into clusters of taxonomic groups. Typically, sequences that are >= 97% identical are grouped into a single OTU, while sequences that are < 97% identical are split into different OTUs. The clustering process compares the sequences against a known database (e.g., Human, Mouse_C57BL). 
**Richness** - total number of taxonomic groups observed in a sample.
**Evenness** - the uniformity of individuals distributed among different species. 
**Chao1 Index** - Estimates the total richness by taking into account rare species.
**Shannon Index** - Accounts for both richness and evenness; it yields higher values when the sample has higher diversity.
**Simpson Index** - Measures the probability that two individuals randomly selected from a sample belong to the same species. 

Additional information on the topic could be found [here](https://cran.r-project.org/web/packages/tabula/vignettes/alpha.html)

### Beta Diversity Analysis
Beta diversity measures and quantifies how different two environments from different samples are. It requires a mathematical distance metric or dissimilarity matrix and statistical validation to produce an insightful output. 

**_Dissimilarity Matrix_** 
Shows how different each pair of samples is; choose the best metrics based on your data type.
**Bray-Curtis Index** - The most common metric for abundance data, and accounts for both the presence of species and their relative quantities. 

**Aitchison Distance** - often used in microbiomes as it handles the compositional bias by using centered log-ratio (clr). 

**UniFrac** - best choice for microbiome/metagenomic studies as unweighted UniFrac considers phylogenetic (evolutionary) relationships while weighted UniFrac considers the abundance of those lineages. 

**_Distance Metric_**
Distance matrices can be projected into 2D or 3D visual pattern. 
**PCoA (Principal Coordinates Analysis)** - The standard method for microbiome and ecological data that plots samples based on custom distance metrics. 

**NMDS (Non-metric Multidimensional Scaling)** - Ranks distances rather than absolute distances, which is highly useful when samples are extremely non-normal. 

**_Statistical Testing_**
Determines if the differences observed in the distance metrics are statistically significant using specific multivariate statistical tests. 
**PERMANOVA (Permutational Multivariate Analysis of Variance)** - The gold standard for determining if groupings (control v.s. treatment) explain a significant portion of the variation in community composition. 

**PERMDISP (Test of Multivariate Homogeneity of Group Dispersions)** - Used alongside PERMANOVA to ensure that the variation within the groups is roughly equal. If unequal, the PERMANOVA results may be skewed by variability rather than true differences in composition. 

### Rarefaction 
Evaluate sequencing depth and control for size biases when measuring active microbial diversity and gene expression using Phyloseq in R. It shows the number of unique transcripts or species identified as a function of the number of sequences sampled. 

**_Key Concepts & Purposes_**
1. Rarefaction curves help determine if the sequencing depth is sufficient to capture the active community's diversity and functional profile. If a curve plateaus, it indicates saturation and that additional sequencing will yield only marginal new information.
2. It adjusts for unequal sequencing depths across different samples by randomly subsampling a standardized number of reads from each sample.
3. It compares active taxonomic richness and functional gene diversity (a-diversity) across multiple samples.
4. Rarefaction curves saturate at different rates based on community activity and transcript complexity: highly active, low-diversity samples saturate quickly, while complex, varied communities would require further sequencing.
5. Rarefying metatranscriptomic data might discard valid data and reduce statistical power. Alternative methods for differential gene expression analysis include DESeq2 or TPM (Transcripts Per Million). 

### Cooccurrence 
Maps the concurrent expression of active genes across multiple microbial species within an environment. It identifies which microbes are actively responding to each other. 

**_Key Concepts & Purposes_**
1. Correlation statistics such as Pearson or Spearman are used to build networks where nodes represent microbial taxa or functional genes, and edges represent significant positive or negative co-expression/cooccurrence.
2. The analysis is most commonly used to identify microbial cross-feeding, where the metabolic byproducts of one species fuel the gene expression and growth of another.
3. It is also commonly used to identify specific microbial communities that are actively expressing genes for breaking down pollutants and analyzing how specific bacterial active networks in the human gut shift under the influence of diseases.  

### CLR-PCA
A combination of Centered Log-Ratio (CLR) transformation and Principal Component Analysis (PCA) to measure active gene expression within a microbial community. It resolves the arbitrary total read counts and raw abundances present in RNA sequenced data by stabilizing the variance across sparse matrices and transforming the data via log-ratio, allowing PCA to accurately reduce dimensionality and identify functional shifts across samples. 

**_Key Concepts & Purposes_**
1. CLR normalizes sequencing depths so that the data can be compared properly.
2. CLR-PCA is mathematically equivalent to generating an ordination based on Aitchison distances.
3. CLR transformation converts raw counts to relative abundances to account for the closed nature of the data. Applying PCA to the transformed data matrix extracts the primary Principal Components (PCs), explaining the most variance in gene expression.

### Core Microbiome Analysis
Identifies the microbial taxa consistently present across multiple samples within a specific environment. It identifies the stable, baseline microbial community. 

**_Three typical primary quantification methods:_**
1. **Occupancy/Prevalence** - The percentage of samples in which a specific microorganism is found.
2. **Relative Abundance** - The proportion of the total microbial community made up by a specific taxon across all samples.
3. **Hybrid** - a metric requiring an organism to exceed a certain minimum abundance and a specific detection threshold. (e.g., 95% prevalence at >= 0.01% abundance).

**_Purposes_**
1. Identify stable microbial interactions that act as an indicator of systemic health and predict or diagnose diseases.
2. Target imbalances during medical interventions by understanding the shared core microbiome across healthy populations. 

### Prevalence Abundance Analysis
Prevalence measures how often a specific microbe appears across a sample population, while abundance measures how much of that microbe is present. These metrics help identify core taxa, understand spatial distributions, and detect significant differences between experimental or disease groups. 

**_Key Concepts_**
* **Prevalence** - The percentage of samples in a dataset where a specific taxon is detected above a defined limit. This excludes low-quality reads and only analyzes the core microbiome. 
* **Abundance** - Evaluated as _relative abundance_ (the proportion a specific taxon represents the total microbial community) or _absolute abundance_ (the exact number or physical mass of cells). 
* **Prevalence Filtering** - Filters out taxa present in fewer than 10% to 50% of samples to reduce dimensionality of the dataset and computational noise before running statistical tests. 
* **Differential Abundance Analysis** - to find microbes that significantly differ between groups, algorithms that account for compositionality and zero-inflation, such as ANCOM-BC, LinDA (linear regression models on centered log-ratio transformed data aimed to correct biases with FDR control), and MaAsLin2 (multivariate analyses that adjust for metadata and confounding variables), are used. 

### Proportionality Analysis
A statistical method for identifying pairs of microbial taxa that maintain constant relative ratios across samples. It resolves the issue of compositionality where sequencing data yields only relative abundance, indicating a misleading negative correlation between unrelated taxa. 

Contrary to correlation metrics such as Pearson or Spearman, proportionality metrics are more suitable for relative abundance data and build highly accurate microbial association networks. 

**_Three Primary Proportionality Metrics_**
* **Phi** - A measure of proportional variances. If the ratio of two taxa remains constant, their log-ratio variances approach zero.
* **Rho** - A concordance-based metric ranging from -1 to 1. It serves a similar purpose to a correlation coefficient but is adjusted for compositional data.
* **Theta** - A proportionality metric that specifically evaluates the similarity between the relative abundances of two taxa.

**_Requirements & Considerations_** 
* A Microbiome dataset has many zeros, which is the reason it is imperative to filter or impute extreme sparsity or structural zeros before evaluation.

### Rank Abundance 
A graphical tool to visualize the composition of a microbial community. It plots the relative abundance of each microbial species on the y-axis against its abundance rank (from most to least) on the x-axis. 

In almost all healthy microbiomes, the rank abundance curve depicts a pattern of a long tail of rare microbes, with just one or a few highly abundant species. The most abundant constructs a large percentage of the total community, while the majority of taxa exist in very low numbers. 

**_How to Analyze_**
* **Species Richness** - Total number of distinct microbial species in a sample; indicated by the length of the curve.
* **Species Evenness** - The slope of the curve reveals how evenly the microbes are distributed. A steep curve indicates a community dominated by a few species, while a shallow slope indicates a balanced community.

### Taxa-Taxa Correlation
Measures how the abundances of different microorganisms vary together across samples. It is used to map microbial networks, identify cooccurrence or co-exlusion, and discover keystone species. It shares the same statistical challenges as proportionality analysis.

**_Interpretations:_**
* **Cooccurrence** - Suggests that taxa thrive in similar environmental conditions, cross-feed, or rely on the same metabolic byproducts.
* **Co-exclusion** - Suggests competitive exclusion, where taxa fight for the same limited nutrients or produce antimicrobial compounds that inhibit the growth of others.
* **Keystone Taxa** - Highly connected network hubs that exert massive influence over the entire microbiome's structure and function, despite the occasional low relative abundance. 

### Outlier Detection
Identifying samples or bacterial taxa that significantly deviate from the majority. It is a crucial preprocessing step for quality control and detecting abnormal microbial states or extreme variations in specific bacterial species. 

Conventional outlier methods often fall short when applied to microbiome data as it is highly dimensional, sparse, and zero-inflated. Thus, specialized techniques are used to detect samples heavily skewed by other factors. 

**_Some Common Methods:_**
* **CLOUD (Non-parametric Detection Test)_** - Evaluates conformity to determine if the microbiome data deviates from a healthy reference distribution.
* **PCoA/NMDS (Principal Coordinates Analysis/Non-metric Multidimensional Scaling Plots)** - Typically uses Bray-Curtis or UniFrac distance matrices to flag samples that fall far outside established experimental clusters.
* **Alpha Diversity Filters** - Evaluates within-sample diversity metrics (e.g., Shannon index or OUTs) and filters subjects with unusually low or high richness, which might indicate contamination or DNA extraction failures. 

### Sample Clustering Dendrogram 
A tree-like diagram that groups samples based on the similarity of their microbial communities. It visually identifies clusters where similar bacteria cluster closely together on the same branch. 

**_How to Analyze:_**
* **Leaves (Tips)** - Each individual sample.
* **Branches** - Lines connecting samples and create a node when joined with another branch.
* **Cluster** - Groups of branches that stem from a single, lower node. Samples within the same cluster share similar taxonomic profiles.
* **Branch Length/Height** - The vertical or horizontal scale representing the distance or dissimilarity between samples. The shorter the branch, the more similarities the bacterial composition in the microbiomes shares; the longer the branch, the more differences are present in the bacterial composition.

**_Common Analysis Tools:_**
1. **Distance Metric** - Measure how different the samples are with a distance matrix. Common metrics include: 
   * **Bray-Curtis** - Focuses on abundance data (how many of each bacterium are present).
   * **Jaccard** - Focuses on presence/absence (which bacteria are present, ignoring counts).
   * **UniFrac (Weighted/Unweighted)** - Incorporates phylogenetic relationships (how closely related the bacteria are to one another).
2. **Clustering Algorithm** - Groups the samples in the distance matrix hierarchically. Common algorithms include Ward's method, average linkage, or complete linkage.


Typically, dendrograms are plotted alongside _heatmaps_ or _metadata bars_. 
* **Heatmaps** - Displays the abundance of specific bacterial taxa (on the phylum or genus level), depicting exactly which bacteria are driving the clusters. 
* **Metadata Bars** - Color-coded bars indicating traits such as disease state, diet, or treatment, etc. 

### Cladogram
Displays taxonomic changes in the gut microbiome between groups. Common tools include [LEfSe](https://bioconductor.org/packages//release/bioc/vignettes/lefser/inst/doc/lefser.html) and [GraPhlAn](https://github.com/biobakery/graphlan/wiki). 

**LEfSe (Linear Discriminant Analysis Effect Size)** - discovers metagenomic biomarkers that identify, rank, and visualize features differentially abundant between two or more biological groups. It combines statistical significance with biological relevance to find features that explain differences between conditions. 
    * **Kruskal-Wallis Rank-Sum Test** - Detects features with significant differential abundance. 
    * **Wilcoxon Rank-Sum Test** - Evaluates biological consistency among subclasses. 
    * **Linear Discriminant Analysis (LDA)** - Estimates the effect size of each differentially abundant feature. LDA scores indicate the magnitude of differences. 

**GraPhlAn (Graphical Phylogenetic Analysis)** - generates high-quality, circular visualizations of taxonomic and phylogenetic trees that display evolutionary relationships alongside metadata, species abundance, biomarkers, and metabolic functions. It requires an input of a tree from pipelines like LEfSe. 

**_Key Features & Functions:_**
* Visualizes multi-level annotations (taxonomic ranks from Phylum to Strain).
* Displays abundance and phenotypic data using varying clade sizes, ring annotations, and gradients.
* Visually highlights bacterial markers that are different between experimental conditions.
* **Taxonomic Trees** - Illustrates the composition of microbial communities with nested rings.
* **Functional Ontologies** - Displays enriched metabolic pathways or gene functions.
* **Phylogenetic Structures** - Shows evolutionary distance and relationships. 

**_Rectangular v.s. Circular Cladogram:_**

**Rectangular Cladogram** - Best for small to medium-sized datasets. The horizontally extended branches are highly readable and facilitate tracing lineages, comparing specific taxa side-by-side, and reading long textual labels. 

**Circular Cladogram** - Best for large, complex datasets as it maximizes space efficiency that allows entire large-scale phylogenies to fit cleanly onto a single screen or page. 
