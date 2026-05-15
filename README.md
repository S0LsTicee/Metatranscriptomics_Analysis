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
    kneaddata_database --download human_genome bowtie2 $DIR/databases/kneaddata_db/
* Human Transcriptome
    kneaddata_database --download human_transcriptome bowtie2 $DIR/databases/kneaddata_db/
* SILVA RNA
    kneaddata_database --download ribosomal_RNA bowtie2 $DIR/databases/kneaddata_db/
* Specific Host Genome, such as mouse (mouse_C57BL) or cat (cat_genome)
    kneaddata_database --download mouse_C57BL bowtie2 $DIR/databases/kneaddata_db/
    kneaddata_database --download cat_genome bowtie2 $DIR/databases/kneaddata_db/

Download the raw SRA files from NCBI and put them in the same folder designated for raw SRA files. SRA reads are either paired-end or single-end, with paired-end reads sequencing both ends of a DNA fragment, providing higher alignment confidence and enabling detection of structural rearrangements, while single-end reads sequence only one end, enabling a faster, simpler analysis process. Metagenomics is usually paired-end, but to double-check, read the metadata for the SRA file using 
```
sra-stat --meta SRRXXXXXXX.sra
```
## Convert SRA to Fastq Files
Move to the folder designated to store fastq files using ```cd``` and run either of the following codes depending on whether the SRA file is paired-end or single-end. 
**Paired-end**
```
fastq-dump --split-files --gzip $DIR/raw_sra/SRRXXXXXXXX/SRRXXXXXXXX.sra --outdir $DIR/fastq/
```
**Single-end**
```
fastq-dump --gzip $DIR/raw_sra/SRR12345678/SRR12345678.sra --outdir $DIR/fastq/
```
**Batch Conversion for Paired-end**
```
for sra_file in $DIR/raw_sra/*/*.sra; do
    fastq-dump --split-files --gzip "$sra_file" --outdir $DIR/fastq/
done
```
**Check whether the fastq files exist**
```
ls $DIR/fastq
```

## Run KneadData
Move to the folder designated to store KneadData outputs using ```cd```. For paired-end reads, which is likely the case, run the following code (```--log``` is optional): 
```
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
```
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

**Threads** - determines the number of CPU cores used for parallel processing, affecting runtime and performance. Increasing the thread count decreases runtime and allows tools like BowTie2 and Trimmomatic to run simultaneously. The number of threads should not be arbitrary, but based on your device's CPU. To check your device's CPU, do the following and use a thread count less than the total CPU to ensure functionality of other services:
_Window's_
```
```

_Mac iOS_ 
```
```

### Run MetaPhlAn
Move to the folder designated to store MetaPhlAn outputs using ```cd```. 
