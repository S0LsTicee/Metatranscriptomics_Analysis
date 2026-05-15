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

### Some Key Concepts:
**nproc** - argument that defines the number of processor cores (threads) to use for parallel processing. --nproc 4 or higher is recommended for faster profiling of large metagenomic datasets. However, always check your device's CPU core and lower the number if MetaPhlAn crashes or runs on indefinitely. 

Merge all MetaPhlAn output files into a single file in a table format to facilitate further analysis with the metadata using R. 
```ruby
merge_metaphlan_tables.py \
    $DIR/metaphlan_output/*_profile.txt \
    -o $DIR/metaphlan_output/merged_profiles.txt
```
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

Manage the RStudio server, such as when the RStudio server lags due to an overwhelming history input:
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


### Beta Diversity Analysis
