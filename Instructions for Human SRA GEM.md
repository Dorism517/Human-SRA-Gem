# Lab Overview
The experimental experince of RNA-seq processing can be extremely complex, and take up a lot of space depending on the data being analyzed. The use of the GEMmaker program is critical, as it allows the creation of a complex workflow in massive-scale RNA-seq analyses. The GEMmaker program will construct GEMs, Gene Expression count Matrices, which can show expression levels, unique reads, and various other useful experimental data. In order to do this, it is necessary to use SRA data. The SRA, Sequence Read Archive, is the **largest publicy available repository of high-throughput sequenceing data**. In this lab, we will aim to use human SRA experiments and the human genome to create a GEM, and then process the GEM and analyze our findings.

# Learning Objectives
- [ ] Understand the concept of GEM creation, including the necessary software to create a GEM
- [ ] Learn the processes that can be done to GEMs in order to analyze their data
- [ ] Analyze a downstream workflow of your created GEM

# Instructions for Human SRA GEM
**The overall goal of this project is to create a GEM from Human SRA files, and then normalize and analyze it**

First, research SRA experiments from the human genome within the NCBI sequence read archive. Some SRA experiments may have more than one run, all from the same study, which makes selections easier.

The Identified study below (SRX000001) resulted in 10 different runs. **For simplicity, the first 3 runs will be used to make the GEM**, however further research can be done using all runs.

| Study | Run Identifier | Description |
| :-----------: | :-----------: | :-----------: |
| SRX000001 | SRR000021 | Paired-end mapping reveals extensive structural variation in the human genome |
| SRX000001 | SRR000026 | Paired-end mapping reveals extensive structural variation in the human genome |
| SRX000001 | SRR000027 | Paired-end mapping reveals extensive structural variation in the human genome |
| SRX000001 | SRR000034 | Paired-end mapping reveals extensive structural variation in the human genome |
| SRX000001 | SRR000045 | Paired-end mapping reveals extensive structural variation in the human genome |
| SRX000001 | SRR000054 | Paired-end mapping reveals extensive structural variation in the human genome |
| SRX000001 | SRR000059 | Paired-end mapping reveals extensive structural variation in the human genome |
| SRX000001 | SRR000063 | Paired-end mapping reveals extensive structural variation in the human genome |
| SRX000001 | SRR000065 | Paired-end mapping reveals extensive structural variation in the human genome |

### Creating a GEM

1. Next, in order to utilize the SRA experiment into a GEM, there are a few programs that must be installed:
    
    1. Install **Nextflow**; *a workflow manager that can manage pipelines/workflows of different complexities*

        -sudo apt update #Updates software repositories.
        
        -sudo apt install default-jre #Install Java
        
        -wget -qO- https://get.nextflow.io | bash #Compile nextflow program into one location
        
        -sudo cp nextflow /usr/local/bin #Adds nextflow to your PATH, change PATH command based on your path
        
        -nextflow #tests nextflow to confirm correct installation
        
    2. Install **Singularity**; *allows you to run containerized code, which contain software environments to run a code*
    
        -sudo apt update #Update the Linux software repositories
        
        -sudo apt-get install -y \build-essential \libseccomp-dev \pkg-config \squashfs-tools \cryptsetup #Install dependencies
        
        -export VERSION=1.16.6 OS=linux ARCH=amd64 #Install GO software, make sure to change based on operating system and versions released
        
        -wget -O /tmp/go${VERSION}.${OS}-${ARCH}.tar.gz \ https://dl.google.com/go/go${VERSION}.${OS}-${ARCH}.tar.gz
        
        -sudo tar -C /usr/local -xzf /tmp/go${VERSION}.${OS}-${ARCH}.tar.gz
        
        -echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
        
        -echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
        
        -source ~/.bashrc #load the changes to your path
        
        -go #check if GO software is successfully installed
    3. Install **singularity container software**; *virtual containers of software environments*
        
        -git clone https://github.com/sylabs/singularity.git
        
        -cd singularity
        
        -./mconfig
        
        -make -C builddir
        
        -sudo make -C builddir install
        
        -source ~/.bashrc #load changes to your path
        
        -singularity #test if singularity is successfully installed

2. Now that the necessary software is installed, GEMmaker can utilize these programs to run. It's important to create a working directory for GEMmaker, and test it. 
    1. Create a working directory, for example called GEMmaker_runs. For this labs example purposes, this directory will be located in ~/Desktop/classroom/myfiles/GEMmaker_runs
    
        -cd ~/Desktop/classroom/myfiles
        
        -mkdir GEMmaker_runs
        
        -cd GEMmaker_runs
        
    2. Now, within your GEM working directory, run a test on GEMmaker. *The goal of this is to see a count of TMP gene expression matrix in the results section. As it is a fake genome and run just to test the program, the resulting GEMs will be very small, which can be seen using the command* ls -l *to list the line count/size.*
        
        -nextflow run systemsgenetics/gemmaker -profile test,singularity
        
        -cd results #check results for the expected files
        
        -cd ..
        
        -cd reports #look for a QC analysis of the SRA files

3. Now, let's build the *Homo Sapiens* GEM. To start, the complete human cDNA genome must be downloaded. 
    1. Download the human genome: wget Homo_sapiens.GRCh38.cdna.all.fa.gz
    
        -gunzip Homo_sapiens.GRCh38.cdna.all.fa #unzip the file
        
    2. Use singularity to index the file, and name the indexed file with *.indexed* for differentiation purposes 
    
        -singularity exec -B ${PWD} https://depot.galaxyproject.org/singularity/kallisto:0.46.2--h4f7b962_1 kallisto index -i Homo_sapiens.GRCh38.cdna.all.fa.indexed Homo_sapiens.GRCh38.cdna.all.fa
        
4. Input in the selected SRA run experiment identifiers.

    - nano SRAs.txt #create a text file

    - SRR000021
    SRR000021
    SRR000021 #input the SRR identifiers of the chosen individual runs
    
5. Build the GEM. GEMmaker used NCBI to access the needed datasets and runs nextflow. Results will be seen in the **results directory**. If wanting to process the other runs, the process can be changed by altering the experiment ID's.

    -nextflow run systemsgenetics/gemmaker -profile singularity \
--pipeline kallisto \
--kallisto_index_path Homo_sapiens.GRCh38.cdna.all.fa.indexed \
--sras SRAs.txt

6. Locate your created GEM

    -cd results
    
    -cd GEMs
    
    -ls
    
    -head GEMmaker.GEM.TPM.human.txt

## Processing a GEM
### Now that the GEM has been created using the selected SRA experiment runs, it is time to prep our GEM in order to analyze it.

1.Create and activate a GEM prep environment

    -conda create -n gemprep python=3.6 matplotlib mpi4py numpy pandas r scikit-learn seaborn
    
    -source activate gemprep #activate the created GEMprep environment, to deactivate hit CTRL+D
    
2. Verify that the environment has changed from the *(base)* environment to the *(gemprep)* environment. This should be visible at the beginning of the coding line.
3. Now that the GEMprep environment has been activated, clone the GEMprep program into your working directory

    -cd ~/Desktop/classroom/myfiles/GEMmaker_runs

    -git clone https://github.com/SystemsGenetics/GEMprep

    -cd GEMprep
    
4. Now, use the cloned GEM information to prepare your GEM for analysis. *If we were comparing multiple gems, then there would be an additional merge step* python ~/Desktop/classroom/myfiles/GEMmaker_runs/GEMprep/bin/merge.py {gem1} {gem2} *using the GEMprep software*

    -python ~/Desktop/classroom/myfiles/GEMmaker_runs/GEMprep/bin/normalize.py GEMmaker.GEM.TPM.human.txt --log2 #Log2 transform your GEM
    
    -python ~/Desktop/classroom/myfiles/GEMmaker_runs/GEMprep/bin/normalize.py GEMmaker.GEM.TPM.human.txt --quantile #Quantile normalize
    
## Analyzing a GEM
