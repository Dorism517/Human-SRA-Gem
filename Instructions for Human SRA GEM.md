# Lab Overview
The experimental experince of RNA-seq processing can be extremely complex, and take up a lot of space depending on the data being analyzed. The use of the GEMmaker program is critical, as it allows the creation of a complex workflow in massive-scale RNA-seq analyses. The GEMmaker program will construct GEMs, Gene Expression count Matrices, which can show expression levels, unique reads, and various other useful experimental data. In order to do this, it is necessary to use SRA data. The SRA, Sequence Read Archive, is the **largest publicy available repository of high-throughput sequenceing data**. In this lab, we will aim to use human SRA experiments and the human genome to create GEMs, and then process the GEMs and analyze our findings. We will be using DNA-seq for our analysis, to identify differentially expressed genes in the Autism Purkinje neurons vs. the control Purkinje neurons, to attempt to identify genes that may contribute to the defects to Purkinje neurons that result in ASD behavior.

# Learning Objectives
- [ ] Understand the concept of GEM creation, including the necessary software to create a GEM
- [ ] Learn the processes that can be done to GEMs in order to analyze their data
- [ ] Analyze a downstream workflow of your created GEM
- [ ] Learn about austism spectrum disorder and the differences in Purkinje neurons between ASD patients and control patients

# Instructions for Human SRA GEM
**The overall goal of this project is to create a GEM from Human SRA files, and then normalize and analyze it**

First, research SRA experiments from the human genome within the NCBI sequence read archive. Some SRA experiments may have more than one run, all from the same study, which makes selections easier.

The Identified study below (PRJNA869106) had approximately 12 autism SRA uploads and 24 control uploads. **For simplicity, the table below lists the first 5 runs from each conditions which will be used to make the GEMs**, however further research can be done using all runs. The study involves looking at Purkinje neurons, which are located in the cerebral cortex of the brain. These neurons, when defected, have been discovered to cause system-wide autism spectrum disorder behavioral presentation. Using the transcriptomic analyzed human postmortem Purkinje neurons, we will assemble a control GEM and a autism GEM and use downstream analysis to compare the expression of the two.

| Project Identifier | Run Identifier | Description |
| :-----------: | :-----------: | :-----------: |
| PRJNA869106 | SRX17045810 | Control; Transcriptomic analysis of isolated human postmortem Purkinje neurons implicates developmental organization/connectivity, extracellular matrix organization, calcium ion response, immune function and signaling alterations in autism spectrum disorders |
| PRJNA869106 | SRX17045808 | Control; Transcriptomic analysis of isolated human postmortem Purkinje neurons implicates developmental organization/connectivity, extracellular matrix organization, calcium ion response, immune function and signaling alterations in autism spectrum disorders |
| PRJNA869106 | SRX17045807 | Control; Transcriptomic analysis of isolated human postmortem Purkinje neurons implicates developmental organization/connectivity, extracellular matrix organization, calcium ion response, immune function and signaling alterations in autism spectrum disorders |
| PRJNA869106 | SRX17045804 | Control; Transcriptomic analysis of isolated human postmortem Purkinje neurons implicates developmental organization/connectivity, extracellular matrix organization, calcium ion response, immune function and signaling alterations in autism spectrum disorders |
| PRJNA869106 | SRX17045802 | Control; Transcriptomic analysis of isolated human postmortem Purkinje neurons implicates developmental organization/connectivity, extracellular matrix organization, calcium ion response, immune function and signaling alterations in autism spectrum disorders |
| PRJNA869106 | SRX17045798 | Autism; Transcriptomic analysis of isolated human postmortem Purkinje neurons implicates developmental organization/connectivity, extracellular matrix organization, calcium ion response, immune function and signaling alterations in autism spectrum disorders |
| PRJNA869106 | SRX17045799 | Autism; Transcriptomic analysis of isolated human postmortem Purkinje neurons implicates developmental organization/connectivity, extracellular matrix organization, calcium ion response, immune function and signaling alterations in autism spectrum disorders |
| PRJNA869106 | SRX17045800 | Autism; Transcriptomic analysis of isolated human postmortem Purkinje neurons implicates developmental organization/connectivity, extracellular matrix organization, calcium ion response, immune function and signaling alterations in autism spectrum disorder |
| PRJNA869106 | SRX17045801 | Autism; Transcriptomic analysis of isolated human postmortem Purkinje neurons implicates developmental organization/connectivity, extracellular matrix organization, calcium ion response, immune function and signaling alterations in autism spectrum disorder |
| PRJNA869106 | SRX17045812 | Autism; Transcriptomic analysis of isolated human postmortem Purkinje neurons implicates developmental organization/connectivity, extracellular matrix organization, calcium ion response, immune function and signaling alterations in autism spectrum disorder |


## 
## Creating the Autism and Control GEMs

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

### Now, let's build the Autism GEM. To start, the complete human cDNA genome must be downloaded. 
1. Download the human genome: wget Homo_sapiens.GRCh38.cdna.all.fa.gz
    
    -gunzip Homo_sapiens.GRCh38.cdna.all.fa #unzip the file
        
    2. Use singularity to index the file, and name the indexed file with *.indexed* for differentiation purposes 
    
        -singularity exec -B ${PWD} https://depot.galaxyproject.org/singularity/kallisto:0.46.2--h4f7b962_1 kallisto index -i Homo_sapiens.GRCh38.cdna.all.fa.indexed Homo_sapiens.GRCh38.cdna.all.fa
        
2. Input in the selected SRA run experiment identifiers.

    - nano SRAsAutism.txt #create a text file

    - SRX17045798
    SRX17045799
    SRX17045800
    SRX17045801
    SRX17045812 #input the SRR identifiers of the chosen individual runs
    
3. Build the GEM. GEMmaker used NCBI to access the needed datasets and runs nextflow. Results will be seen in the **results directory**. If wanting to process the other runs, the process can be changed by altering the experiment ID's.

    -nextflow run systemsgenetics/gemmaker -profile singularity \
--pipeline kallisto \
--kallisto_index_path Homo_sapiens.GRCh38.cdna.all.fa.indexed \
--sras SRAsAutism.txt

4. Locate your created GEM

    -cd results
    
    -cd GEMs
    
    -ls
    
    -head GEMmaker.GEM.AUTISM.human.txt
    
### Now, let's build the Control GEM.
1. Input in the selected SRA run experiment identifiers.

    - nano SRAsControl.txt #create a text file

    - SRX17045810
    SRX17045808
    SRX17045807
    SRX17045804
    SRX17045802 #input the SRR identifiers of the chosen individual runs

2. Build the GEM. GEMmaker used NCBI to access the needed datasets and runs nextflow. Results will be seen in the **results directory**. If wanting to process the other runs, the process can be changed by altering the experiment ID's. We won't need to download the *Homo Sapiens* genome again as we already have it to use from the Autism GEM.

    -nextflow run systemsgenetics/gemmaker -profile singularity \
--pipeline kallisto \
--kallisto_index_path Homo_sapiens.GRCh38.cdna.all.fa.indexed \
--sras SRAsControl.txt

3. Locate your created GEM

    -cd results
    
    -cd GEMs
    
    -ls
    
    -head GEMmaker.GEM.CONTROL.human.txt

## 
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
    
4. Now, use the cloned GEM information to prepare your GEM for analysis. *If we were comparing multiple gems, then there would be an additional merge step* python ~/Desktop/classroom/myfiles/GEMmaker_runs/GEMprep/bin/merge.py {gem1} {gem2} *using the GEMprep software*. The below steps prep the GEMs for any analysis that could be chosen, but different analyses could require more preprocessing.

    -python ~/Desktop/classroom/myfiles/GEMmaker_runs/GEMprep/bin/merge.py GEMmaker.GEM.AUTISM.human.txt GEMmaker.GEM.CONTROL.human.txt merged-autism-control-gem.txt
    
    -python ~/Desktop/classroom/myfiles/GEMmaker_runs/GEMprep/bin/normalize.py merged-autism-control-gem.txt --log2 #Log2 transform your GEM
    
    -python ~/Desktop/classroom/myfiles/GEMmaker_runs/GEMprep/bin/normalize.py merged-autism-control-gem.txt --quantile #Quantile normalize
    
 
    
## 
## Analyzing a GEM
### Now, we will perform a downstream analysis on our GEM. The analytic technique we will be using is DESeq2, in order to find DEGS, differentially expressed genes.

### Processing our GEM for DESeq2 analysis
1. First, we will convert RSEM RNA-seq values into integers

    cat merged-autism-control-gem.txt | sed -e 's/\.[0-9]*//g' -e 's/ *$//' > merged-autism-control-gem-integer.txt

2. Next, we will remove duplicate gene rows and dashes, and convert to comma seperated again

    -cat merged-autism-control-gem-integer.txt | awk '!a[$1]++' > merged-autism-control-gem-integer-unique.txt #removes duplicates
    
    -cat merged-autism-control-gem-integer-unique.txt | sed 's/-/_/g' > merged-autism-control-gem-clean.txt
    
    -cat merged-autism-control-gem-clean.txt | sed 's/\s/,/g' >  merged-autism-control-gem-clean.csv

3. Next, convert to dashes, add necessary DESeq2 labels and header row, in order for the R code to recognize certain areas in the GEM

    -cat merged-autism-control.labels.txt| sed 's/-/_/g' >  merged-autism-control-dash.labels.txt #convert to dashes
    
    -cat merged-autism-control-dash.labels.txt | sed 's/\s/,/g' | sed 's/$/,HUMAN_AUTISM_CONTROL/' > merged-autism-control.comparison.tmp #Add labels
    
    -cat merged-autism-control.comparison.tmp | sed 's/sample,label,HUMAN_AUTISM_CONTROL/Sample,Group,Comparison/' > merged-autism-control.comparison.csv #Add header row




### Starting our analysis
1. First, we will install the DESeq2 using the following **R code in the R program**.

    -if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager") 
BiocManager::install("DESeq2")

    -a #updates all downloaded software, useful when operating within a directory which will save progres, prevents further repeat downloading and saves time

    -library(DESeq2)
    
2. Now, we will set our working directory. This will enable all generated results through R code to be sent to the correct directory with the previous work.

    -setwd("~/Desktop/classroom/myfiles/GEMmaker_runs")

3. Now, we will enter a number of functions to prep our DESeq2 run and define our terms

    #enter a function to extract a sub-matric of counts for each group
    -subgem <- function(gem, anot, group ){
  datalist = list()
  subanot = subset(anot, Comparison == group)
  for (id in subanot$Sample) {
    ind = which(colnames(gem) == id)
    genes = gem[0]
    exp = gem[,ind]
    datalist[[id]] <- exp
  }
  subcounts = cbind(genes, datalist)
  return(subcounts)
}

    #Enter a function to extract a subset of the sample annotation matrix for each group
    -subanot <- function(anot, group){
  datalist = list()
  print(str(group))
  subanot = subset(anot, Comparison == group)
  print(str(subanot))
  return(subanot)
}

    #Enter function to run DESeq2
    -run_deseq <- function(counts, annotation){
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = annotation,
                                design = ~ Group)
 #Filter and normalize genes with low total counts across all samples (Edit as needed): 
  dds <- dds[rowSums(counts(dds)) >= 50,]
  dds <- DESeq(dds)
  norm = fpm(dds)
#Sort the columns in the FPM data frame (Edit as needed): 
  conditionA = which(annotation[2] == "AUTISM_HUMAN")  #EDIT
  conditionB = which(annotation[2] == "CONTROL_HUMAN")  #EDIT
  norm = subset(norm, select=c(conditionA, conditionB))
  print(str(norm))
#Retrieve the results (Edit groupIDs as needed): 
  res <- results(dds, contrast=c("Group", "AUTISM_HUMAN", "CONTROL_HUMAN"))
  print(summary(res))
  res <- cbind(res, norm) # Add FPM values to results for easy visualization
  resultsNames(dds)
  return(res)

    -#Enter a function to print results to a file:
main <- function(countfile, anotfile, outfile){
  outname = outfile
  counts = read.delim(countfile, sep=',', header=TRUE, row.names='Hugo_Symbol') #EDIT
  samples = read.delim(anotfile, sep=',', row.names = NULL, check.names=FALSE)
  groups = unique(samples$Comparison)
  for (t in groups){
    subcounts = subgem(counts, samples, t)
    subannotation = subanot(samples, t)
    results = run_deseq(subcounts, subannotation)
    #Filter and sort results table
    f_results = subset(results, padj < 0.05)
    o_results = f_results[order(f_results$padj),]
    write.csv(o_results, outname, row.names = TRUE)
  }
  return(results)
}

4. Run analysis with the input file (The created merged Autism and control GEM, and the comparison GEM created)
    
    -main('merged-autism-control-gem.csv', 'merged-autism-control-gem.comparison.csv', 'human-autism-control-degs-tab-delim.csv')

    -Head human-autism-control-degs.csv | awk -F ',' '{print $1,$2,$3,$4,$5,$6,$7}' OFS=, > human-autism-control-degs-tab-delim.csv #print the specific needed columns and convert to tab-deliminated
    
5. The various listed genes in the first column of the outputted file can be compared using information such as the resulting log2FoldChange scores to determine the change in expression.


##
## Congratulations! You have successfully created a GEM from prexisting RNA-seq data, completed multiple processing steps, and analyzed the results of a comparison between a control and Autism condition for Purkinje cells!
