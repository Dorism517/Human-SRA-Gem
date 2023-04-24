# Instructions for Human SRA GEM
**The overall goal of this project is to create a GEM from Human SRA files, and then normalize and analyze it**

1. First, research SRA experiments from the human genome within the NCBI sequence read archive. Some SRA experiments may have more than one run, all from the same study, which makes selections easier.
2. The Identified study below (SRX000001) resulted in 10 different runs. For simplicity, the first 3 runs will be used to make the GEM, however further research can be done using all runs.

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
3. Next, in order to utilize the SRA experiment into a GEM, first the human genome must be downloaded and indexed. The genome can be indexed using the kallisto mapping software (make sure singularity container software is also installed so the software doesn't have to be installed locally)
4. Download the human genome: wget cdna.all.fa.gz
