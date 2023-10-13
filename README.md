## CM-Drug
![资源 12@600x](https://github.com/cloudsummer/CM-Drug/assets/24847317/3881a629-728e-41d1-9963-09eaf29d3ad9)


- [1. Download the docker image](#Download-the-docker-image)
- [2. Run the docker image](#Run-the-docker-image)
- [3. Download the file](#Dowload-the-file)
- [4. Running the code](#Running-the-code)
    

We utilized NEXTFLOW for the processing pipeline of bulk RNA-seq raw data of the datasets in ICBcomb, and the software used has been encapsulated within the Docker image: xynicoo/rnaseq:n3-fastpMqc.

### Download the docker image
 ![资源 13@600x](https://github.com/cloudsummer/CM-Drug/assets/24847317/99e2bff0-685b-4264-a270-05856d663909)
The environments, dependencies, and toolkits required for our workflow have all been encapsulated into a Docker image. We utilized the Docker image to perform calculations.


In our workflow
Just simply run the following code on a server with Docker and NextFlow installed:

(Docker version we used is 20.10.21, build 20.10.21-0ubuntu1~22.04.3; the NextFlow version we used is 20.07.1 build 5412) 

```
#the docker image was made by me (the author of ICBcomb)
docker pull xynicoo/rnaseq:n3-fastpMqc
```

### Run the docker image

Modify the file of "parameters_of_RNAseq_workflow"

All software parameters are preconfigured in the "parameters_of_RNAseq_workflow" file. If you need to modify the runtime parameters of the software, you can make changes to this file.

Software detail in the docker image "xynicoo/rnaseq:n3-fastpMqc":

- FastQC (version 0.11.9) was used for data quality control (QC).

- Fastp (version 0.23.1) was employed for adapter sequence removal and trimming to obtain high-quality clean reads. 
 
- Clean reads were mapped to the human reference genome GRCh38 or the mouse reference genome GRCm39 by HISAT2 (version 2.2.1).
 
- SAMtools (version 1.16) was used to convert the “.sam” file into a “.bam” file.
 
- StringTie (version 2.2.1) was used to estimate the abundance of transcripts for each sample.
 
- FeatureCounts (version 2.0.3) was used to calculate gene expression and get the raw counts (reads matrix).

### Download the file

Modify the file "nextflow.config"

To run a Nextflow configuration file and specify parameters such as the path to the fastq files, reference genome path, user UID, and other relevant settings.

### Running the code

Running the following code will initiate background processing, and save the log in "NF.log":

```$ nohup nextflow ./parameters_of_RNAseq_workflow -with-docker xynicoo/rnaseq:n3-fastpMqc -c nextflow.config >> NF.log 2>&1 &```


