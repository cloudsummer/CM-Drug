## CM-Drug
![资源 12@600x](https://github.com/cloudsummer/CM-Drug/assets/24847317/3881a629-728e-41d1-9963-09eaf29d3ad9)


- [1. Download the docker image](#Download-the-docker-image)
- [2. Run the docker image](#Run-the-docker-image)
- [3. Download the file](#Download-the-file)
- [4. Running the code](#Running-the-code)
    

We utilized NEXTFLOW for the processing pipeline of bulk RNA-seq raw data of the datasets in ICBcomb, and the software used has been encapsulated within the Docker image: xynicoo/rstudio:4.3-5.

### Download the docker image
 ![资源 13@600x](https://github.com/cloudsummer/CM-Drug/assets/24847317/99e2bff0-685b-4264-a270-05856d663909)
The environments, dependencies, and toolkits required for our workflow have all been encapsulated into a Docker image. We utilized the Docker image to perform calculations.

Just simply run the following code on a server with Docker and NextFlow installed:

(Docker version we used is 20.10.21, build 20.10.21-0ubuntu1~22.04.3) 

```
#the docker image was made by me (the author of CM-Drug)
docker pull xynicoo/rstudio:4.3-5
```

### Run the docker image

All software parameters are preconfigured in the "parameters_of_RNAseq_workflow" file. If you need to modify the runtime parameters of the software, you can make changes to this file.

```
docker run -d --rm -p {YOUR-PORT}:8787 \
                -v /home/{YOUR-USER-NAME}:{THE PATH INSIDE THE DOCKER CONTAINTER YOU WANT TO MAP} \
                -e USER={YOUR-USER-NAME} -e PASSWORD={YOUR-PASSWORD} \
                -e USERID={YOUR-USERID (UID)} -e GROUPID={YOUR-GROUPID (GID)} -e ROOT=TRUE \
                --name {YOUR-DOCKER-CONTAINER} xynicoo/rstudio:4.3-5
```

To facilitate the use of Docker images, we have provided a shell script "docker_run_script" for starting Docker images. This shell script is provided as an example with the username XXX and R-studio server port 9999. Users can modify it according to their needs.

Modify the file "docker_run_script", and run the following code in the shell:
```
bash docker_run_script up 9999
```

if only the two lines of feedback shown in the following image are displayed, without any errors, it indicates that the Docker image has been successfully run.

<img width="400" alt="image" src="https://github.com/cloudsummer/CM-Drug/assets/24847317/f8333d97-2d0c-4456-a322-ac370d694bbf">



### Download the file

The size of the data for the input file is huge. Due to GitHub's file size limitations, we have uploaded the data to the Data Drive. You can download from following link:

https://drive.google.com/drive/folders/1sXLo5w_yuQDcS0XMMIsgpLCO9FliiGUf?usp=drive_link

### Running the code

After downloading the data, you will receive three folders representing three phases of LINCS data. Each time you run the code, copy one of the folders to your path and use 'setwd()' to set it as the working directory. This working directory corresponds to {work directory} in the code.

"CM-Drug-1.R" is used to process Phase 1 data, "CM-Drug-2.R" is used to process Phase 2 data, and "CM-Drug-New.R" is used to process the 2020 data.




