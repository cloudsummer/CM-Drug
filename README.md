## CM-Drug
![资源 15@600x](https://github.com/cloudsummer/CM-Drug/assets/24847317/7c7daebb-974b-4262-919f-fcf56e82f440)


- [1. Download the docker image](#Download-the-docker-image)
- [2. Run the docker image](#Run-the-docker-image)
- [3. Download the file](#Download-the-file)
- [4. Running the code](#Running-the-code)

  "CM-Drug-New.R" is used to process the 2020 New data, "CM-Drug-1.R" is used to process Phase 1 data, "CM-Drug-2.R" is used to process Phase 2 data.
- [5. Demo](#Demo)
    

### Download the docker image
![资源 16@600x](https://github.com/cloudsummer/CM-Drug/assets/24847317/4fdaa516-e77b-4fdd-870e-f5b5707ab563)



The environments, dependencies, and toolkits required for our workflow have all been encapsulated into a Docker image. We utilized the Docker image to perform calculations.

Just simply run the following code on a server with Docker installed:

(Docker version we used is 20.10.21, build 20.10.21-0ubuntu1~22.04.3) 

```
#the docker image was made by the author of CM-Drug
docker pull xynicoo/rstudio:4.3-7
```

### Run the docker image

All software parameters are preconfigured in the "parameters_of_RNAseq_workflow" file. If you need to modify the runtime parameters of the software, you can make changes to this file.

```
docker run -d --rm -p {YOUR-PORT}:8787 \
                -v /home/{YOUR-USER-NAME}:{THE PATH INSIDE THE DOCKER CONTAINTER YOU WANT TO MAP} \
                -e USER={YOUR-USER-NAME} -e PASSWORD={YOUR-PASSWORD} \
                -e USERID={YOUR-USERID (UID)} -e GROUPID={YOUR-GROUPID (GID)} -e ROOT=TRUE \
                --name {YOUR-DOCKER-CONTAINER} xynicoo/rstudio:4.3-7
```

To facilitate the use of Docker images, we have provided a shell script "docker_run_script" for starting Docker images. This shell script is provided as an example with the username XXX and R-studio server port 9998. Users can modify it according to their needs.

Modify the file "docker_run_script", and run the following code in the shell:
```
bash docker_run_script up 9998
```

if only the two lines of feedback shown in the following image are displayed, without any errors, it indicates that the Docker image has been successfully run.

<img width="407" alt="image" src="https://github.com/cloudsummer/CM-Drug/assets/24847317/afe1e273-0e3f-4f8e-8322-df08c05f4e67">



### Download the file

The size of the data for the input file is huge. Due to GitHub's file size limitations, we have uploaded the data to Google Drive. You can download it from the following link:

https://drive.google.com/drive/folders/1sXLo5w_yuQDcS0XMMIsgpLCO9FliiGUf?usp=drive_link

### Running the code

After downloading the data, you will receive three folders ("New", "1" and "2"), representing three phases of data.

<img width="315" alt="image" src="https://github.com/cloudsummer/CM-Drug/assets/24847317/4d36ed29-8cfd-4632-842f-5902c4beb658">


Each time you run the code, copy one of the folders to your path and set it as the working directory. This working directory corresponds to {work directory} in the code.

"CM-Drug-New.R" is used to process the 2020 New data, "CM-Drug-1.R" is used to process Phase 1 data, "CM-Drug-2.R" is used to process Phase 2 data.

After building the Docker image and preparing the files, you can run the code. 

### Demo

Below is the video we created, demonstrating the code execution process using CM-Drug-New as an example.

(Click on the image below to be redirected to the corresponding YouTube video link.)

[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/51E8Q5B4m3I/0.jpg)](https://www.youtube.com/watch?v=51E8Q5B4m3I)



