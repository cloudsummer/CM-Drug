## CM-Drug
![资源 15@600x](https://github.com/cloudsummer/CM-Drug/assets/24847317/19d9cbc4-47e1-481c-8e31-c7ea3330b999)


- [1. Download the docker image](#Download-the-docker-image)
- [2. Run the docker image](#Run-the-docker-image)
- [3. Download the file](#Download-the-file)
- [4. Running the code](#Running-the-code)
    

### Download the docker image
![资源 16@600x](https://github.com/cloudsummer/CM-Drug/assets/24847317/c2e47ed2-07bc-4298-904b-5df33fb3b580)

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

To facilitate the use of Docker images, we have provided a shell script "docker_run_script" for starting Docker images. This shell script is provided as an example with the username XXX and R-studio server port 9999. Users can modify it according to their needs.

Modify the file "docker_run_script", and run the following code in the shell:
```
bash docker_run_script up 9999
```

if only the two lines of feedback shown in the following image are displayed, without any errors, it indicates that the Docker image has been successfully run.

<img width="400" alt="image" src="https://github.com/cloudsummer/CM-Drug/assets/24847317/f8333d97-2d0c-4456-a322-ac370d694bbf">



### Download the file

The size of the data for the input file is huge. Due to GitHub's file size limitations, we have uploaded the data to Google Drive. You can download it from the following link:

https://drive.google.com/drive/folders/1sXLo5w_yuQDcS0XMMIsgpLCO9FliiGUf?usp=drive_link

### Running the code

After downloading the data, you will receive three folders ("New", "1" and "2"), representing three phases of data.

<img width="309" alt="image" src="https://github.com/cloudsummer/CM-Drug/assets/24847317/1c967e97-1ef1-4445-9af1-ef7ee5850e2e">

Each time you run the code, copy one of the folders to your path and set it as the working directory. This working directory corresponds to {work directory} in the code.

"CM-Drug-1.R" is used to process Phase 1 data, "CM-Drug-2.R" is used to process Phase 2 data, and "CM-Drug-New.R" is used to process the 2020 data.


https://github.com/cloudsummer/CM-Drug/assets/24847317/5da42804-6bf3-425d-9bba-6387494dfc9d

[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/51E8Q5B4m3I/0.jpg)](https://youtu.be/51E8Q5B4m3I)



