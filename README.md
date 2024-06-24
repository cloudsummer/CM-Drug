## CM-Drug

![313563468-837b6b4e-b59d-4568-8224-e52965d5634b](https://github.com/cloudsummer/CM-Drug.1/assets/24847317/f211bceb-7480-4c60-80c0-190fb0359b9f)


- [1. Example of the docker image structure](#Example-of-the-docker-image-structure)
- [2. Example to run the docker image](#Example-to-run-the-docker-image)
- [3. Download the file](#Download-the-file)
- [4. Running the code](#Running-the-code)
- [5. Demo](#Demo)
    

### Example of the docker image structure

![313563494-c2999480-31af-4736-8e42-2809e69a354b](https://github.com/cloudsummer/CM-Drug.1/assets/24847317/72079658-0148-4066-a954-a5a73034029e)


You can encapsulate the environments, dependencies, and toolkits required for the workflow into a Docker image. Each time you change the work environment, you just need to utilize the Docker image to perform calculations. Here, we provide a Docker image with commonly used packages integrated with this project. The image name is xynicoo/rstudio:4.3-7.


The Docker version we used is 20.10.21, build 20.10.21-0ubuntu1~22.04.3. You can use a different version of Docker or switch to a similar tool like Singularity images, which can be used on a cluster and does not require root privileges, depending on your needs and the server conditions you are using.

### Example to run the docker image


```
docker run -d --rm -p {YOUR-PORT}:8787 \
                -v /home/{YOUR-USER-NAME}:{THE PATH INSIDE THE DOCKER CONTAINTER YOU WANT TO MAP} \
                -e USER={YOUR-USER-NAME} -e PASSWORD={YOUR-PASSWORD} \
                -e USERID={YOUR-USERID (UID)} -e GROUPID={YOUR-GROUPID (GID)} -e ROOT=TRUE \
                --name {YOUR-DOCKER-CONTAINER} {Docker-Image-Name}
```

To facilitate the use of Docker images, we have provided a shell script "docker_run_script" for starting Docker images. This shell script is provided as an example with the username XXX and R-studio server port 9999. Users can modify it according to their needs.

Modify the file "docker_run_script", and run the following code in the shell:
```
bash docker_run_script up 9999
```

if only the two lines of feedback shown in the following image are displayed, without any errors, it indicates that the Docker image has been successfully run.

<img width="464" alt="313563540-8cdb32b0-05eb-4ac8-b187-fbb067271d9d" src="https://github.com/cloudsummer/CM-Drug.1/assets/24847317/23d18678-876a-49f7-b4f2-8b766df8146a">




### Download the file

The size of the data for the input file is huge. Due to GitHub's file size limitations, we have uploaded the data to Google Drive. You can download it from the following link:

[https://drive.google.com/drive/folders/1JdFmlEAJA66ZzyUP2ZjrIiTKFBGOrRLi?usp=drive_link](https://drive.google.com/drive/folders/1JdFmlEAJA66ZzyUP2ZjrIiTKFBGOrRLi?usp=drive_link)


### Running the code

After downloading the data, you will receive the folder named "data", which contains the data.

<img width="118" alt="image" src="https://github.com/cloudsummer/CM/assets/24847317/7c097fa7-851a-4742-9e34-37f7ace72fa7">

If you run the code, copy the folders ("data") to your path and set it as the working directory. This working directory corresponds to {work directory} in the code.

"CM-Drug.R" is the code.

(The data from LINCS2020 is the most recent, published in 2021. The new dataset not only contains a significantly larger amount of data but has also provided more abundant meta-information on compounds. It is recommended to utilize the new data.)

If you use Docker, after creating the Docker container, preparing the files, and setting the working directory, you can run the code. 

### Demo

Below is the video we created, demonstrating the code execution process.

(Click on the image below to be redirected to the corresponding YouTube video link.)

[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/51E8Q5B4m3I/0.jpg)](https://www.youtube.com/watch?v=51E8Q5B4m3I)



