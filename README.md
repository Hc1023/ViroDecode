# ViroDecode

Helper Scripts running [ViroDecode](https://hub.docker.com/repository/docker/sisih/ViroDecode).

## Using Docker (Require root privilege)

Run ViroDecode via the Docker image. For your convenience, we have created a script, [run_virodecode_docker.py](https://github.com/Hc1023/ViroDecode/blob/main/run_virodecode_docker.py), which handles Docker-related issues such as mounting files for you and is cross-platform compatible. 

1. Install Docker and Python on your machine.
2. Pull the image by Docker.
```sh
sudo docker pull sisih/viraldecode:1
```
3. Run the python script.
```sh
python run_virodecode_docker.py -i sisih/viraldecode:1 -1 fq1 -2 fq2
```

## Using Singularity (Recommended for HPC)

For security reasons, many high performance computing resources disallow Docker usage. Instead, many of these platforms support Singularity, which is very similar to Docker. Thus, if you are using a compute resource that supports Singularity, you can use the [run_virodecode_singularity.py](https://github.com/Hc1023/ViroDecode/blob/main/run_virodecode_singularity.py) script.

1. Load Singularity and Python on your computer cluster. This must be done each time you log in (or you can add it to your `.bashrc` for convenience). 
2. Pull the image by Singularity.
```sh
singularity pull docker://sisih/viraldecode:1
```
3. Run the python script.
```sh
python run_virodecode_singularity.py -i virodecode_1.sif -1 fq1 -2 fq2
```