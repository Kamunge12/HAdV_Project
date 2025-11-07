Started by setting my conda enviroment for HPC
```
# downloaded the conda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# installed it to my local enviroment
bash Miniconda3-latest-Linux-x86_64.sh
# then initialized conda
source ~/miniconda3/bin/activate
conda init bash
# reloading conda
source ~/.bashrc
# created my bio enviroment
conda create -n bio python=3.10
# initialized the new enviroment
conda init
# activated my bio enviroment
conda activate bio
```

Installing the packages
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

installed toulligqc
```
conda install -c bioconda toulligqc
```

installed fastplong
```
conda install -c bioconda fastplong
```

installed Burrows-Wheeler transformation (bwa)
```
conda install -c bioconda bwa
```



