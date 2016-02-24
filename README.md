# grin
GRIN Epilepsy Study

### Snakemake setup
Do this once:
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
# type "yes"
# allow this to install into your home
conda config --add channels bioconda
conda config --add channels r
<<<<<<< HEAD
=======
conda config --add channels asmeurer
>>>>>>> master
conda create --name grinenv --file requirements.txt
source activate grinenv
```

In R. (We use the latest VariantFiltering)
```
source("http://bioconductor.org/biocLite.R")
library(BiocInstaller)
install.packages("devtools")
library(devtools)
<<<<<<< HEAD
install_github("rcastelo/VariantFiltering")
=======
install_github("VariantFiltering", "rcastelo")
>>>>>>> master
```

### To run a trio
```
source activate grinenv
snakemake GRCh37/gemini/PRG_MAE_1.gemini.db
```

### To generate a VCF of all trios, separately and combined
```
source activate grinenv
snakemake
```


## To run on Respublica
- `--drmaa` is not allowed on Respublica yet, use `-c qsub`
```
source activate grinenv
snakemake -c "qsub -l h_vmem=40G -l mem_free=40G" -j 
```
