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
conda config --add channels asmeurer
conda create --name grinenv --file requirements.txt
source activate grinenv
```

In R. (We use the latest VariantFiltering)
```
source("http://bioconductor.org/biocLite.R")
library(BiocInstaller)
install.packages("devtools")
library(devtools)
install_github("VariantFiltering", "rcastelo")
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
snakemake -j --cluster-config cluster.json -c "qsub -l h_vmem={cluster.mem} -l mem_free={cluster.mem} -pe smp {threads}"
```
