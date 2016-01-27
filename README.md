# grin
GRIN Epilepsy Study

### Snakemake setup
Do this once:
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
# type "yes"
# allow this to install into your home
conda create -n snakeenv python=3.4
source activate snakeenv
conda install -c bioconda snakemake
conda install pandas
conda install drmaa
conda config --add channels r
conda install -c https://conda.anaconda.org/bioconda rpy2
conda install -c https://conda.anaconda.org/bioconda bioconductor-variantannotation
conda install -c https://conda.anaconda.org/r r-dplyr
conda install -c https://conda.anaconda.org/bioconda bioconductor-biocgenerics
conda install -c https://conda.anaconda.org/bioconda bioconductor-variantannotation
```

### To run a trio
```
source activate snakeenv
snakemake GRCh37/gemini/PRG_MAE_1.gemini.db
```

### To generate a VCF of all trios, separately and combined
```
source activate snakeenv
snakemake
```


## To run on Respublica
- `--drmaa` is not allowed on Respublica yet, use `-c qsub`
```
source activate snakeenv
snakemake -c "qsub -l h_vmem=40G -l mem_free=40G" -j 
```
