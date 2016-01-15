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
```

### To run a trio
```
source activate snakeenv
snakemake gemini/PRG_MAE_1.db
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
snakemake --config outdir="." -c "qsub -l h_vmem=40G -l mem_free=40G" -j 
```
