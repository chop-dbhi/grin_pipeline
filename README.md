# grin
GRIN Epilepsy Study

### Snakemake setup
Do this once:
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
# type "yes"
conda create --name grinenv --file requirements.txt python=3.5
conda create --name geminienv --file python2_requirements.txt
source activate grinenv
```

### To run a trio
```
source activate grinenv
snakemake GRCh38/analysis/1KGCEU.trio.com.filtered.ad.de.nm.snpeff.models.html
```

### To generate a VCF of all trios, separately and combined
```
source activate grinenv
snakemake
```

### Setup on Respublica
Set a TMPDIR in your `~/.bash_profile`:
```
export TMPDIR=/mnt/lustre/users/YOURUSERNAME/scratch
```

Develop an appropriate config:
```
ln -s configs/config.respublica.yaml configs/config.yaml
cp configs/localconfig.sample.yaml localconfig.yaml
```

Choose a sample_table.txt
```
cd metadata
ln -s sample_table_full.txt sample_table.txt
```

Get some files:
```
mkdir fastq
cp /mnt/isilon/cbmi/variome/leipzig/GRIN/fastq/NA*gz fastq
```

- `--drmaa` is not allowed on Respublica yet, use `-c qsub`

### Run on Respublica
```
source activate grinenv
snakemake -j 300 --notemp --cluster-config configs/cluster.yaml -c "qsub -V -l h_vmem={cluster.h_vmem} -l mem_free={cluster.mem_free} -l m_mem_free={cluster.m_mem_free} -pe smp {threads}"
```

#### Screenshots
IGV screenshots are generated for each variant identified in the inheritance models. Unfortunately IGV cannot run in headless mode, and will require an X-window session to be started and IGV to be running on that machine. Each screenshot take 2-3 seconds to generate.

### Development
It's best to run in dev with --notemp since intermediates can take a long time to produce

To update the requirements file (after installing some new package):
```
conda list --explicit > requirements.txt
```

To update your conda environment with a new requirements file:
```
conda install --name grinenv --file  requirements.txt
```
You might find this will trigger Snakemake to want to remake downstream files, because executables are listed as input (this ensures they actually exist). To remedy this you can postdate any offending executables if you are confident their updates do not affect results.
```
touch -d 20160101 /home/leipzigj/miniconda3/envs/grinenv/bin/*
touch -h -d 20160101 /home/leipzigj/miniconda3/envs/grinenv/bin/*
touch -d 20150101 tools/GenomeAnalysisTK-3.6.tar.bz2 
```

### Setting up a new project

Create a table like this one:

[Sample table pros](https://github.research.chop.edu/BiG/grin/blob/master/metadata/sample_table_pros.txt)

Then create a sample.yaml file based off of this template:

[Sample yaml file template](https://github.research.chop.edu/BiG/grin/blob/master/configs/localconfig.sample.yaml)

(Deanne to write more on this)


