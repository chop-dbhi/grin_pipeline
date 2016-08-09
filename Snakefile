import sys
import glob
import re
import pandas
import yaml
import subprocess
import configparser
import shutil
from snakemake.utils import R
from functools import cmp_to_key
from util import varsub
from util import lists
"""
run on respublica

source activate grinenv

snakemake -j 300 --cluster-config configs/cluster.yaml -c "qsub -V -l h_vmem={cluster.g_vmem} -l mem_free={cluster.mem_free} -l mem_free_l={cluster.m_mem_free} -pe smp {threads}"

steps run in the order:

rule symlinks    # make symlinks in fastq to original fastq files
rule validateBam # validate oldbams (may be skipped)
rule extract     # extrat fastq from oldbams files
rule make_yaml   # make QC files (invoke rule fastqc, may be skipped)
rule make_sams   # make sam files (invoke rule align)
rule make_bams   # make bam fiels (invoke rule sam_to_bam)
rule merge_lanes # get sorted bam files (invoke rule novosortbam)
rule make_merged # get sorted merged bam file (invoke rule merge_lanes)
rule family_vcfs # make family.vcfs
rule triovcfs    # make trio.phased (invoke rule run_phase_by_transmission)
rule vepvcfs     # run vep
rule xbrowse     # get files for xbrowse
"""

shell.prefix("source ~/.bash_profile;") 

configfile: "configs/baseconfig.yaml"
configfile: "localconfig.yaml"     # copied and customized from configs/localconfig.sample.yaml

varsub(config)  # substitute $isilon variable

ENV3 = os.path.join(updir(shutil.which("conda"),3),config['python3_environment'],'bin') + '/'
ENV2 = os.path.join(updir(shutil.which("conda"),3),config['python2_environment'],'bin') + '/'

#hg37/hg38
freeze = config['freeze']

SLINK = "{{SLINK}}"

workdir: config['projdir']

include:
    "rules/targets.rules"
include:
    "rules/utilities.rules"
include:
    "rules/sequence.rules"
include:
    "rules/alignment.rules"
include:
    "rules/variantcalling.rules"
include:
    "rules/annotation.rules"
include:
    "rules/analysis.rules"
include:
    "rules/reporting.rules"



#### Internal
onsuccess:
    print("Workflow finished, no error")
    shell("mail -s 'workflow finished' "+config['admins']+" < {log}")

onerror:
    print("An error occurred")
    shell("mail -s 'an error occurred' "+config['admins']+" < {log}")

