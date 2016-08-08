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

freeze = config['freeze']

lastctg = config['lastctg'][freeze]


ENV3 = os.path.join(updir(shutil.which("conda"),3),config['python3_environment'],'bin') + '/'
ENV2 = os.path.join(updir(shutil.which("conda"),3),config['python2_environment'],'bin') + '/'

#hg37/hg38
freeze = config['freeze']

SLINK = "{{SLINK}}"


dlocs = dict()

workdir: config['projdir']

include:
    "rules/targets.rules"
include:
    "rules/utilities.rules"
    
### Sequence QC
include:
    "rules/sequence.rules"
RISAWSNAMES = ['E01621','E01623','E01622']
RISAWSFASTQ =[s for s in SAMPLELANES for r in RISAWSNAMES if r in s]
DOWNLOADDIR = "kiel"
DOWNLOADS = glob.glob(DOWNLOADDIR + "/*/fastq/*/*/*fastq.gz")
EXTRACT = [config['datadirs']['fastq'] + "/"+os.path.basename(name).rsplit(".")[0]+"_"+pair+".fastq.gz" for name in glob.glob(config['datadirs']['oldbams'] + "/*bam") for pair in ["R1","R2"]]
FASTQS = glob.glob(config['datadirs']['fastq'] + "/*.gz")
FASTQCS = [config['datadirs']['fastqc'] + "/" + re.sub("\.fastq.gz$", "_fastqc.zip", os.path.basename(name)) for name in FASTQS]
#FamilyID       Subject Mother  Father  Sex     Affected_status Not_in_Varbank
#Trio_SL        C2952   C2953   C2954   f       EOEE
#ISR_#45        E08320                          f       Focal Epilepsy  x
sample_table = pandas.read_table(config['sample_table'])
MANIFESTSAMPLES = list(set(list(sample_table['Subject'])+list(sample_table['Mother'].dropna())+list(sample_table['Father'].dropna())))
# ['E0974_GCTACGC_L006_R1', 'E0975_CGAGGCT_L006_R1', 'E0977_GTAGAGG_L006_R1', 'E0975_CGAGGCT_L006_R2', 'E0977_GTAGAGG_L006_R2', 'E0974_GCTACGC_L006_R2']
ALLPAIRNAMES = set([os.path.basename(name).split(os.extsep)[0] for name in FASTQS])
PAIRNAMESINSAMPLETABLE = [name for name in ALLPAIRNAMES for sample in MANIFESTSAMPLES if name.startswith(sample)] # some files may not be manifested
SAMPLESONDISK = [sample for name in ALLPAIRNAMES for sample in MANIFESTSAMPLES if name.startswith(sample)]
MISSINGSAMPLES = [sample for sample in MANIFESTSAMPLES if sample not in SAMPLESONDISK]
MANIFESTEDPAIRS = [name for name in ALLPAIRNAMES if name in PAIRNAMESINSAMPLETABLE]
UNMANIFESTEDPAIRS = [name for name in ALLPAIRNAMES if name not in PAIRNAMESINSAMPLETABLE]
SAMPLELANES = set([name.rsplit("_",maxsplit=1)[0] for name in PAIRNAMESINSAMPLETABLE]) # ['E0974_GCTACGC_L006', 'E0975_CGAGGCT_L006', 'E0977_GTAGAGG_L006']
EXISTINGSAMPLES = set([name.split("_",maxsplit=1)[0] for name in SAMPLELANES])

### Alignment
include:
    "rules/alignment.rules"
SAMS = [config['process_dir'][freeze] + config['results']['sams'] + "/" + name + ".sam" for name in SAMPLELANES]
BAMS = [config['process_dir'][freeze] + config['results']['bams'] + "/" + name + ".bam" for name in SAMPLELANES]
SBAMS = [config['process_dir'][freeze] + config['results']['bams'] + "/" + name + ".sorted.bam" for name in SAMPLELANES]
MBAMS = [config['process_dir'][freeze] + config['results']['bams'] + "/" + name + ".sorted.merged.bam" for name in EXISTINGSAMPLES]
DBAIS = [config['process_dir'][freeze] + config['results']['picard'] + "/" + name + ".rmdup.bai" for name in EXISTINGSAMPLES]
DBAMS = [config['process_dir'][freeze] + config['results']['picard'] + "/" + name + ".rmdup.bam" for name in EXISTINGSAMPLES]
GBAIS = [config['process_dir'][freeze] + config['results']['picard'] + "/" + name + ".group.bai" for name in EXISTINGSAMPLES]
GBAMS = [config['process_dir'][freeze] + config['results']['picard'] + "/" + name + ".group.bam" for name in EXISTINGSAMPLES]
RBAMS = [config['process_dir'][freeze] + config['results']['realigned'] + "/" + name + ".bam" for name in EXISTINGSAMPLES]
INDELS = config['process_dir'][freeze] + config['results']['realigned'] + "/indels.list"
LISTS = [config['process_dir'][freeze] + config['results']['lists'] + "/" + name + ".list" for name in EXISTINGSAMPLES]
TABLES = [config['process_dir'][freeze] + config['results']['recalibrated'] + "/" + name + ".table" for name in EXISTINGSAMPLES]
RECBAMS = [config['process_dir'][freeze] + config['results']['recalibrated'] + "/" + name + ".bam" for name in EXISTINGSAMPLES]
POSTTABLES = [config['process_dir'][freeze] + config['results']['postrecalibrated'] + "/" + name + ".table" for name in EXISTINGSAMPLES]
GATKPDFS = [config['process_dir'][freeze] + config['results']['pdfs'] + "/" + name + ".pdf" for name in EXISTINGSAMPLES]

### Variant Calling
include:
    "rules/variantcalling.rules"
GVCFS = [config['process_dir'][freeze] + config['results']['gvcfs'] + "/" + name + ".gvcf" for name in EXISTINGSAMPLES]
GVCFSLIST = ' '.join(["--variant " + config['process_dir'][freeze] + config['results']['gvcfs'] + "/" + name + ".gvcf" for name in EXISTINGSAMPLES])

### Annotation
include:
    "rules/annotation.rules"
ANNOVARDBS = [config['annovardbdir'] + "/" + freeze + "_" + db + ".installed" for db in config['annovardbs']]
ANNOVAR_PROTOCOLS = ','.join(config['annovardbs'])

### Analysis
include:
    "rules/analysis.rules"
COMPLETETRIOSFAMIDS = sorted(list(set([row['FamilyID']+'_'+row['Subject'] for index, row in sample_table.iterrows() if all([row[member] in EXISTINGSAMPLES for member in ['Mother','Father','Subject']])]))) # a quad produces two trios
TRIOVCFS = [config['landing_dir'][freeze] + config['results']['vcfs'] + "/" + trio + ".trio.phased.vcf" for trio in COMPLETETRIOSFAMIDS]
COMPLETEFAMILYFAMIDS = set([row['FamilyID'] for index, row in sample_table.iterrows() if all([row[member] in EXISTINGSAMPLES for member in ['Mother','Father','Subject']])]) # quads are one family
FAMILYVCFS = [config['landing_dir'][freeze] + config['results']['vcfs'] + "/" + trio + ".family.vcf" for trio in COMPLETEFAMILYFAMIDS]
COMVCFS = [config['landing_dir'][freeze] + config['results']['vcfs'] + "/" + trio + ".family.com.filtered.vcf" for trio in COMPLETEFAMILYFAMIDS]
VEPVCFS = [config['landing_dir'][freeze] + config['results']['vep'] + "/" + trio + ".family.com.filtered.vep.vcf" for trio in COMPLETEFAMILYFAMIDS]
INCOMPLETEFAMILIES = set([row['FamilyID'] for index, row in sample_table.iterrows() if any([row[member] not in EXISTINGSAMPLES and not pandas.isnull(row[member]) for member in ['Mother','Father','Subject']])])
TRIOGEMS = [config['landing_dir'][freeze] + config['results']['gemini'] + "/" + trio + ".gemini.db" for trio in COMPLETEFAMILYFAMIDS]
ANALYSISREADY = [config['landing_dir'][freeze] + config['results']['vcfs'] + "/" + trio + ".trio.phased.com.filtered.ad.de.nm.snpeff.noask.vcf.bgz" for trio in COMPLETETRIOSFAMIDS]
RDATA         = [config['landing_dir'][freeze] + config['results']['analysis'] + "/" + trio + ".trio.phased.com.filtered.ad.de.nm.snpeff.noask." + model + ".RData" for trio in COMPLETETRIOSFAMIDS for model in ['denovo','arhomo']]
ANALYSES      = [config['landing_dir'][freeze] + config['results']['analysis'] + "/" + trio + ".trio.phased.com.filtered.ad.de.nm.snpeff.noask.models.html" for trio in COMPLETETRIOSFAMIDS]


inlude:
    "rules/reporting.rules"

#### Internal
onsuccess:
    print("Workflow finished, no error")
    shell("mail -s 'workflow finished' "+config['admins']+" < {log}")

onerror:
    print("An error occurred")
    shell("mail -s 'an error occurred' "+config['admins']+" < {log}")

