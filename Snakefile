import glob
import re
import pandas
import yaml
import subprocess
import configparser
from snakemake.utils import R
from functools import cmp_to_key
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
rule xbrowse    # get files for xbrowse
"""

shell.prefix("source ~/.bash_profile;") 

configfile: "configs/baseconfig.yaml"
configfile: "configs/config.yaml"
configfile: "test.yaml"

ENV3 = '{condaenv}/'.format(condaenv=config['python3_environment'])
ENV2 = '{condaenv}/'.format(condaenv=config['python2_environment'])

genome = config['buildve']

SLINK = "{{SLINK}}"

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

# some files may not be manifested
PAIRNAMESINSAMPLETABLE = [name for name in ALLPAIRNAMES for sample in MANIFESTSAMPLES if name.startswith(sample)]
SAMPLESONDISK = [sample for name in ALLPAIRNAMES for sample in MANIFESTSAMPLES if name.startswith(sample)]
MISSINGSAMPLES = [sample for sample in MANIFESTSAMPLES if sample not in SAMPLESONDISK]
MANIFESTEDPAIRS = [name for name in ALLPAIRNAMES if name in PAIRNAMESINSAMPLETABLE]
UNMANIFESTEDPAIRS = [name for name in ALLPAIRNAMES if name not in PAIRNAMESINSAMPLETABLE]

# pair up
# ['E0974_GCTACGC_L006', 'E0975_CGAGGCT_L006', 'E0977_GTAGAGG_L006']
SAMPLELANES = set([name.rsplit("_",maxsplit=1)[0] for name in PAIRNAMESINSAMPLETABLE])

EXISTINGSAMPLES = set([name.split("_",maxsplit=1)[0] for name in SAMPLELANES])

# a quad produces two trios
COMPLETETRIOSFAMIDS = sorted(list(set([row['FamilyID']+'_'+row['Subject'] for index, row in sample_table.iterrows() if all([row[member] in EXISTINGSAMPLES for member in ['Mother','Father','Subject']])])))
TRIOVCFS = [config['types'][genome] + config['results']['vcfs'] + "/" + trio + ".trio.phased.vcf" for trio in COMPLETETRIOSFAMIDS]

# quads are one family
COMPLETEFAMILYFAMIDS = set([row['FamilyID'] for index, row in sample_table.iterrows() if all([row[member] in EXISTINGSAMPLES for member in ['Mother','Father','Subject']])])
FAMILYVCFS = [config['types'][genome] + config['results']['vcfs'] + "/" + trio + ".family.vcf" for trio in COMPLETEFAMILYFAMIDS]
VEPVCFS = [config['types'][genome] + config['results']['vep'] + "/" + trio + ".family.com.filtered.vep.vcf" for trio in COMPLETEFAMILYFAMIDS]

# VEPVCFS = glob.glob(config['types'][genome] + config['results']['vcfs'] + "/*.vcf")
# VEPVCFS = [re.sub("/vcfs/", "/vep/", name) for name in VEPVCFS]
# VEPVCFS = [re.sub("vcf$", "com.filtered.vep.vcf", name) for name in VEPVCFS]

#for name in VEPVCFS:
#   print(name)
#quit()

INCOMPLETEFAMILIES = set([row['FamilyID'] for index, row in sample_table.iterrows() if any([row[member] not in EXISTINGSAMPLES and not pandas.isnull(row[member]) for member in ['Mother','Father','Subject']])])
TRIOGEMS = [config['types'][genome] + config['results']['gemini'] + "/" + trio + ".gemini.db" for trio in COMPLETEFAMILYFAMIDS]

ANALYSISREADY = [config['types'][genome] + config['results']['vcfs'] + "/" + trio + ".trio.phased.com.filtered.ad.de.nm.snpeff.noask.vcf.bgz" for trio in COMPLETETRIOSFAMIDS]
RDATA         = [config['types'][genome] + config['results']['analysis'] + "/" + trio + ".trio.phased.com.filtered.ad.de.nm.snpeff.noask." + model + ".RData" for trio in COMPLETETRIOSFAMIDS for model in ['denovo','arhomo']]
ANALYSES      = [config['types'][genome] + config['results']['analysis'] + "/" + trio + ".trio.phased.com.filtered.ad.de.nm.snpeff.noask.models.html" for trio in COMPLETETRIOSFAMIDS]

SAMS = [config['types'][genome] + config['results']['sams'] + "/" + name + ".sam" for name in SAMPLELANES]
BAMS = [config['types'][genome] + config['results']['bams'] + "/" + name + ".bam" for name in SAMPLELANES]
SBAMS = [config['types'][genome] + config['results']['bams'] + "/" + name + ".sorted.bam" for name in SAMPLELANES]
MBAMS = [config['types'][genome] + config['results']['bams'] + "/" + name + ".sorted.merged.bam" for name in EXISTINGSAMPLES]
DBAIS = [config['types'][genome] + config['results']['picard'] + "/" + name + ".rmdup.bai" for name in EXISTINGSAMPLES]
DBAMS = [config['types'][genome] + config['results']['picard'] + "/" + name + ".rmdup.bam" for name in EXISTINGSAMPLES]
GBAIS = [config['types'][genome] + config['results']['picard'] + "/" + name + ".group.bai" for name in EXISTINGSAMPLES]
GBAMS = [config['types'][genome] + config['results']['picard'] + "/" + name + ".group.bam" for name in EXISTINGSAMPLES]
RBAMS = [config['types'][genome] + config['results']['realigned'] + "/" + name + ".bam" for name in EXISTINGSAMPLES]
LISTS = [config['types'][genome] + config['results']['lists'] + "/" + name + ".list" for name in EXISTINGSAMPLES]
TABLES = [config['types'][genome] + config['results']['recalibrated'] + "/" + name + ".table" for name in EXISTINGSAMPLES]
RECBAMS = [config['types'][genome] + config['results']['recalibrated'] + "/" + name + ".bam" for name in EXISTINGSAMPLES]
POSTTABLES = [config['types'][genome] + config['results']['postrecalibrated'] + "/" + name + ".table" for name in EXISTINGSAMPLES]
PDFS = [config['types'][genome] + config['results']['pdfs'] + "/" + name + ".pdf" for name in EXISTINGSAMPLES]
GVCFS = [config['types'][genome] + config['results']['gvcfs'] + "/" + name + ".gvcf" for name in EXISTINGSAMPLES]
GVCFSLIST = ' '.join(["--variant " + config['types'][genome] + config['results']['gvcfs'] + "/" + name + ".gvcf" for name in EXISTINGSAMPLES])

ANNOVARDBS = [config['annovardbdir'] + "/" + genome + "_" + db + ".installed" for db in config['annovardbs']]

ANNOVAR_PROTOCOLS = ','.join(config['annovardbs'])


INDELS = config['types'][genome] + config['results']['realigned'] + "/indels.list"

dlocs = dict()

workdir: config['projdir']

rule all:
    input: 
        trios = TRIOVCFS,
        analysis = ANALYSES,
        phased = config['types'][genome] + config['results']['vcfs'] + "/joint.family.vcf" # must run after all gvcf files created; will create joint.vcf if not already
#include TRIOGEMS for gemini (GRCh37 only)

#this is useful for extracting sequences from bams
rule extract:
    input: EXTRACT

rule rdata:
    input: RDATA

rule xbrowse:
    input: config['types'][genome] + config['results']['vep'] + "/project.yaml", config['types'][genome] + config['results']['vep'] + "/samples.txt", config['types'][genome] + config['results']['vep'] + "/samples.ped"

rule vepvcfs:
    input: VEPVCFS

rule triovcfs:
    input: TRIOVCFS

rule analysisready:
    input: ANALYSISREADY

rule analyses:
    input: ANALYSES

rule Rdeps:
    run:
        R("""
        update.packages(ask=FALSE,repos=c("http://cran.rstudio.com"))
        source("http://bioconductor.org/biocLite.R")
        library(devtools)
        install_github("zhezhangsh/GtUtility")
        install_github("zhezhangsh/awsomics")
        biocLite(c("BSgenome.Hsapiens.UCSC.hg38",
        "TxDb.Hsapiens.UCSC.hg38.knownGene",
        "SNPlocs.Hsapiens.dbSNP.20120608",
        "MafDb.ALL.wgs.phase1.release.v3.20101123",
        "MafDb.ESP6500SI.V2.SSA137",
        "MafDb.ExAC.r0.3.sites",
        "phastCons100way.UCSC.hg19",
        "PolyPhen.Hsapiens.dbSNP131",
        "SNPlocs.Hsapiens.dbSNP144.GRCh38",
        "SIFT.Hsapiens.dbSNP137",
        "org.Hs.eg.db"),suppressUpdates=TRUE)
        install_github("rcastelo/VariantFiltering")
        """)

# this is a utility to put things in the correct order in case something upstream gets touched
rule catchup:
    params:
        picard = config['types'][genome] + config['results']['picard'],
        lists = config['types'][genome] + config['results']['lists'],
        realigned = config['types'][genome] + config['results']['realigned'],
        recalibrated = config['types'][genome] + config['results']['recalibrated'],
        postrecalibrated = config['types'][genome] + config['results']['postrecalibrated'],
        gvcfs = config['types'][genome] + config['results']['gvcfs'],
        vcfs = config['types'][genome] + config['results']['vcfs'],
        analysis = config['types'][genome] + config['results']['analysis']
    shell:
        """
        touch {params.picard}/*rmdup.bam
        touch {params.picard}/*txt
        sleep 2
        touch {params.picard}/*rmdup.bai
        sleep 2
        touch {params.picard}/*group.bam
        sleep 2
        touch {params.picard}/*group.bai
        sleep 2
        touch {params.lists}/*
        sleep 2
        touch {params.realigned}/*bam
        sleep 2
        touch {params.realigned}/*bai
        sleep 2
        touch {params.recalibrated}/*table
        sleep 2
        touch {params.recalibrated}/*bam
        sleep 2
        touch {params.recalibrated}/*bai
        sleep 2
        touch {params.postrecalibrated}/*bam
        sleep 2
        touch {params.postrecalibrated}/*bai
        sleep 2
        touch {params.gvcfs}/*
        sleep 2
        touch {params.vcfs}/*
        sleep 2
        touch {params.analysis}/*
        """

rule sample_concordance:
    output:
        ok="samplesondisk.txt", ms="missingsamples.txt", ump="unmanifestedpairs.txt", ic="incompletefamilies.txt"
    run:
        assert(len(SAMPLELANES) == len(PAIRNAMESINSAMPLETABLE)/2)
        print("Received Pairs (on disk): {0}".format(len(ALLPAIRNAMES)))
        print("Unmanifested Pairs (on disk, not in sample table): {0}".format(len(UNMANIFESTEDPAIRS)))
        f = open(output.ump, 'w')
        for pair in sorted(UNMANIFESTEDPAIRS):
            f.write("{0}\n".format(pair))
        print("Existing Samples (on disk, in sample table): {0}".format(len(EXISTINGSAMPLES)))
        print("Missing Samples (in sample table, not on disk): {0}".format(len(MISSINGSAMPLES)))
        f = open(output.ms, 'w')
        for sample in sorted(MISSINGSAMPLES):
            f.write("{0}\n".format(sample))
        print("Manifested Pairs (in sample table): {0}".format(len(PAIRNAMESINSAMPLETABLE)))
        f = open(output.ok, 'w')
        for sample in sorted(MANIFESTEDPAIRS):
            f.write("{0}\n".format(sample))
        print("Fastqs: {0} Lanes in sample table {1}".format(len(FASTQS), len(PAIRNAMESINSAMPLETABLE)))
        print("Complete Families (trios or quads): {0}".format(len(COMPLETEFAMILYFAMIDS)))
        print("Incomplete Families (files missing): {0}".format(len(INCOMPLETEFAMILIES)))
        f = open(output.ic, 'w')
        for index, row in sample_table.iterrows():
            if row['FamilyID'] in INCOMPLETEFAMILIES:
                f.write("{0}".format(row['FamilyID']))
                for member in ['Mother','Father','Subject']:
                    if row[member] not in EXISTINGSAMPLES and not pandas.isnull(row[member]):
                        f.write("\t{0}".format(row[member]))
                f.write("\n")

rule mkdirs:
    run:
        for adir in config['datadirs']:
            makedir(config['datadirs'][adir])

        for adir in config['types'][genome] + config['results'][genome]:
            makedir(config['datadirs'][genome][adir])

rule dummy:    # just to test the python codes above
    input:  workflow.basedir + "/Snakefile"

    run:
        for file in FASTQCS:
            print(file)
        #check_gvcfs(GVCFS)
        #for file in VEPVCFS:
        #    print(file)
        #for file in FAMILYVCFS:
        #    print(file)
        #for file in SAMPLELANES:
        #    print(file)
        #for file in EXISTINGSAMPLES:
        #    print(file)

rule target_lists:
    input: LISTS

rule print_reads:
    input: RECBAMS

rule join_gvcfs:
    input: config['types'][genome] + config['results']['vcfs'] + "/joint.vcf"

rule make_gvcfs:
    input: GVCFS

rule make_bams:
    input: BAMS

rule make_sams:
    input: SAMS

rule make_merged:
    input: MBAMS

rule pdfs:
    input: PDFS

rule recalibrate:
    input: POSTTABLES

rule combine_indels:
    input: INDELS

rule realign:     # use combined indel list
    input: RBAMS  # realign bams

rule add_group:    # must be run before create_target
    input: GBAIS   # realign bais

rule make_bais:
    input: DBAIS

rule rmdupbams:
    input: DBAMS

rule printbams:
    run:
        print(BAMS)

rule printtrios:
    run:
        print(EXISTINGSAMPLES)
        print("complete trios {0}".format(COMPLETETRIOSFAMIDS))
        print("gems: {0}".format(TRIOGEMS))

rule gzip:
    input: "{sample}.{ext,(fq|fastq)}"
    output: "{sample}.{ext}.gz"
    shell:
        """
        gzip -c {input} > {output}
        """

#### Sequence ####
# make symlinks
# do this once
# remove _001
rule symlinks:
     input: DOWNLOADS
     run:
        if not os.path.exists('fastq'):
            shell("mkdir fastq")
        for file in DOWNLOADS:
            name=os.path.basename(file).replace('_001','')
            if re.search('-L2_', name):
                name = re.sub('-L2_', '_', name)
            else:
                name = re.sub('^E', 'E0', name)
            os.symlink('../' + file, 'fastq/' + name)

rule validateBam:
    input:
        bam = config['datadirs']['oldbams'] + "/{sample}.bam"
    output:
        report = config['datadirs']['oldbams'] + "/{sample}.report.txt"
    shell:
        """
        picard ValidateSamFile I={input.bam} O={output.report}
        """

### Extract reads from older BAM files
rule extractreads:
    input:
        java = ENV3 + config['tools']['java'],
        bam = config['datadirs']['oldbams'] + "/{sample}.bam"
    output:
        pair1 = config['datadirs']['fastq'] + "/{sample}_R1.fastq",
        pair2 = config['datadirs']['fastq'] + "/{sample}_R2.fastq",
        unpaired = config['datadirs']['fastq'] + "/{sample}_U.fastq"
    log:
        config['datadirs']['log'] + "/{sample}.extract.log"
    params:
        picard = config['jars']['picard']['path'],
        md = config['jars']['picard']['samtofastq'],
        opts = config['tools']['opts']['med'],
        metrics = config['types'][genome] + config['results']['picard']
    shell:
        """
        {input.java} {params.opts} -jar {params.picard} \
        {params.md} \
        INPUT={input.bam} \
        FASTQ={output.pair1} \
        SECOND_END_FASTQ={output.pair2} \
        UNPAIRED_FASTQ={output.unpaired} \
        INCLUDE_NON_PF_READS=TRUE \
        VALIDATION_STRINGENCY=LENIENT 2> {log}
        """

### QC ####
rule fastqc: 
    input: 
        pair1 = config['datadirs']['fastq'] + "/{sample}1.fastq.gz",
        pair2 = config['datadirs']['fastq'] + "/{sample}2.fastq.gz",
        seq2qc = ENV3 + config['tools']['seq2qc']
    log: 
        config['datadirs']['log'] + "/{sample}.fastqc.log" 
    output: 
        pair1 = config['datadirs']['fastqc'] + "/{sample}1_fastqc.zip",
        pair2 = config['datadirs']['fastqc'] + "/{sample}2_fastqc.zip",
        html1 = config['datadirs']['fastqc'] + "/{sample}1_fastqc.html",
        html2 = config['datadirs']['fastqc'] + "/{sample}2_fastqc.html"
    params:
        qcdir = config['datadirs']['fastqc']
    # how to run qsub?
    shell: "{input.seq2qc} -o {params.qcdir} {input.pair1} {input.pair2} 2> {log}" 

### Mirror to Cavatica
# waiting on boto fix to accept cbttc bucket name with dots in it
#import configparser
#awsconfig = configparser.ConfigParser()
#awsconfig.read("/home/leipzigj/.aws/credentials")
#from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
#S3 = S3RemoteProvider(access_key_id=awsconfig['default']['aws_access_key_id'], secret_access_key=awsconfig['default']['aws_secret_access_key'], aws_s3_calling_format='boto.s3.connection.OrdinaryCallingFormat')

rule copy_to_cavatica:
    input:
        config['datadirs']['fastq'] + "/{filename}"
    output:
        "aws/{filename}.sent"
        #S3.remote("cbttc.seq.data/Ingo_project/{filename}")
    threads: 1
    shell:
        """
        /home/leipzigj/miniconda3/envs/grinenv/bin/aws s3 cp {input} s3://cbttc.seq.data/Ingo_project/{wildcards.filename}
        touch {output}
        """

rule mirror_cavatica:
    input:
        expand("aws/{filename}_{pair}.fastq.gz.sent",filename=SAMPLELANES,pair=['R1','R2'])
#        S3.remote(expand("cbttc.seq.data/Ingo_project/{filename}.fastq.gz",filename=SAMPLELANES))

rule copy_to_risaws:
     input:
        config['datadirs']['fastq'] + "/{filename}"
     output:
        "risaws/{filename}.sent"
     threads: 1
     shell:
        """
        /home/leipzigj/miniconda3/envs/grinenv/bin/aws s3 cp s3://wuxi-demo-trios.s3.amazonaws.com/{wildcards.filename}
        touch {output}
        """

RISAWSNAMES = ['E01621','E01623','E01622']
RISAWSFASTQ =[s for s in SAMPLELANES for r in RISAWSNAMES if r in s]

rule mirror_risaws:
     input:
        expand("risaws/{filename}_{pair}.fastq.gz.sent",filename=RISAWSFASTQ,pair=['R1','R2'])

#### Alignment ####
rule align:
    input:
        pair1 = config['datadirs']['fastq'] + "/{sample}_R1.fastq.gz",
        pair2 = config['datadirs']['fastq'] + "/{sample}_R2.fastq.gz",
        align = ENV3 + config['tools']['align']
    output:
        sam = config['types'][genome] + config['results']['sams'] + "/{sample}.sam" # may be set to temp
    threads:
        12
    log: 
        config['datadirs']['log'] + "/{sample}.novoalign.log"
    params:
        refidx = config['refidx'][genome]
    shell:
        """
        {input.align} -c {threads} -a -k -d {params.refidx} -o SAM -f {input.pair1} {input.pair2} 1> {output.sam} 2> {log}
        """

rule sam_to_bam:
    input:
        sam = config['types'][genome] + config['results']['sams'] + "/{sample}.sam",
        samtools = ENV3 + config['tools']['samtools']
    output:
        bam = config['types'][genome] + config['results']['bams'] + "/{sample,[^.]+}.bam"
    threads:
        12   # also depends on -j
    shell:
        """
        {input.samtools} view -@ {threads} -bS {input.sam} > {output.bam}
        """

# novosort creates index
rule novosortbam:
    input:
        bam = config['types'][genome] + config['results']['bams'] + "/{sample}.bam",
        sort = ENV3 + config['tools']['sortbam']
    output:
        sorted = config['types'][genome] + config['results']['bams'] + "/{sample}.sorted.bam",
    threads:
        12
    shell:
        """
        {input.sort} -m 14g -t . --removeduplicates --keeptags -i -o {output.sorted} {input.bam}
        """

## here we create individual realign target lists, one from each bam file
## later, the individual lists will be combined into a single list
## an alternative ways is to directly create a signle list from the bam files

rule target_list: # create individual realign target list
    input:  # deduced bams
        bai = config['types'][genome] + config['results']['picard'] + "/{sample}.group.bai",
        bam = config['types'][genome] + config['results']['picard'] + "/{sample}.group.bam",
        java = ENV3 + config['tools']['java']
    output:
        samplelist = config['types'][genome] + config['results']['lists'] + "/{sample}.list"
    log:
        config['datadirs']['log'] + "/{sample}.target_list.log"
    params:
        jar = config['jars']['gatk'],
        opts = config['tools']['opts']['high'],
        ref = config['ref'][genome],
        knownsites = config['known'][genome]
    #threads:
    #    24
    shell:
        """
        {input.java} {params.opts} -jar {params.jar} \
        -T RealignerTargetCreator \
        -R {params.ref} \
        -I {input.bam} \
        -known {params.knownsites} \
        -o {output.samplelist} 2> {log}
        """
        # -nt {threads} \

def makedir(adir):
    if not os.path.exists(adir):
        shell("mkdir -p " + adir)

def cmp(a,b):
    global dlocs
    if dlocs[a][0] == dlocs[b][0]:
       if dlocs[a][1] > dlocs[b][1]:
           return 1
       else:
           return -1
    else:
       if dlocs[a][0] > dlocs[b][0]:
           return 1
       else:
           return -1

def combine():
    global dlocs
    "combine indels in the individual lists into a single list"
    indels = dict()
    for name in LISTS:
        fin = open(name, "r")
        for line in fin.readlines():
            line = re.sub('\n', '', line)
            chr, loc = line.split(':')
            if chr in indels:
                if loc not in indels[chr]:
                    locs = loc.split('-')
                    if len(locs) == 1:
                        locs.append(locs[0])
                    indels[chr][loc] = [int(locs[0]), int(locs[1])]
            else:
                indels[chr] = dict()
                locs = loc.split('-')
                if len(locs) == 1:
                    locs.append(locs[0])
                indels[chr][loc] = [int(locs[0]), int(locs[1])]
                    
    # how to sniff version?
    chrs = list(indels.keys())  # works for both version 2 and 3
    chrs.sort()
    with open(INDELS, "w") as out:
        for chr in chrs:
            dlocs = indels[chr]
            ks = list(dlocs.keys())
            # ks.sort(key=lamda x:dlocs[x[0])
            ks.sort(key=cmp_to_key(cmp))
            for loc in ks:
                out.write(chr + ':' + loc + "\n") 

# combine indels from individual target alignment lists into a single list
rule combine_lists:
    input: LISTS
    output: INDELS
    run:
         combine()

rule realign_target:   # with one combined list file
    input:  # deduced bams
        #list = INDELS,
        list = config['types'][genome] + config['results']['lists'] + "/{sample}.list",
        dbam = config['types'][genome] + config['results']['picard'] + "/{sample}.group.bai",
        java = ENV3 + config['tools']['java']
    output:
        rbam = config['types'][genome] + config['results']['realigned'] + "/{sample}.bam"
    params:
        jar = config['jars']['gatk'],
        opts = config['tools']['opts']['med'],
        ref = config['ref'][genome],
        known = config['known'][genome],
        result = config['types'][genome] + config['results']['picard']
    shell:
        """
        {input.java} {params.opts} -jar {params.jar} \
        -T IndelRealigner \
        -R {params.ref} \
        -I {params.result}/{wildcards.sample}.group.bam \
        -targetIntervals {input.list} \
        -known {params.known} \
        -o {output.rbam}
        """

# Base recalibration (not be confused with variant recalibration)
# http://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr
# https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php
rule generate_recalibration_table:
    input:
        bam = config['types'][genome] + config['results']['realigned'] + "/{sample}.bam",
        java = ENV3 + config['tools']['java']
    output:
        table = config['types'][genome] + config['results']['recalibrated'] + "/{sample}.table"
    log:
        config['datadirs']['log'] + "/{sample}.generate_recalibration_table.log"
    params:
        jar = config['jars']['gatk'],
        opts = config['tools']['opts']['med'],
        ref = config['ref'][genome],
        known = config['known'][genome]
    shell:
        """
        {input.java} {params.opts} -jar {params.jar} \
        -T BaseRecalibrator \
        -R {params.ref} \
        -I {input.bam} \
        -knownSites {params.known} \
        -o {output.table} 2> {log}
        """

# https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_readutils_PrintReads.php
# Note that when PrintReads is used as part of the Base Quality Score Recalibration workflow
rule recalibrate_bam:
    input:
        table = config['types'][genome] + config['results']['recalibrated'] + "/{sample}.table",
        bam = config['types'][genome] + config['results']['realigned'] + "/{sample}.bam",
        java = ENV3 + config['tools']['java']
    output:
        bam = config['types'][genome] + config['results']['recalibrated'] + "/{sample}.bam"
    log:
        config['datadirs']['log'] + "/{sample}.recalibrate_bam.log"
    params:
        jar = config['jars']['gatk'],
        opts = config['tools']['opts']['med'],
        ref = config['ref'][genome]
    threads:
        8
    shell:
        """
        {input.java} {params.opts} -jar {params.jar} \
        -T PrintReads \
        -nct {threads} \
        -R {params.ref} \
        -I {input.bam} \
        -BQSR {input.table} \
        -o {output.bam} 2> {log}
        """

rule post_recalibrated_table:
    input:
        table = config['types'][genome] + config['results']['recalibrated'] + "/{sample}.table",
        bam = config['types'][genome] + config['results']['realigned'] + "/{sample}.bam",
        java = ENV3 + config['tools']['java']
    output:
        table = config['types'][genome] + config['results']['postrecalibrated'] + "/{sample}.table",
    params:
        jar = config['jars']['gatk'],
        opts = config['tools']['opts']['med'],
        ref = config['ref'][genome],
        known = config['known'][genome]
    shell:
        """
        {input.java} {params.opts} -jar {params.jar} \
        -T BaseRecalibrator \
        -R {params.ref} \
        -I {input.bam} \
        -knownSites {params.known} \
        -BQSR {input.table} \
        -o {output.table}
        """

# https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_bqsr_AnalyzeCovariates.php
# generates plots for visualizing the quality of a recalibration run
rule analyze_bqsr:
    input:
        before = config['types'][genome] + config['results']['recalibrated'] + "/{sample}.table",
        after = config['types'][genome] + config['results']['postrecalibrated'] + "/{sample}.table",
        java = ENV3 + config['tools']['java']
    output:
        pdf = config['types'][genome] + config['results']['pdfs'] + "/{sample}.pdf"
    params:
        jar = config['jars']['gatk'],
        opts = config['tools']['opts']['med'],
        ref = config['ref'][genome]
    shell:
        """
        {input.java} {params.opts} -jar {params.jar} \
        -T AnalyzeCovariates \
        -R {params.ref} \
        -before {input.before} \
        -after {input.after} \
        -plots {output.pdf}
        """

# merge lanes
# E01188-L2_S26_L005.sorted.bam E01188-L2_S26_L006.sorted.bam > E01188.sorted.merged.bam
def get_all_sorted_bams(samplename):
    bams = [bam for bam in SBAMS if bam.startswith(samplename) and bam.endswith(".sorted.bam")]
    return(bams)

rule merge_lanes:
    input: bams = lambda wildcards: get_all_sorted_bams(wildcards.sample), samtools = ENV3 + config['tools']['samtools']
    output: "{sample}.sorted.merged.bam"
    threads:
        1
    run:
        if len(input.bams)>1:
            shell("{input.samtools} merge {output} {input.bams}")
        else:
            shell("cp {input.bams} {output}")

rule depth_of_coverage:
    input:
        bam = config['types'][genome] + config['results']['recalibrated'] + "/{sample}.bam",
        java = ENV3 + config['tools']['java']
    output:
        "{sample}.DoC"
    params:
        jar = config['jars']['gatk'],
        opts = config['tools']['opts']['med'],
        ref = config['ref'][genome]
    shell:
        """
        {input.java} {params.opts} -jar {params.jar} \
        -T DepthOfCoverage \
        -I {input.bam} \
        -R {params.ref} \
        -L {params.targets} \
        -l INFO \
        -dt BY_SAMPLE \
        --omitDepthOutputAtEachBase \
        --omitLocusTable \
        --minBaseQuality 0 \
        --minMappingQuality 20 \
        --start 1 \
        --stop 5000 \
        --nBins 200 \
        --includeRefNSites \
        -o {output}
        """

rule mark_duplicates:
    input:
        bam = config['types'][genome] + config['results']['bams'] + "/{sample}.sorted.merged.bam",
        java = ENV3 + config['tools']['java']
    output:
        bam = config['types'][genome] + config['results']['picard'] + "/{sample}.rmdup.bam"
    log:
        config['datadirs']['log'] + "/{sample}.markdups.log"
    params:
        picard = config['jars']['picard']['path'],
        md = config['jars']['picard']['markdups'],
        opts = config['tools']['opts']['med'],
        metrics = config['types'][genome] + config['results']['picard']
    shell:
        # will (and need the permision to) create a tmp directory
        # with the name of login under specified tmp directory
        # Exception in thread "main" net.sf.picard.PicardException: Exception creating temporary directory.
        """
        {input.java} {params.opts} -jar {params.picard} \
        {params.md} \
        INPUT={input.bam} \
        OUTPUT={output.bam} \
        METRICS_FILE={params.metrics}/{wildcards.sample}.txt 2> {log}
        """

# samtools index seems more reliable than picard
# in terms of returning an exit code
rule make_index:
    input:
        bam = config['types'][genome] + config['results']['picard'] + "/{sampleandext}.bam",
        samtools = ENV3 + config['tools']['samtools']
    output:
        bai = config['types'][genome] + config['results']['picard'] + "/{sampleandext}.bai"
    shell:
        """
        {input.samtools} index {input.bam} {output.bai}
        """


rule add_readgroup:
    input:
        bam = config['types'][genome] + config['results']['picard'] + "/{sample}.rmdup.bam",
        bai = config['types'][genome] + config['results']['picard'] + "/{sample}.rmdup.bai",
        java = ENV3 + config['tools']['java']
    output:
        bam = config['types'][genome] + config['results']['picard'] + "/{sample}.group.bam"
    log:
        config['datadirs']['log'] + "/{sample}.add_readgroup.log"
    params:
        picard = config['jars']['picard']['path'],
        opts = config['tools']['opts']['med'],
        rg = config['jars']['picard']['readgroups']
    shell:
        """
        {input.java} {params.opts} -jar {params.picard} \
        {params.rg} \
        I={input.bam} \
        O={output.bam} \
        PL=illumina \
        LB={wildcards.sample} \
        PU={wildcards.sample} \
        SM={wildcards.sample} 2> {log}
        """

#### Variant Calling ####
rule make_gvcf:
    input:
        bam = config['types'][genome] + config['results']['recalibrated'] + "/{sample}.bam",
        java = ENV3 + config['tools']['java']
    output:
        gvcf = config['types'][genome] + config['results']['gvcfs'] + "/{sample}.gvcf"
    params:
        jar = config['jars']['gatk'],
        opts = config['tools']['opts']['med'],
        ref = config['ref'][genome]
    shell:
        """
        {input.java} {params.opts} -jar {params.jar} \
        -T HaplotypeCaller \
        -R {params.ref} \
        -I {input.bam} \
        --emitRefConfidence GVCF \
        --variant_index_type LINEAR \
        --variant_index_parameter 128000 \
        --genotyping_mode DISCOVERY \
        -stand_emit_conf 10 \
        -stand_call_conf 30 \
        -o {output.gvcf}
        """

# sanity check on gvcf files

def check_gvcfs(gvcfs):
    for file in gvcfs:
        cmd = "tail -1 " + file
        # print(cmd)
        (output, error) = call_command(cmd)
        line = output.decode()
        if not re.search("^chrEBV.*\n", line):
            print(file + ' is incomplete')
            quit()
    
def call_command(command):
    process = subprocess.Popen(command.split(' '),
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    return process.communicate()

# if supplid a subject return either only one trio for that subject
# if supplied 'family' as the subject, then all members, including siblings are returned
# the family argument is really for sanity checking if a subject is given
def gvcf_samples_in_family(family,subject):
    # joint means all existing samples
    if family == 'joint':
        return [GVCFS,GVCFSLIST]
    else:
        if subject == 'family':
            rows = sample_table.loc[(sample_table['FamilyID'] == family)]
        else:
            rows = sample_table.loc[(sample_table['FamilyID'] == family) & (sample_table['Subject'] == subject)]
        samples = list(rows['Subject'].dropna())+list(rows['Mother'].dropna())+list(rows['Father'].dropna())
        assert(len(samples)>0)
        if subject != 'family':
            assert(len(samples)==3)
        gvcfs = [config['types'][genome] + config['results']['gvcfs'] + "/" + name + ".gvcf" for name in set(samples)]
        gvcfslist = ' '.join(["--variant " + config['types'][genome] + config['results']['gvcfs'] + "/" + name + ".gvcf" for name in set(samples)])
        return [gvcfs, gvcfslist]

# make sure the family is done first
# if family size is 3, then you're done, just copy
# otherwise make both trios
rule trio_vcfs:
    input:
        gvcfs = lambda wildcards: gvcf_samples_in_family(wildcards.family,wildcards.subject)[0],
        family = config['types'][genome] + config['results']['vcfs'] + "/{family}.family.vcf",
        familygvcfs = lambda wildcards: gvcf_samples_in_family(wildcards.family,'family')[0],
        java = ENV3 + config['tools']['java']
    output:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{family}_{subject}.trio.vcf"
    log:
        config['datadirs']['log'] + "/{family}_{subject}.trio.vcf.log"
    params:
        jar = config['jars']['gatk'],
        ref = config['ref'][genome],
        gvcfslist = lambda wildcards: gvcf_samples_in_family(wildcards.family,wildcards.subject)[1],
        opts = config['tools']['opts']['med'],
        db = config['dbsnp']
    threads: 8
    run:
        assert(len(input.gvcfs)==3)
        check_gvcfs(input.gvcfs) 
        check_gvcfs(input.familygvcfs)
        # use the family count to determine course of action
        if len(input.familygvcfs)==3:
            shell("cp {input.family} {output.vcf}")
        else:
            shell("""
                {input.java} {params.opts} -jar {params.jar} \
                -T GenotypeGVCFs \
                --disable_auto_index_creation_and_locking_when_reading_rods \
                --dbsnp {params.db} \
                -nt {threads} \
                -R {params.ref} \
                {params.gvcfslist} \
                -o {output.vcf} 2> {log}
            """)

rule family_vcfs:
    input:
        gvcfs = lambda wildcards: gvcf_samples_in_family(wildcards.family,'family')[0],
        java = ENV3 + config['tools']['java']
    output:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{family}.family.vcf",
        idx = config['types'][genome] + config['results']['vcfs'] + "/{family}.family.vcf.idx"
    log:
        config['datadirs']['log'] + "/{family}.family.vcf.log"
    params:
        jar = config['jars']['gatk'],
        ref = config['ref'][genome],
        gvcfslist = lambda wildcards: gvcf_samples_in_family(wildcards.family,'family')[1],
        opts = config['tools']['opts']['med'],
        db = config['dbsnp']
    threads: 8
    run:
        check_gvcfs(input.gvcfs) 
        shell("""
        {input.java} {params.opts} -jar {params.jar} \
        -T GenotypeGVCFs \
        --disable_auto_index_creation_and_locking_when_reading_rods \
        --dbsnp {params.db} \
        -nt {threads} \
        -R {params.ref} \
        {params.gvcfslist} \
        -o {output.vcf} 2> {log}
        """)

def samples_in_family(family):
    rows = sample_table.loc[sample_table['FamilyID'] ==  family]
    for row in rows:
        samples += [row['Subject'],row['Mother'],row['Father']]

# convert the GRIN sample table into a GATK compliant 6-column pedfile:
# http://gatkforums.broadinstitute.org/gatk/discussion/37/pedigree-analysis
# For these tools, the PED files must contain only the first 6 columns from the PLINK format PED file, and no alleles, like a FAM file in PLINK.
# Family ID
# Individual ID
# Paternal ID
# Maternal ID
# Sex (1=male; 2=female; other=unknown)
# Phenotype
rule sample_table_to_pedfile:
    input:
        config['sample_table']
    params:
        existingsamples = EXISTINGSAMPLES
    output:
        config['pedfile']
    run:
        st = pandas.read_table("{0}".format(input))
        complete_trios = [index for index, row in st.iterrows() if all([row[member] in params.existingsamples for member in ['Mother','Father','Subject']])]
        trios = st.loc[complete_trios]
        trios['Sex']=trios['Sex'].replace(['M','F'],[1,2])
        trios['Affected_status']=trios['Affected_status'].replace('unaffected',1)
        trios['Affected_status']=trios['Affected_status'].replace('[^1].+', 2, regex=True)
        ped = trios[[0,1,3,2,4,5]]
        ped = ped.fillna(0)
        moms = [mom for mom in list(ped['Mother'].dropna()) if mom not in list(ped['Subject'])]
        dads = [dad for dad in list(ped['Father'].dropna()) if dad not in list(ped['Subject'])]
        #add rows for unaffected parents
        momfams=ped[ped['Mother'].isin(moms)]['FamilyID']
        momdf = pandas.DataFrame(momfams)
        momdf['Subject']=moms
        momdf['Father']=0
        momdf['Mother']=0
        momdf['Sex']=2
        momdf['Affected_status']=1
        
        dadfams=ped[ped['Father'].isin(dads)]['FamilyID']
        daddf = pandas.DataFrame(dadfams)
        daddf['Subject']=dads
        daddf['Father']=0
        daddf['Mother']=0
        daddf['Sex']=1
        daddf['Affected_status']=1
        
        ped = ped.append([momdf,daddf])
        ped = ped.drop_duplicates()
        ped = ped.sort_values(by=['FamilyID','Subject'])
        ped.to_csv("{0}".format(output), sep='\t',index=False)

# R's Varia
rule analysis_pedfile:
    input:
        config['pedfile']
    output:
        config['types'][genome] + config['results']['analysis'] + "/{family}_{subject}.pedfile"
    run:
        globalpedfile = pandas.read_table("{0}".format(input))
        
        probandrow = globalpedfile[(globalpedfile['FamilyID'] == wildcards.family) & (globalpedfile['Subject'] == wildcards.subject)]
        assert(len(probandrow)==1)
        parentalrows = globalpedfile[(globalpedfile['Subject'] == probandrow['Father'].item()) | (globalpedfile['Subject'] == probandrow['Mother'].item())]
        assert(len(parentalrows)==2)
        ped = probandrow.append([parentalrows])
        
        # hashes in the family names confuse R
        ped = ped.replace('#','',regex=True)
        ped.to_csv("{0}".format(output), sep='\t',index=False)
        
        
rule run_phase_by_transmission:
    input:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.trio.vcf",
        ped = config['pedfile'],
        java = ENV3 + config['tools']['java']
    output:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.trio.phased.vcf",
        idx = config['types'][genome] + config['results']['vcfs'] + "/{file}.trio.phased.vcf.idx",
        mvf = config['types'][genome] + config['results']['vcfs'] + "/{file}.mendelian_violations.txt"
    params:
        jar  = config['jars']['gatk'],
        ref = config['ref'][genome],
        opts = config['tools']['opts']['med']
    log: 
        config['datadirs']['log'] + "/{file}.phase_by_transmission.log" 
    shell:
        """
        {input.java} {params.opts} -jar {params.jar} \
        -T PhaseByTransmission \
        -R {params.ref} \
        -V {input.vcf} \
        -ped {input.ped} \
        -mvf {output.mvf} \
        -o {output.vcf} 2> {log}
        """
#        --disable_auto_index_creation_and_locking_when_reading_rods \

# https://www.broadinstitute.org/gatk/guide/article?id=2806
# https://github.com/chapmanb/bcbio-nextgen/blob/master/bcbio/variation/vfilter.py

rule gatk_snps_only:
    input:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.vcf",
        java = ENV3 + config['tools']['java']
    output:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.snps.vcf",
        idx = config['types'][genome] + config['results']['vcfs'] + "/{file}.snps.vcf.idx"
    params:
        jar  = config['jars']['gatk'],
        ref = config['ref'][genome],
        opts = config['tools']['opts']['low']
    log:
        config['datadirs']['log'] + "/{file}.gatk_snps_only.log"
    shell:
        """
        {input.java} {params.opts} -jar {params.jar} \
        -T SelectVariants \
        -R {params.ref} \
        -V {input.vcf} \
        -selectType SNP \
        -o {output.vcf} 2> {log}
        """

rule gatk_indels_only:
    input:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.vcf",
        java = ENV3 + config['tools']['java']
    output:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.indels.vcf",
        idx = config['types'][genome] + config['results']['vcfs'] + "/{file}.indels.vcf.idx"
    params:
        jar  = config['jars']['gatk'],
        ref = config['ref'][genome],
        opts = config['tools']['opts']['low']
    log:
        config['datadirs']['log'] + "/{file}.gatk_indels_only.log"
    shell:
        """
        {input.java} {params.opts} -jar {params.jar} \
        -T SelectVariants \
        -R {params.ref} \
        -V {input.vcf} \
        -selectType INDEL \
        -o {output.vcf} 2> {log}
        """
        
# hard filtration
# this "filters out, not filters for" filterExpression
rule gatk_hard_filtration_snps:
    input:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.snps.vcf",
        java = ENV3 + config['tools']['java']
    output:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.snps.hard.vcf",
        idx = config['types'][genome] + config['results']['vcfs'] + "/{file}.snps.hard.vcf.idx"
    params:
        jar  = config['jars']['gatk'],
        ref = config['ref'][genome],
        opts = config['tools']['opts']['low']
    log:
        "log/{file}.gatk_hard_filtration.log"
    shell:
        """
        {input.java} {params.opts} -jar {params.jar} \
        -R {params.ref} \
        -T VariantFiltration \
        -o {output.vcf} \
        --variant {input.vcf} \
        --filterExpression "QD < 2.0 || MQ < 30.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
        --filterName "GATK3.5-hard-filter" \
        2> {log}
        """

rule gatk_hard_filtration_indels:
    input:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.indels.vcf",
        java = ENV3 + config['tools']['java']
    output:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.indels.hard.vcf",
        idx = config['types'][genome] + config['results']['vcfs'] + "/{file}.indels.hard.vcf.idx"
    params:
        jar  = config['jars']['gatk'],
        ref = config['ref'][genome],
        opts = config['tools']['opts']['low']
    log:
        "log/{file}.gatk_hard_filtration.log"
    shell:
        """
        {input.java} {params.opts} -jar {params.jar} \
        -R {params.ref} \
        -T VariantFiltration \
        -o {output.vcf} \
        --variant {input.vcf} \
        --filterExpression "QD < 2.0 || ReadPosRankSum < -20.0 || FS > 200.0" \
        --filterName "GATK3.5-hard-filter" \
        2> {log}
        """

rule select_passing:
    input:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.{type,(snps|indels)}.hard.vcf",
        java = ENV3 + config['tools']['java']
    output:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.{type,(snps|indels)}.filtered.vcf",
        idx = config['types'][genome] + config['results']['vcfs'] + "/{file}.{type,(snps|indels)}.filtered.vcf"
    params:
        jar  = config['jars']['gatk'],
        ref = config['ref'][genome],
        opts = config['tools']['opts']['low']
    log:
        "log/{file}.select_passing_variants.log"
    shell:
        """
        {input.java} {params.opts} -jar {params.jar} \
        -R {params.ref} \
         -T SelectVariants \
        -o {output.vcf} \
        --variant {input.vcf} \
        --excludeFiltered \
        2> {log}
        """

# """Run GATK CombineVariants to combine variant files.
#
# The default rule combines files with suffixes filteredSNP.vcf and
# filteredINDEL.vcf.
# """
rule gatk_combine_variants:
    input:
        snps = config['types'][genome] + config['results']['vcfs'] + "/{file}.snps.filtered.vcf",
        indels = config['types'][genome] + config['results']['vcfs'] + "/{file}.indels.filtered.vcf",
        java = ENV3 + config['tools']['java']
    output:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.com.filtered.vcf",
        idx = config['types'][genome] + config['results']['vcfs'] + "/{file}.com.filtered.vcf.idx"
    params:
        jar  = config['jars']['gatk'],
        ref = config['ref'][genome],
        opts = config['tools']['opts']['med']
    log:
        "log/{file}.select_passing_variants.log"
    shell:
        """
        {input.java} {params.opts} -jar {params.jar} \
        -R {params.ref} \
        -T CombineVariants \
        --variant  {input.snps} \
        --variant  {input.indels} \
        -o {output.vcf} \
        --assumeIdenticalSamples \
        2> {log}
        """

#### Annotation ####
rule ad_vcf:
    input:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.vcf",
    output:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.ad.vcf"
    params:
        ref = config['ref'][genome]
    shell:
        """
        cat {input} | sed 's/ID=AD,Number=./ID=AD,Number=R/' > {output}
        """

# decomposes multiallelic variants into biallelic in a VCF file.
rule decompose_for_gemini:
    input:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.ad.vcf",
        vt = ENV3 + config['tools']['vt']
    output:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.ad.de.vcf"
    shell:
        """
        {input.vt} decompose -s -o {output} {input.vcf}
        """

rule normalize_for_gemini:
    input:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.ad.de.vcf",
        vt = ENV3 + config['tools']['vt']
    output:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.ad.de.nm.vcf"
    params:
        ref = config['ref'][genome]
    shell:
        """
        {input.vt} normalize -r {params.ref} -o {output} {input.vcf}
        """

rule vcf_qt:
    input:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.ad.de.vcf",
        vt = ENV3 + config['tools']['vt']
    output:
        vcf = config['types'][genome] + config['results']['vtpeek'] + "/{file}.vtpeek.txt"
    params:
        ref = config['ref'][genome]
    shell:
        """
        {input.vt} peek -o {output} {input.vcf}
        """

rule vcf_profile:
    input:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.ad.de.vcf",
        vt = ENV3 + config['tools']['vt'],
        ped = config['pedfile']
    output:
        vcf = config['types'][genome] + config['results']['vtpeek'] + "/{file}.vtmendelprofile.txt"
    params:
        ref = config['ref'][genome]
    shell:
        """
        {input.vt} profile_mendelian -o {output} -p {input.ped} -x mendel {input.vcf}
        """
        
# ud - upstream downstream interval length (in bases)
rule run_snpeff:
    input:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.vcf",
    output:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.snpeff.vcf"
    params:
        snpeff  = config['jars']['snpeff']['path'],
        conf = config['jars']['snpeff']['cnf'],
        opts = config['tools']['opts']['med'],
        database = config['jars']['snpeff']['db'],
        updown = config['jars']['snpeff']['ud'],
        format = config['jars']['snpeff']['format']
    shell:
        """
        source /home/leipzigj/miniconda3/envs/grinenv/bin/activate grinenv 
        {params.snpeff} \
        {params.opts} \
        -c {params.conf} \
        -ud {params.updown} \
        {params.format} \
        {params.database} \
         {input.vcf} > {output.vcf}
        """

#### run VEP  ####

rule for_xbrowse:
    input: VEPVCFS
    output:
         yaml = config['types'][genome] + config['results']['vep'] + "/project.yaml",
         list = config['types'][genome] + config['results']['vep'] + "/samples.txt",
         ped = config['types'][genome] + config['results']['vep'] + "/samples.ped"
    params:
         pedfile = config['pedfile']
    run:
        with open(output.yaml, "w") as out:
            out.write("---\n\n") 
            out.write("project_id: '%s'\n" % (config['xbrowse']['id']))
            out.write("project_name: '%s'\n" % (config['xbrowse']['name']))
            out.write("sample_id_list: 'samples.txt'\n")
            out.write("ped_files:\n")
            out.write("  - 'samples.ped'\n")
            out.write("vcf_files:\n")
            for name in VEPVCFS:
                out.write("  - '%s'\n" % (re.sub('.*\/', '', name)))

        with open(output.list, "w") as out:
            for name in sorted(EXISTINGSAMPLES):
                out.write(name + "\n")

        with open(output.ped, "w") as out:
            fin = open(params.pedfile, "r")
            for line in fin.readlines():
                fields = line.split()
                # print(fields)
                if fields[1] in EXISTINGSAMPLES:
                    out.write(line)

rule run_vep:
    input:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{family}.vcf",
        vep = config['tools']['vep']
    output:
        vep = config['types'][genome] + config['results']['vep'] + "/{family}.vep.vcf"
    params:
        xbrowse = config['xbrowse'],
        vepdir = config['vepdir'],
        vepgen = config['vepgenomes'][genome],
        assgen = config['vep'][genome],
    run:
        cmd = "perl " + input.vep + " \
          --everything --vcf --allele_number --no_stats --cache --offline \
          --force_overwrite --cache_version 84 \
          --dir " + params.vepdir + " \
          --fasta " + params.vepgen + " \
          --assembly " + params.assgen + " \
          --plugin LoF,human_ancestor_fa:" + params.xbrowse + "/data/reference_data/human_ancestor.fa.gz,filter_position:0.05... \
          --plugin dbNSFP," + params.xbrowse + "/data/reference_data/dbNSFP.gz,Polyphen2_HVAR_pred,CADD_phred,SIFT_pred,FATHMM_pred,MutationTaster_pred,MetaSVM_pred \
              -i " + input.vcf + " -o " + output.vep
        shell(cmd)

#### run multiqc  ####

rule run_multiqc:
    input:
        gbams = GBAMS,
        multiqc = ENV3 + config['tools']['multiqc']
    output: config['datadirs']['multiqc'] + '/multiqc_report.html'
    params:
        dirs = config['types'][genome] + config['results']['picard'] + ' fastqc',
        outdir = config['datadirs']['multiqc'] 
    shell:
        """
        {input.multiqc} -o {params.outdir} {params.dirs} # will detect input file types?
        """

#### run annovar  ####

rule table_annovar:
    input:
        ANNOVARDBS,
        avinput = config['types'][genome] + config['results']['vcfs'] + "/{file}.vcf",
        annovar = ENV3 + config['tools']['table_annovar']
    output:
        config['types'][genome] + config['results']['vcfs'] + "/{file}.annovar.vcf"
    params:
        opts = "-buildver "+genome
                +" -protocol "+ANNOVAR_PROTOCOLS
                +" -operation "+config['operations']
                +" -nastring . \
                -out joint \
                -tempdir {config['tmpdir']} \
                -remove \
                -dot2underline \
                -vcfinput",
        dbdir = config['annovardbdir']
    shell:
        """
        {input.annovar} {params.opts} {input.avinput} {params.dbdir}
        """

rule run_annovar:
    input:
        ANNOVARDBS,
        avinput = config['types'][genome] + config['results']['vcfs'] + "/{file}.avinput",
        annovar = ENV3 + config['tools']['annotate_variation']
    output:
        config['types'][genome] + config['results']['vcfs'] + "/annovar.done"
    params:
        opts = "-buildver {genome}",
        dbdir = config['annovardbdir']
    run:
        # gene based annotation
        shell("{input.annovar} -geneanno {params.opts} {input.avinput} {params.dbdir}")

        for db in config['annovardbs']:

            # region based annotation
            shell("{input.annovar} -regionanno -dbtype {db} {params.opts} {input.avinput} {params.dbdir}")

            # filter based annotation
            shell("{input.annovar} -filter -dbtype {db} {params.opts} {input.avinput} {params.dbdir}")

        shell("touch {output}")

rule vcf2avinput:
    input:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.vcf",
        cmd = ENV3 + config['tools']['vcf2avinput']
    output:
        config['types'][genome] + config['results']['vcfs'] + "/{file}.avinput",
    shell:
        "{input.cmd} -format vcf2old {input.vcf} -outfile {output}"

rule install_annovar_db:
    input:
        annovar = ENV3 + config['tools']['annotate_variation']
    output:
        config['annovardbdir'] + "/{genome}_{db}.installed"
    params:
        dbdir = config['annovardbdir'],
    run:
        opts = config['annovaropts'][wildcards.genome][wildcards.db]
        if wildcards.db == 'ALL.sites.2014_10':
            shell("{input.annovar} -buildver {wildcards.genome} {opts} 1000g2014oct {params.dbdir}")
            shell("unzip -d {params.dbdir} {params.dbdir}/{wildcards.genome}_1000g2014oct.zip {wildcards.genome}_{wildcards.db}.txt")
        else:
            shell("{input.annovar} -buildver {wildcards.genome} {opts} {wildcards.db} {params.dbdir}")
        shell("touch {output}")


rule compress_vcf:
    input:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.vcf",
        bgzip = ENV3 + config['tools']['bgzip']
    output:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.vcf.bgz",
    shell:
        """
        {input.bgzip} -c {input.vcf} > {output}
        """

rule tabix:
    input:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.vcf.bgz",
        tabix = ENV3 + config['tools']['tabix']
    output:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.vcf.bgz.tbi",
    shell:
        """
        {input.tabix} -p vcf {input.vcf}
        """

rule gemini_db:
    input:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.trio.phased.com.filtered.ad.de.nm.snpeff.vcf.bgz",
        tbi = config['types'][genome] + config['results']['vcfs'] + "/{file}.trio.phased.com.filtered.ad.de.nm.snpeff.vcf.bgz.tbi",
        ped = config['pedfile'],
        gemini = ENV2 + config['tools']['gemini']
    output:
        config['types'][genome] + config['results']['gemini'] + "/{file}.gemini.db"
    threads:
        3
    shell:
        """
        {ENV2}
        {input.gemini} load --cores {threads} -t snpEff -v {input.vcf} -p {input.ped} {output}
        """

rule testR:
    run:
        #import rpy2.robjects as robjects
        #print(robjects.r)
        (output, error) = call_command("which R")
        print("{0} {1}".format(output,error))
        try:
            import rpy2
            print("rpy2 version is", rpy2.__version__)
        except Exception as e:
            print("no rpy2", e)
        #print("")
        #R("""
        #rnorm(100)
        #""")

rule describeR:
    run:
        R("""
        library(dplyr)
        library(VariantFiltering)
        print(sessionInfo())
        """)

# alter records that have a '*' in the ALT field - these are causing problems in VariantAnnotation
# https://github.com/Bioconductor-mirror/VariantAnnotation/commit/837f1f0c9fdcdc9de7677a9a067963413dfe26e7
rule noasterisk:
    input:
        vcf = "{file}.vcf"
    output:
        vcf = "{file}.noask.vcf"
    shell:
        """
        cat {input} | sed -e 's/\*/N/g' > {output}
        """


#trio.phased.com.filtered.ad.de.nm.snpeff.vcf.bgz
#### Analysis ####
rule variantAnalysisSetupUind:
    input:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{familypro}.vcf.bgz",
    output:
        uind = config['types'][genome] + config['results']['analysis'] + "/{familypro}.uind.RData",
    params:
        bsgenome = config['analysis']['bsgenome'][genome],
        txdb     = config['analysis']['txdb'][genome],
        snpdb    = config['analysis']['snpdb'][genome],
        esp      = config['analysis']['esp'][genome],
        exac     = config['analysis']['exac'][genome],
        sift     = config['analysis']['sift'][genome],
        phylo    = config['analysis']['phylo'][genome]
    run:
        R("""
        library(dplyr)
        library(VariantFiltering)
        uind_param <- VariantFilteringParam(vcfFilenames="{input.vcf}",
                                       bsgenome="{params.bsgenome}", 
                                       txdb="{params.txdb}",   
                                       snpdb="{params.snpdb}",
                                       otherAnnotations=c("{params.esp}",
                                                          "{params.exac}",
                                                          "{params.sift}",
                                                          "{params.phylo}"
                                                         )
        )
        cat("loading unrelated\n")
        uind <- unrelatedIndividuals(uind_param)
        save(uind,file="{output.uind}")
        """)

rule variantAnalysisSetupModel:
    input:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{family}_{pro,\w+}.{ext}.vcf.bgz",
        ped = config['types'][genome] + config['results']['analysis'] + "/{family}_{pro,\w+}.pedfile"
    output:
        result = config['types'][genome] + config['results']['analysis'] + "/{family}_{pro,\w+}.{ext}.{model,(denovo|arhomo|cmpdhet|xlinked)}.RData"
    params:
        bsgenome = config['analysis']['bsgenome'][genome],
        txdb     = config['analysis']['txdb'][genome],
        snpdb    = config['analysis']['snpdb'][genome],
        esp      = config['analysis']['esp'][genome],
        exac     = config['analysis']['exac'][genome],
        sift     = config['analysis']['sift'][genome],
        phylo    = config['analysis']['phylo'][genome]
    run:
        model_lut = {"denovo":"deNovo","arhomo":"autosomalRecessiveHomozygous","cmpdhet":"autosomalRecessiveHeterozygous","xlinked":"xLinked"}
        modelname=model_lut[wildcards.model]
        if wildcards.model == 'xlinked' and getGender(wildcards.pro) == 'F':
            print(wildcards.pro + ' is female. Why are you asking for x-linked variants?')
            quit()
        R("""
        library(dplyr)
        library(VariantFiltering)
        var_param <- VariantFilteringParam(vcfFilenames="{input.vcf}",
                                       pedFilename="{input.ped}", bsgenome="{params.bsgenome}", 
                                       txdb="{params.txdb}",   
                                       snpdb="{params.snpdb}",
                                       otherAnnotations=c("{params.esp}",
                                                          "{params.exac}",
                                                          "{params.sift}",
                                                          "{params.phylo}"
                                                         )
        )
        cat("loading {modelname}\n")
        varresult<-{modelname}(var_param)
        save(varresult,file="{output.result}")
        """)

def getGender(proband):
    res = sample_table.loc[sample_table['Subject'] ==  proband,"Sex"]
    sex = str(res.values[0])
    return sex

# boys only
def xlinked(wildcards):
    if getGender(wildcards.pro)=='M':
        return config['types'][genome] + config['results']['analysis'] + "/{0}_{1}.{2}.xlinked.RData".format(wildcards.family,wildcards.pro,wildcards.ext)
    else:
        return []

rule variantAnalysisModels:
    input:
        denovo = config['types'][genome] + config['results']['analysis'] + "/{family}_{pro,\w+}.{ext}.denovo.RData",
        arhomo = config['types'][genome] + config['results']['analysis'] + "/{family}_{pro,\w+}.{ext}.arhomo.RData",
        cmpdhet = config['types'][genome] + config['results']['analysis'] + "/{family}_{pro,\w+}.{ext}.cmpdhet.RData",
        xlinked = xlinked,
        ped = config['types'][genome] + config['results']['analysis'] + "/{family}_{pro,\w+}.pedfile",
        source = "reports/grin_epilepsy_models.Rmd"
    output:
        html = config['types'][genome] + config['results']['analysis'] + "/{family}_{pro,\w+}.{ext}.models.html"
    params:
        dirpath = config['types'][genome] + config['results']['analysis'],
        outfile = "/{family}_{pro,\w+}.{ext}.denovo.html"
    run:
        R("""
        library(rmarkdown)
        denovo<-get(load('{input.denovo}'))
        arhomo<-get(load('{input.arhomo}'))
        cmpdhet<-get(load('{input.cmpdhet}'))
        if('{input.xlinked}' != ''){{
            xlinked<-get(load('{input.xlinked}'))
        }}
        mytrio<-"{wildcards.family} ({wildcards.pro})"
        ped <-read.table("{input.ped}",header=TRUE)
        rmarkdown::render("{input.source}",output_file="{params.outfile}",output_dir="{params.dirpath}")
        """)

        
    
rule run_denovogear:
    input:
        vcf = config['types'][genome] + config['results']['vcfs'] + "/{file}.trio.vcf",
        ped = config['pedfile']
    output:
        dnm_auto = config['types'][genome] + config['results']['analysis'] + "/{file}.dnm_auto.txt"
    log: 
        config['datadirs']['log'] + "/{file}.dnm_auto.log.log" 
    shell:
        """
        echo "SNP_INDEL CHILD_ID chr pos ref alt maxlike_null pp_null tgt_null(child/mom/dad) snpcode code maxlike_dnm pp_dnm tgt_dnm lookup flag child_rd dad_rd mom_rd child_mq dad_mq mom_mq" > {output.dnm_auto}
        dng dnm auto --vcf {input.vcf} --ped {input.ped} | \
        sed 's/CHILD_ID: //g;s/chr: //g;s/pos: //g;s/ref: //g;s/alt: //g;s/maxlike_null: //g;s/pp_null: //g;s/tgt: //g;s/tgt_null(child\/mom\/dad): //g;s/snpcode: //g;s/code: //g;s/maxlike_dnm: //g;s/pp_dnm: //g;s/tgt_dnm(child\/mom\/dad): //g;s/lookup: //g;s/flag: //g;s/READ_DEPTH child: //g;s/dad: //g;s/mom: //g;s/MAPPING_QUALITY child: //g;s/dad: //g;s/mom: //g' >> {output.dnm_auto}
        """

#### Report ####
# create YAML file used in meta-FastQC report
rule makeyaml:
    output: 
        yaml = "fastqc.yaml",
    params:
        projdir = config['projdir'],
        fastqc = config['datadirs']['fastqc'],
        samplelanes = SAMPLELANES
    run:
        with open(output.yaml, "w") as out:
           idx = 1
           out.write("paired: yes\n") 
           out.write("output: {0}\n".format(params.projdir)) 
           out.write("fastqc:\n") 
           for name in params.samplelanes:
               out.write("  {0}\n".format(name)) 
               out.write("  - {0}/{1}_R1_fastqc.zip\n".format(params.fastqc,name))
               out.write("  - {0}/{1}_R2_fastqc.zip\n".format(params.fastqc,name))

#### Create fastqc summary

rule fastqc_summary:
    input: yaml = 'summary_fastqc.yaml'
    output: html = 'summary_fastqc.html'
    params: projdir = config['projdir']
    run: 
        R("""
        PROJECT_HOME<-"{params.projdir}";
        path.out<-"{params.projdir}/fastqc/summary";
        fn.yaml<-"{params.projdir}/summary_fastqc.yaml";
        knitr::knit("reports/summary_fastqc.Rmd")
        rmarkdown::render('summary_fastqc.md', output_format='html_document')
        """)


rule make_yaml:
    input: FASTQCS
    output: yaml = 'summary_fastqc.yaml'
    run:
        # print(FASTQCS)
        NAMES = [re.sub("\_R1_fastqc\.zip$", "", os.path.basename(name)) for name in FASTQCS if re.search("_R1_fastqc\.zip$", name)]
        # print(NAMES)
        with open(output.yaml, "w") as out:
            out.write("paired: yes\n")
            out.write("output: " + config['summary_fastqc'] + "\n");
            out.write("fastqc:\n");
            for name in NAMES:
                out.write('  ' + name + ":\n");
                out.write('  - ' + config['datadirs']['fastqc'] + '/' + name + "_R1_fastqc.zip\n");
                out.write('  - ' + config['datadirs']['fastqc'] + '/' + name + "_R2_fastqc.zip\n");


#### Create Markdown index of FastQC report files
rule makemd:
    output: 
        md = config['datadirs']['website'] + "/fastqc.md",
    run:
        note = """
The latest version (v0.11.4) of FastQC was downloaded from the following URL:
[http://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc](http://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc)
and installed as &lt;isilon&gt;/bin/fastqc.

"""
        with open(output.md, "w") as out:
           idx = 1
           out.write("## FastQC Reports\n") 
           out.write(note)
           for name in BASENAMES:
               out.write(" " + str(idx) + ".") 
               out.write(" **" + name + "**") 
               out.write(" [R1]({{SLINK}}/fastqc/{0}_R1_fastqc.html)".format(name))
               out.write(" [R2]({{SLINK}}/fastqc/{0}_R2_fastqc.html)".format(name))
               out.write("\n\n")
               idx += 1

rule siteindex:
    #input: ANALYSES,COMPLETETRIOSFAMIDS,ANALYSISREADY
    output:
        config['datadirs']['website'] + "/index.md"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("""
#### De novo analysis reports
""")
            for s,p in zip(ANALYSES,COMPLETETRIOSFAMIDS):
                outfile.write("> [`{0}`]({1}/{2})\n\n".format(p, SLINK, s))
            outfile.write("""
#### Annotated trio vcfs
""")
            for s in ANALYSISREADY:
                outfile.write("> [`{0}`]({1}/{2})\n\n".format(s, SLINK, s))

            # this link won't be treated as a child page
            #outfile.write("[fastqc summary]({{SLINK}}/summary_fastqc.html)\n")
            # must add a child index to the mybic project admin page to treat this link as a child page
            outfile.write('<a href="summary_fastqc.html">fastqc summary</a>' + "\n")
            outfile.write("<p>\n")
            #outfile.write("[multiqc report]({{SLINK}}/" + config['datadirs']['multiqc'] + "/multiqc_report.html)\n")
            outfile.write('<a href="multiqc_report.html">multiqc report</a>' + "\n")


#### Internal
onsuccess:
    print("Workflow finished, no error")
    shell("mail -s 'workflow finished' "+config['admins']+" < {log}")

onerror:
    print("An error occurred")
    shell("mail -s 'an error occurred' "+config['admins']+" < {log}")


