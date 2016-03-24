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
source activate snakeenv
snakemake -j -c "qsub -l h_vmem=40G -l mem_free=40G" 
"""

configfile: "configs/baseconfig.yaml"
configfile: "configs/config.yaml"

ENV3 = '{condaenv}/'.format(condaenv=config['python3_environment'])
ENV2 = '{condaenv}/'.format(condaenv=config['python2_environment'])

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
UNMANIFESTEDPAIRS = [name for name in ALLPAIRNAMES if name not in PAIRNAMESINSAMPLETABLE]

# pair up
# ['E0974_GCTACGC_L006', 'E0975_CGAGGCT_L006', 'E0977_GTAGAGG_L006']
SAMPLELANES = set([name.rsplit("_",maxsplit=1)[0] for name in PAIRNAMESINSAMPLETABLE])

EXISTINGSAMPLES = set([name.split("_",maxsplit=1)[0] for name in SAMPLELANES])

# a quad produces two trios
COMPLETETRIOSFAMIDS = sorted(list(set([row['FamilyID']+'_'+row['Subject'] for index, row in sample_table.iterrows() if all([row[member] in EXISTINGSAMPLES for member in ['Mother','Father','Subject']])])))
TRIOVCFS = [config['datadirs']['vcfs'] + "/" + trio + ".trio.phased.vcf" for trio in COMPLETETRIOSFAMIDS]

# quads are one family
COMPLETEFAMILYFAMIDS = set([row['FamilyID'] for index, row in sample_table.iterrows() if all([row[member] in EXISTINGSAMPLES for member in ['Mother','Father','Subject']])])
FAMILYVCFS = [config['datadirs']['vcfs'] + "/" + trio + ".family.vcf" for trio in COMPLETEFAMILYFAMIDS]

INCOMPLETEFAMILIES = set([row['FamilyID'] for index, row in sample_table.iterrows() if any([row[member] not in EXISTINGSAMPLES and not pandas.isnull(row[member]) for member in ['Mother','Father','Subject']])])
TRIOGEMS = [config['datadirs']['gemini'] + "/" + trio + ".gemini.db" for trio in COMPLETEFAMILYFAMIDS]

ANALYSISREADY = [config['datadirs']['vcfs'] + "/" + trio + ".trio.phased.com.filtered.ad.de.nm.snpeff.vcf.bgz" for trio in COMPLETETRIOSFAMIDS]
RDATA         = [config['datadirs']['analysis'] + "/" + trio + ".trio.phased.com.filtered.ad.de.nm.snpeff.denovo.RData" for trio in COMPLETETRIOSFAMIDS]
ANALYSES      = [config['datadirs']['analysis'] + "/" + trio + ".trio.phased.com.filtered.ad.de.nm.snpeff.denovo.html" for trio in COMPLETETRIOSFAMIDS]

SAMS = [config['datadirs']['sams'] + "/" + name + ".sam" for name in SAMPLELANES]
BAMS = [config['datadirs']['bams'] + "/" + name + ".bam" for name in SAMPLELANES]
SBAMS = [config['datadirs']['bams'] + "/" + name + ".sorted.bam" for name in SAMPLELANES]
MBAMS = [config['datadirs']['bams'] + "/" + name + ".sorted.merged.bam" for name in EXISTINGSAMPLES]
DBAIS = [config['datadirs']['picard'] + "/" + name + ".rmdup.bai" for name in EXISTINGSAMPLES]
DBAMS = [config['datadirs']['picard'] + "/" + name + ".rmdup.bam" for name in EXISTINGSAMPLES]
GBAIS = [config['datadirs']['picard'] + "/" + name + ".group.bai" for name in EXISTINGSAMPLES]
GBAMS = [config['datadirs']['picard'] + "/" + name + ".group.bam" for name in EXISTINGSAMPLES]
RBAMS = [config['datadirs']['realigned'] + "/" + name + ".bam" for name in EXISTINGSAMPLES]
LISTS = [config['datadirs']['lists'] + "/" + name + ".list" for name in EXISTINGSAMPLES]
TABLES = [config['datadirs']['recalibrated'] + "/" + name + ".table" for name in EXISTINGSAMPLES]
RECBAMS = [config['datadirs']['recalibrated'] + "/" + name + ".bam" for name in EXISTINGSAMPLES]
POSTTABLES = [config['datadirs']['postrecalibrated'] + "/" + name + ".table" for name in EXISTINGSAMPLES]
PDFS = [config['datadirs']['pdfs'] + "/" + name + ".pdf" for name in EXISTINGSAMPLES]
GVCFS = [config['datadirs']['gvcfs'] + "/" + name + ".gvcf" for name in EXISTINGSAMPLES]
GVCFSLIST = ' '.join(["--variant " + config['datadirs']['gvcfs'] + "/" + name + ".gvcf" for name in EXISTINGSAMPLES])

ANNOVARDBS = [config['annovardbdir'] + "/" + config['buildve'] + "_" + db + ".installed" for db in config['annovardbs']]

ANNOVAR_PROTOCOLS = ','.join(config['annovardbs'])


INDELS = config['datadirs']['realigned'] + "/indels.list"

dlocs = dict()

workdir: config['projdir']

rule all:
    input: 
        trios = TRIOVCFS,
        analysis = ANALYSES,
        phased = config['datadirs']['vcfs'] + "/joint.family.vcf" # must run after all gvcf files created; will create joint.vcf if not already
#include TRIOGEMS for gemini (GRCh37 only)

#this is useful for extracting sequences from bams
rule extract:
    input: EXTRACT

rule rdata:
    input: RDATA

rule triovcfs:
    input: TRIOVCFS

rule analysisready:
    input: ANALYSISREADY

rule analyses:
    input: ANALYSES

# this is a utility to put things in the correct order in case something upstream gets touched
rule catchup:
    params:
        picard = config['datadirs']['picard'],
        lists = config['datadirs']['lists'],
        realigned = config['datadirs']['realigned'],
        recalibrated = config['datadirs']['recalibrated'],
        postrecalibrated = config['datadirs']['postrecalibrated'],
        gvcfs = config['datadirs']['gvcfs'],
        vcfs = config['datadirs']['vcfs'],
        analysis = config['datadirs']['analysis']
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
        ms="missingsamples.txt", ump="unmanifestedpairs.txt", ic="incompletefamilies.txt"
    run:
        assert(len(SAMPLELANES) == len(PAIRNAMESINSAMPLETABLE)/2)
        print("Received Pairs (on disk): {0}".format(len(ALLPAIRNAMES)))
        print("Unmanifested Pairs (on disk, not in sample table): {0}".format(len(UNMANIFESTEDPAIRS)))
        f = open(output.ump, 'w')
        for pair in UNMANIFESTEDPAIRS:
            f.write("{0}\n".format(pair))
        print("Existing Samples (on disk, in sample table): {0}".format(len(EXISTINGSAMPLES)))
        print("Missing Samples (in sample table, not on disk): {0}".format(len(MISSINGSAMPLES)))
        f = open(output.ms, 'w')
        for sample in MISSINGSAMPLES:
            f.write("{0}\n".format(sample))
        print("Manifested Pairs (in sample table): {0}".format(len(PAIRNAMESINSAMPLETABLE)))
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
            shell("mkdir -p " + config['datadirs'][adir])

rule dummy:    # just to test the python codes above
    input:  workflow.basedir + "/Snakefile"

    run:
        check_gvcfs(GVCFS)

rule target_lists:
    input: LISTS

rule print_reads:
    input: RECBAMS

rule join_gvcfs:
    input: config['datadirs']['vcfs'] + "/joint.vcf"

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
        """
        fastq=os.path.basename(input).replace('_001','')
        os.symlink(input,fastq)
        """

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
        metrics = config['datadirs']['picard']
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

rule mirror:
    input:
        expand("aws/{filename}_{pair}.fastq.gz.sent",filename=SAMPLELANES,pair=['R1','R2'])
#        S3.remote(expand("cbttc.seq.data/Ingo_project/{filename}.fastq.gz",filename=SAMPLELANES))

#### Alignment ####
rule align:
    input:
        pair1 = config['datadirs']['fastq'] + "/{sample}_R1.fastq.gz",
        pair2 = config['datadirs']['fastq'] + "/{sample}_R2.fastq.gz",
        align = ENV3 + config['tools']['align']
    output:
        sam = config['datadirs']['sams'] + "/{sample}.sam" # may be set to temp
    threads:
        12
    log: 
        config['datadirs']['log'] + "/{sample}.novoalign.log"
    params:
        refidx = config['refidx']
    shell:
        """
        {input.align} -c {threads} -a -k -d {params.refidx} -o SAM -f {input.pair1} {input.pair2} 1> {output.sam} 2> {log}
        """

rule sam_to_bam:
    input:
        sam = config['datadirs']['sams'] + "/{sample}.sam",
        samtools = ENV3 + config['tools']['samtools']
    output:
        bam = config['datadirs']['bams'] + "/{sample,[^.]+}.bam"
    threads:
        12   # also depends on -j
    shell:
        """
        {input.samtools} view -@ {threads} -bS {input.sam} > {output.bam}
        """

# novosort creates index
rule novosortbam:
    input:
        bam = config['datadirs']['bams'] + "/{sample}.bam",
        sort = ENV3 + config['tools']['sortbam']
    output:
        sorted = config['datadirs']['bams'] + "/{sample}.sorted.bam",
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
        bai = config['datadirs']['picard'] + "/{sample}.group.bai",
        bam = config['datadirs']['picard'] + "/{sample}.group.bam",
        java = ENV3 + config['tools']['java']
    output:
        samplelist = config['datadirs']['lists'] + "/{sample}.list"
    log:
        config['datadirs']['log'] + "/{sample}.target_list.log"
    params:
        jar = config['jars']['gatk'],
        opts = config['tools']['opts']['high'],
        ref = config['ref'],
        knownsites = config['known']
    threads:
        24
    shell:
        """
        {input.java} {params.opts} -jar {params.jar} \
        -T RealignerTargetCreator \
        -nt {threads} \
        -R {params.ref} \
        -I {input.bam} \
        -known {params.knownsites} \
        -o {output.samplelist} 2> {log}
        """

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
        list = config['datadirs']['lists'] + "/{sample}.list",
        dbam = config['datadirs']['picard'] + "/{sample}.group.bai",
        java = ENV3 + config['tools']['java']
    output:
        rbam = config['datadirs']['realigned'] + "/{sample}.bam"
    params:
        jar = config['jars']['gatk'],
        opts = config['tools']['opts']['med'],
        ref = config['ref'],
        known = config['known']
    shell:
        """
        {input.java} {params.opts} -jar {params.jar} \
        -T IndelRealigner \
        -R {params.ref} \
        -I {config[datadirs][picard]}/{wildcards.sample}.group.bam \
        -targetIntervals {input.list} \
        -known {params.known} \
        -o {output.rbam}
        """

# Base recalibration (not be confused with variant recalibration)
# http://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr
# https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php
rule generate_recalibration_table:
    input:
        bam = config['datadirs']['realigned'] + "/{sample}.bam",
        java = ENV3 + config['tools']['java']
    output:
        table = config['datadirs']['recalibrated'] + "/{sample}.table"
    log:
        config['datadirs']['log'] + "/{sample}.generate_recalibration_table.log"
    params:
        jar = config['jars']['gatk'],
        opts = config['tools']['opts']['med'],
        ref = config['ref'],
        known = config['known']
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
        table = config['datadirs']['recalibrated'] + "/{sample}.table",
        bam = config['datadirs']['realigned'] + "/{sample}.bam",
        java = ENV3 + config['tools']['java']
    output:
        bam = config['datadirs']['recalibrated'] + "/{sample}.bam"
    log:
        config['datadirs']['log'] + "/{sample}.recalibrate_bam.log"
    params:
        jar = config['jars']['gatk'],
        opts = config['tools']['opts']['med'],
        ref = config['ref']
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
        table = config['datadirs']['recalibrated'] + "/{sample}.table",
        bam = config['datadirs']['realigned'] + "/{sample}.bam",
        java = ENV3 + config['tools']['java']
    output:
        table = config['datadirs']['postrecalibrated'] + "/{sample}.table",
    params:
        jar = config['jars']['gatk'],
        opts = config['tools']['opts']['med'],
        ref = config['ref'],
        known = config['known']
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
        before = config['datadirs']['recalibrated'] + "/{sample}.table",
        after = config['datadirs']['postrecalibrated'] + "/{sample}.table",
        java = ENV3 + config['tools']['java']
    output:
        pdf = config['datadirs']['pdfs'] + "/{sample}.pdf"
    params:
        jar = config['jars']['gatk'],
        opts = config['tools']['opts']['med'],
        ref = config['ref']
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
        bam = config['datadirs']['recalibrated'] + "/{sample}.bam",
        java = ENV3 + config['tools']['java']
    output:
        "{sample}.DoC"
    params:
        jar = config['jars']['gatk'],
        opts = config['tools']['opts']['med'],
        ref = config['ref']
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
        bam = config['datadirs']['bams'] + "/{sample}.sorted.merged.bam",
        java = ENV3 + config['tools']['java']
    output:
        bam = config['datadirs']['picard'] + "/{sample}.rmdup.bam"
    log:
        config['datadirs']['log'] + "/{sample}.markdups.log"
    params:
        picard = config['jars']['picard']['path'],
        md = config['jars']['picard']['markdups'],
        opts = config['tools']['opts']['med'],
        metrics = config['datadirs']['picard']
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
        bam = config['datadirs']['picard'] + "/{sampleandext}.bam",
        samtools = ENV3 + config['tools']['samtools']
    output:
        bai = config['datadirs']['picard'] + "/{sampleandext}.bai"
    shell:
        """
        {input.samtools} index {input.bam} {output.bai}
        """


rule add_readgroup:
    input:
        bam = config['datadirs']['picard'] + "/{sample}.rmdup.bam",
        bai = config['datadirs']['picard'] + "/{sample}.rmdup.bai",
        java = ENV3 + config['tools']['java']
    output:
        bam = config['datadirs']['picard'] + "/{sample}.group.bam"
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
        bam = config['datadirs']['recalibrated'] + "/{sample}.bam",
        java = ENV3 + config['tools']['java']
    output:
        gvcf = config['datadirs']['gvcfs'] + "/{sample}.gvcf"
    params:
        jar = config['jars']['gatk'],
        opts = config['tools']['opts']['high'],
        ref = config['ref']
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
        gvcfs = [config['datadirs']['gvcfs'] + "/" + name + ".gvcf" for name in set(samples)]
        gvcfslist = ' '.join(["--variant " + config['datadirs']['gvcfs'] + "/" + name + ".gvcf" for name in set(samples)])
        return [gvcfs, gvcfslist]

# make sure the family is done first
# if family size is 3, then you're done, just copy
# otherwise make both trios
rule trio_vcfs:
    input:
        gvcfs = lambda wildcards: gvcf_samples_in_family(wildcards.family,wildcards.subject)[0],
        family = config['datadirs']['vcfs'] + "/{family}.family.vcf",
        familygvcfs = lambda wildcards: gvcf_samples_in_family(wildcards.family,'family')[0],
        java = ENV3 + config['tools']['java']
    output:
        vcf = config['datadirs']['vcfs'] + "/{family}_{subject}.trio.vcf"
    log:
        config['datadirs']['log'] + "/{family}_{subject}.trio.vcf.log"
    params:
        jar = config['jars']['gatk'],
        ref = config['ref'],
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
        vcf = config['datadirs']['vcfs'] + "/{family}.family.vcf",
        idx = config['datadirs']['vcfs'] + "/{family}.family.vcf.idx"
    log:
        config['datadirs']['log'] + "/{family}.family.vcf.log"
    params:
        jar = config['jars']['gatk'],
        ref = config['ref'],
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
        config['datadirs']['analysis'] + "/{family}_{subject}.pedfile"
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
        vcf = config['datadirs']['vcfs'] + "/{file}.trio.vcf",
        ped = config['pedfile'],
        java = ENV3 + config['tools']['java']
    output:
        vcf = config['datadirs']['vcfs'] + "/{file}.trio.phased.vcf",
        idx = config['datadirs']['vcfs'] + "/{file}.trio.phased.vcf.idx",
        mvf = config['datadirs']['vcfs'] + "/{file}.mendelian_violations.txt"
    params:
        jar  = config['jars']['gatk'],
        ref = config['ref'],
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
        vcf = config['datadirs']['vcfs'] + "/{file}.vcf",
        java = ENV3 + config['tools']['java']
    output:
        vcf = config['datadirs']['vcfs'] + "/{file}.snps.vcf",
        idx = config['datadirs']['vcfs'] + "/{file}.snps.vcf.idx"
    params:
        jar  = config['jars']['gatk'],
        ref = config['ref'],
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
        vcf = config['datadirs']['vcfs'] + "/{file}.vcf",
        java = ENV3 + config['tools']['java']
    output:
        vcf = config['datadirs']['vcfs'] + "/{file}.indels.vcf",
        idx = config['datadirs']['vcfs'] + "/{file}.indels.vcf.idx"
    params:
        jar  = config['jars']['gatk'],
        ref = config['ref'],
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
        vcf = config['datadirs']['vcfs'] + "/{file}.snps.vcf",
        java = ENV3 + config['tools']['java']
    output:
        vcf = config['datadirs']['vcfs'] + "/{file}.snps.hard.vcf",
        idx = config['datadirs']['vcfs'] + "/{file}.snps.hard.vcf.idx"
    params:
        jar  = config['jars']['gatk'],
        ref = config['ref'],
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
        vcf = config['datadirs']['vcfs'] + "/{file}.indels.vcf",
        java = ENV3 + config['tools']['java']
    output:
        vcf = config['datadirs']['vcfs'] + "/{file}.indels.hard.vcf",
        idx = config['datadirs']['vcfs'] + "/{file}.indels.hard.vcf.idx"
    params:
        jar  = config['jars']['gatk'],
        ref = config['ref'],
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
        vcf = config['datadirs']['vcfs'] + "/{file}.{type,(snps|indels)}.hard.vcf",
        java = ENV3 + config['tools']['java']
    output:
        vcf = config['datadirs']['vcfs'] + "/{file}.{type,(snps|indels)}.filtered.vcf",
        idx = config['datadirs']['vcfs'] + "/{file}.{type,(snps|indels)}.filtered.vcf"
    params:
        jar  = config['jars']['gatk'],
        ref = config['ref'],
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
        snps = config['datadirs']['vcfs'] + "/{file}.snps.filtered.vcf",
        indels = config['datadirs']['vcfs'] + "/{file}.indels.filtered.vcf",
        java = ENV3 + config['tools']['java']
    output:
        vcf = config['datadirs']['vcfs'] + "/{file}.com.filtered.vcf",
        idx = config['datadirs']['vcfs'] + "/{file}.com.filtered.vcf.idx"
    params:
        jar  = config['jars']['gatk'],
        ref = config['ref'],
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
        vcf = config['datadirs']['vcfs'] + "/{file}.vcf",
    output:
        vcf = config['datadirs']['vcfs'] + "/{file}.ad.vcf"
    params:
        ref = config['ref']
    shell:
        """
        cat {input} | sed 's/ID=AD,Number=./ID=AD,Number=R/' > {output}
        """

# decomposes multiallelic variants into biallelic in a VCF file.
rule decompose_for_gemini:
    input:
        vcf = config['datadirs']['vcfs'] + "/{file}.ad.vcf",
        vt = ENV3 + config['tools']['vt']
    output:
        vcf = config['datadirs']['vcfs'] + "/{file}.ad.de.vcf"
    shell:
        """
        {input.vt} decompose -s -o {output} {input.vcf}
        """

rule normalize_for_gemini:
    input:
        vcf = config['datadirs']['vcfs'] + "/{file}.ad.de.vcf",
        vt = ENV3 + config['tools']['vt']
    output:
        vcf = config['datadirs']['vcfs'] + "/{file}.ad.de.nm.vcf"
    params:
        ref = config['ref']
    shell:
        """
        {input.vt} normalize -r {params.ref} -o {output} {input.vcf}
        """

rule vcf_qt:
    input:
        vcf = config['datadirs']['vcfs'] + "/{file}.ad.de.vcf",
        vt = ENV3 + config['tools']['vt']
    output:
        vcf = config['datadirs']['vtpeek'] + "/{file}.vtpeek.txt"
    params:
        ref = config['ref']
    shell:
        """
        {input.vt} peek -o {output} {input.vcf}
        """

rule vcf_profile:
    input:
        vcf = config['datadirs']['vcfs'] + "/{file}.ad.de.vcf",
        vt = ENV3 + config['tools']['vt'],
        ped = config['pedfile']
    output:
        vcf = config['datadirs']['vtpeek'] + "/{file}.vtmendelprofile.txt"
    params:
        ref = config['ref']
    shell:
        """
        {input.vt} profile_mendelian -o {output} -p {input.ped} -x mendel {input.vcf}
        """
        
# ud - upstream downstream interval length (in bases)
rule run_snpeff:
    input:
        vcf = config['datadirs']['vcfs'] + "/{file}.vcf",
        java = ENV3 + config['tools']['java']
    output:
        vcf = config['datadirs']['vcfs'] + "/{file}.snpeff.vcf"
    params:
        jar  = config['jars']['snpeff']['path'],
        conf = config['jars']['snpeff']['cnf'],
        opts = config['tools']['opts']['med'],
        database = config['jars']['snpeff']['db'],
        updown = config['jars']['snpeff']['ud'],
        format = config['jars']['snpeff']['format']
    shell:
        """
        {input.java} {params.opts} -jar {params.jar} \
        -c {params.conf} \
        -t {params.database} \
        -ud {params.updown} \
        {params.format} \
         {input.vcf} > {output.vcf}
        """

#### run VEP  ####

rule run_vep:
    input: "{file}.vcf.gz"
    output: "{file}.vep.vcf.gz"
    run:
        shell("""
          perl ./vep/ensembl-tools-release-78/scripts/variant_effect_predictor/variant_effect_predictor.pl \
              --everything --vcf --allele_number --no_stats --cache --offline \
              --dir ./vep_cache/ --force_overwrite --cache_version 78 \
              --fasta ./vep_cache/homo_sapiens/78_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
              --assembly GRCh37 --tabix \
              --plugin LoF,human_ancestor_fa:./loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15 \
              --plugin dbNSFP,./reference_data/dbNSFP/dbNSFPv2.9.gz,Polyphen2_HVAR_pred,CADD_phred,SIFT_pred,FATHMM_pred,MutationTaster_pred,MetaSVM_pred \
              -i {wildcards.file}.vcf.gz -o {wildcards.file}.vep.vcf.gz
        """)

#### run multiqc  ####

rule run_multiqc:
    input:
        gbams = GBAMS,
        multiqc = ENV3 + config['tools']['multiqc']
    output: config['datadirs']['multiqc'] + '/multiqc_report.html'
    params:
        dirs = config['datadirs']['picard'] + ' fastqc',
        outdir = config['datadirs']['multiqc'] 
    shell:
        """
        {input.multiqc} -o {params.outdir} {params.dirs} # will detect input file types?
        """

#### run annovar  ####

rule table_annovar:
    input:
        ANNOVARDBS,
        avinput = config['datadirs']['vcfs'] + "/{file}.vcf",
        annovar = ENV3 + config['tools']['table_annovar']
    output:
        config['datadirs']['vcfs'] + "/{file}.annovar.vcf"
    params:
        opts = "-buildver "+config['buildve']
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
        avinput = config['datadirs']['vcfs'] + "/{file}.avinput",
        annovar = ENV3 + config['tools']['annotate_variation']
    output:
        config['datadirs']['vcfs'] + "/annovar.done"
    params:
        opts = "-buildver {config['buildve']}",
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
        vcf = config['datadirs']['vcfs'] + "/{file}.vcf",
        cmd = ENV3 + config['tools']['vcf2avinput']
    output:
        config['datadirs']['vcfs'] + "/{file}.avinput",
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
        vcf = config['datadirs']['vcfs'] + "/{file}.vcf",
        bgzip = ENV3 + config['tools']['bgzip']
    output:
        vcf = config['datadirs']['vcfs'] + "/{file}.vcf.bgz",
    shell:
        """
        {input.bgzip} -c {input.vcf} > {output}
        """

rule tabix:
    input:
        vcf = config['datadirs']['vcfs'] + "/{file}.vcf.bgz",
        tabix = ENV3 + config['tools']['tabix']
    output:
        vcf = config['datadirs']['vcfs'] + "/{file}.vcf.bgz.tbi",
    shell:
        """
        {input.tabix} -p vcf {input.vcf}
        """

rule gemini_db:
    input:
        vcf = config['datadirs']['vcfs'] + "/{file}.trio.phased.com.filtered.ad.de.nm.snpeff.vcf.bgz",
        tbi = config['datadirs']['vcfs'] + "/{file}.trio.phased.com.filtered.ad.de.nm.snpeff.vcf.bgz.tbi",
        ped = config['pedfile'],
        gemini = ENV2 + config['tools']['gemini']
    output:
        config['datadirs']['gemini'] + "/{file}.gemini.db"
    threads:
        3
    shell:
        """
        {ENV2}
        {input.gemini} load --cores {threads} -t snpEff -v {input.vcf} -p {input.ped} {output}
        """

rule testR:
    run:
        R("""
        rnorm(100)
        """)

rule describeR:
    run:
        R("""
        library(dplyr)
        library(VariantFiltering)
        cat(.libPaths())
        print(sessionInfo())
        """)

#trio.phased.com.filtered.ad.de.nm.snpeff.vcf.bgz
#### Analysis ####
rule variantAnalysisSetupUind:
    input:
        vcf = config['datadirs']['vcfs'] + "/{familypro}.vcf.bgz",
    output:
        uind = config['datadirs']['analysis'] + "/{familypro}.uind.RData",
    params:
        rlibrary = config['analysis']['rlibrary'],
        bsgenome = config['analysis']['bsgenome'],
        txdb     = config['analysis']['txdb'],
        snpdb    = config['analysis']['snpdb'],
        esp      = config['analysis']['esp'],
        exac     = config['analysis']['exac'],
        sift     = config['analysis']['sift'],
        phylo    = config['analysis']['phylo']
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

rule variantAnalysisSetupDeNovo:
    input:
        vcf = config['datadirs']['vcfs'] + "/{family}_{pro,\w+}.{ext}.vcf.bgz",
        ped = config['datadirs']['analysis'] + "/{family}_{pro,\w+}.pedfile"
    output:
        denovo = config['datadirs']['analysis'] + "/{family}_{pro,\w+}.{ext}.denovo.RData"
    params:
        rlibrary = config['analysis']['rlibrary'],
        bsgenome = config['analysis']['bsgenome'],
        txdb     = config['analysis']['txdb'],
        snpdb    = config['analysis']['snpdb'],
        esp      = config['analysis']['esp'],
        exac     = config['analysis']['exac'],
        sift     = config['analysis']['sift'],
        phylo    = config['analysis']['phylo']
    run:
        R("""
        library(dplyr)
        library(VariantFiltering)
        denovo_param <- VariantFilteringParam(vcfFilenames="{input.vcf}",
                                       pedFilename="{input.ped}", bsgenome="{params.bsgenome}", 
                                       txdb="{params.txdb}",   
                                       snpdb="{params.snpdb}",
                                       otherAnnotations=c("{params.esp}",
                                                          "{params.exac}",
                                                          "{params.sift}",
                                                          "{params.phylo}"
                                                         )
        )
        cat("loading denovo\n")
        denovo<-deNovo(denovo_param)
        save(denovo,file="{output.denovo}")
        """)

rule variantAnalysisAll:
    input:
        uind = config['datadirs']['analysis'] + "/{family}_{pro,\w+}.{ext}.uind.RData",
        denovo = config['datadirs']['analysis'] + "/{family}_{pro,\w+}.{ext}.denovo.RData",
        ped = config['datadirs']['analysis'] + "/{family}_{pro,\w+}.pedfile",
        source = "reports/grin_epilepsy.Rmd"
    output:
        html = config['datadirs']['analysis'] + "/{family}_{pro,\w+}.{ext}.all.html"
    params:
        rlibrary = config['analysis']['rlibrary']
    run:
        R("""
        library(rmarkdown)
        uind<- get(load('{input.uind}'))
        mytrio<-"{wildcards.trio} ({wildcards.pro})"
        ped <-read.table("{input.ped}",header=TRUE)
        rmarkdown::render("{input.source}",output_file="{output.html}")
        """)

rule variantAnalysisDeNovo:
    input:
        denovo = config['datadirs']['analysis'] + "/{family}_{pro,\w+}.{ext}.denovo.RData",
        ped = config['datadirs']['analysis'] + "/{family}_{pro,\w+}.pedfile",
        source = "reports/grin_epilepsy_denovo.Rmd"
    output:
        html = config['datadirs']['analysis'] + "/{family}_{pro,\w+}.{ext}.denovo.html"
    params:
        rlibrary = config['analysis']['rlibrary']
    run:
        R("""
        .libPaths( c( .libPaths(), "{params.rlibrary}") )
        library(rmarkdown)
        denovo<- get(load('{input.denovo}'))
        mytrio<-"{wildcards.trio} ({wildcards.pro})"
        ped <-read.table("{input.ped}",header=TRUE)
        rmarkdown::render("{input.source}",output_file="{output.html}")
        """)
        
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


