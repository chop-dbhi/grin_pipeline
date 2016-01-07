import glob
import re
import pandas
from functools import cmp_to_key
"""
run on respublica
source activate snakeenv
snakemake -j -c "qsub -l h_vmem=40G -l mem_free=40G" 
"""

configfile: "config.yaml"

SEQLIST = "files.lst"
SAMPLELIST = "samples.lst"
SLINK = "GRIN/fastqc/"

DOWNLOADDIR = "kiel"
DOWNLOADS = glob.glob(DOWNLOADDIR + "/*/fastq/*/*/*fastq.gz")

FASTQS = glob.glob(config['datadirs']['fastq'] + "/*.gz")

#FamilyID       Subject Mother  Father  Sex     Affected_status Not_in_Varbank
#Trio_SL        C2952   C2953   C2954   f       EOEE
#ISR_#45        E08320                          f       Focal Epilepsy  x
sample_table = pandas.read_table(config['sample_table'],index_col=1)

SAMPLES = list(sample_table.index)

BASENAMES = [re.sub("\_R1\.fastq\.gz$", "", os.path.basename(name)) for name in FASTQS if re.search("_R1\.fastq\.gz$", name)]

# BASENAMES.remove("CNB01-001B_ATCACG_L008")

BASENAMES.sort()

SAMS = [config['datadirs']['sams'] + "/" + name + ".sam" for name in BASENAMES]
BAMS = [config['datadirs']['bams'] + "/" + name + ".bam" for name in BASENAMES]
DBAIS = [config['datadirs']['picard'] + "/" + name + ".rmdup.bai" for name in BASENAMES]
DBAMS = [config['datadirs']['picard'] + "/" + name + ".rmdup.bam" for name in BASENAMES]
GBAIS = [config['datadirs']['picard'] + "/" + name + ".group.bai" for name in BASENAMES]
GBAMS = [config['datadirs']['picard'] + "/" + name + ".group.bam" for name in BASENAMES]
RBAMS = [config['datadirs']['realigned'] + "/" + name + ".bam" for name in BASENAMES]
LISTS = [config['datadirs']['lists'] + "/" + name + ".list" for name in BASENAMES]
TABLES = [config['datadirs']['recalibrated'] + "/" + name + ".table" for name in BASENAMES]
RECBAMS = [config['datadirs']['recalibrated'] + "/" + name + ".bam" for name in BASENAMES]
POSTTABLES = [config['datadirs']['postrecalibrated'] + "/" + name + ".table" for name in BASENAMES]
PDFS = [config['datadirs']['pdfs'] + "/" + name + ".pdf" for name in BASENAMES]
GVCFS = [config['datadirs']['gvcfs'] + "/" + name + ".gvcf" for name in BASENAMES]
GVCFSLIST = ' '.join(["--variant " + config['datadirs']['gvcfs'] + "/" + name + ".gvcf" for name in BASENAMES])
SBAMS = [re.sub("\.bam$", ".sorted.bam", name) for name in BAMS]

INDELS = config['datadirs']['realigned'] + "/indels.list"

dlocs = dict()

workdir: config['projdir']

rule all:
    input: 
        lists = LISTS,          # becore combine lists
        indels = INDELS,        # must run after all lists are created
        rbam = RBAMS,
        phased = config['datadirs']['gvcfs'] + "/phased.vcf" # must run after all gvcf files created; will create joint.vcf if not already

rule basehead:
    run:
        print(FASTQS)

rule dummy:    # just to test the python codes above
    input: "Snakefile"

rule target_lists:
    input: LISTS

rule print_reads:
    input: RECBAMS

rule join_gvcfs:
    input: config['datadirs']['gvcfs'] + "/joint.vcf"

rule make_gvcfs:
    input: GVCFS

rule make_bams:
    input: BAMS

rule make_sams:
    input: SAMS

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

rule sortbams:
    input: DBAMS

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

### QC ####
rule fastqc: 
    input: 
        pair1 = config['datadirs']['samples'] + "/{sample}1.fastq.gz",
        pair2 = config['datadirs']['samples'] + "/{sample}2.fastq.gz",
        seq2qc = config['tools']['seq2qc']
    log: 
        "logs/{sample}.log" 
    output: 
        pair1 = config['datadirs']['fastqc'] + "/{sample}1_fastqc.zip",
        pair2 = config['datadirs']['fastqc'] + "/{sample}2_fastqc.zip",
        html1 = config['datadirs']['fastqc'] + "/{sample}1_fastqc.html",
        html2 = config['datadirs']['fastqc'] + "/{sample}2_fastqc.html"
    params:
        qcdir = config['datadirs']['fastqc']
    # how to run qsub?
    shell: "{input.seq2qc} -o {params.qcdir} {input.pair1} {input.pair2} 2> {log}" 

#### Alignment ####
rule align:
    input:
        pair1 = config['datadirs']['samples'] + "/{sample}_R1.fastq.gz",
        pair2 = config['datadirs']['samples'] + "/{sample}_R2.fastq.gz",
        makesam = config['tools']['makesam']
    output:
        sam = config['datadirs']['sams'] + "/{sample}.sam" # may be set to temp
    threads:
        12   # also depends on -j
    params:
        refseqidx = config['refseqidx']
    shell:
        """
        {input.makesam} -c {threads} -a -k -d {params.refseqidx} -o SAM -f {input.pair1} {input.pair2} > {output.sam}
        """

rule sam_to_bam:
    input:
        sam = config['datadirs']['sams'] + "/{sample}.sam",
        sam2bam = config['tools']['samtools']
    output:
        bam = config['datadirs']['bams'] + "/{sample}.bam"
    threads:
        12   # also depends on -j
    shell:
        """
        {input.sam2bam} view -@ 12 -bS {input.sam} > {output.bam}
        """

# novosort creates index
rule novosortbam:
    input:
        bam = config['datadirs']['bams'] + "/{sample}.bam",
        sort = config['tools']['sortbam']
    output:
        sorted = config['datadirs']['bams'] + "/{sample}.sorted.bam",
    shell:
        """
        {input.sort} -m 14g -t . --removeduplicates --keeptags -i -o {output.sorted} {input.bam}
        """

## here we create individual realign target lists, one from each bam file
## later, the individual lists will be combined into a single list
## an alternative ways is to directly create a signle list from the bam files

rule target_list: # create individual realign target list
    input:  # deduced bams
        bai = config['datadirs']['picard'] + "/{sample}.group.bai", # required, so make sure it's created
        bam = config['datadirs']['picard'] + "/{sample}.group.bam",
        java = config['tools']['java']
    output:
        list = config['datadirs']['lists'] + "/{sample}.list"
    params:
        jar = config['jars']['gatk'],
        javaopts = config['tools']['javaopts'],
        refseq = config['refseq']
    shell:
        """
        {input.java} {params.javaopts} -jar {params.jar} \
        -T RealignerTargetCreator \
        -R {params.refseq} \
        -I {input.bam} \
        -known {config[siv]} \
        -o {output.list}
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

rule make_group_index:
    input:
        bam = config['datadirs']['picard'] + "/{sample}.group.bam",
        java = config['tools']['java']
    output:
        bai = config['datadirs']['picard'] + "/{sample}.group.bai"
    params:
        javaopts = config['tools']['javaopts'],
        jar = config['jars']['buildbamindex']
    shell:
        """
        {input.java} {params.javaopts} -jar {params.jar} \
        INPUT={input.bam}
        """

# combine indels from individual target alignment lists into a single list
rule combine_lists:
    input: LISTS
    output: INDELS
    run:
         combine()

rule realign_target:   # with one combined list file
    input:  # deduced bams
        list = INDELS,
        dbam = config['datadirs']['picard'] + "/{sample}.group.bai",
        java = config['tools']['java']
    output:
        rbam = config['datadirs']['realigned'] + "/{sample}.bam"
    params:
        jar = config['jars']['gatk'],
        javaopts = config['tools']['javaopts'],
        refseq = config['refseq']
    shell:
        """
        {input.java} {params.javaopts} -jar {params.jar} \
        -T IndelRealigner \
        -R {params.refseq} \
        -I {config[datadirs][picard]}/{wildcards.sample}.group.bam \
        -targetIntervals {input.list} \
        -known {config[siv]} \
        -o {output.rbam}
        """

# Base recalibration (not be confused with variant recalibration)
# http://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr
# https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php
rule generate_recalibration table:
    input:
        bam = config['datadirs']['realigned'] + "/{sample}.bam",
        java = config['tools']['java']
    output:
        table = config['datadirs']['recalibrated'] + "/{sample}.table"
    params:
        jar = config['jars']['gatk'],
        javaopts = config['tools']['javaopts'],
        refseq = config['refseq']
    shell:
        """
        {input.java} {params.javaopts} -jar {params.jar} \
        -T BaseRecalibrator \
        -R {params.refseq} \
        -I {input.bam} \
        -knownSites {config[siv]} \
        -o {output.table}
        """

# https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_readutils_PrintReads.php
# Note that when PrintReads is used as part of the Base Quality Score Recalibration workflow
rule recalibrate_bam:
    input:
        table = config['datadirs']['recalibrated'] + "/{sample}.table",
        bam = config['datadirs']['realigned'] + "/{sample}.bam",
        java = config['tools']['java']
    output:
        bam = config['datadirs']['recalibrated'] + "/{sample}.bam"
    params:
        jar = config['jars']['gatk'],
        javaopts = config['tools']['javaopts'],
        refseq = config['refseq']
    shell:
        """
        {input.java} {params.javaopts} -jar {params.jar} \
        -T PrintReads \
        -R {params.refseq} \
        -I {input.bam} \
        -BQSR {input.table} \
        -o {output.bam}
        """

rule post_recalibrated_table:
    input:
        table = config['datadirs']['recalibrated'] + "/{sample}.table",
        bam = config['datadirs']['realigned'] + "/{sample}.bam",
        java = config['tools']['java']
    output:
        table = config['datadirs']['postrecalibrated'] + "/{sample}.table",
    params:
        jar = config['jars']['gatk'],
        javaopts = config['tools']['javaopts'],
        refseq = config['refseq']
    shell:
        """
        {input.java} {params.javaopts} -jar {params.jar} \
        -T BaseRecalibrator \
        -R {params.refseq} \
        -I {input.bam} \
        -knownSites {config[siv]} \
        -BQSR {input.table} \
        -o {output.table}
        """

# https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_bqsr_AnalyzeCovariates.php
# generates plots for visualizing the quality of a recalibration run
rule analyze_bqsr:
    input:
        before = config['datadirs']['recalibrated'] + "/{sample}.table",
        after = config['datadirs']['postrecalibrated'] + "/{sample}.table",
        java = config['tools']['java']
    output:
        pdf = config['datadirs']['pdfs'] + "/{sample}.pdf"
    params:
        jar = config['jars']['gatk'],
        javaopts = config['tools']['javaopts'],
        refseq = config['refseq']
    shell:
        """
        source ~/.bashrc
        module load R/3.2.2
        {input.java} {params.javaopts} -jar {params.jar} \
        -T AnalyzeCovariates \
        -R {params.refseq} \
        -before {input.before} \
        -after {input.after} \
        -plots {output.pdf}
        """

# merge lanes
# E01188-L2_S26_L005.sorted.bam E01188-L2_S26_L006.sorted.bam > E01188.sorted.merged.bam
def get_all_sorted_bams(samplename):
    return glob.glob("{0}*.sorted.bam".format(samplename))

rule merge_lanes:
    input: bams = lambda wildcards: get_all_sorted_bams(wildcards.sample), samtools = config['tools']['samtools']
    output: "{sample}.sorted.merged.bam"
    threads:
        12
    shell:
        "{input.samtools} merge -@ 12 {output} {input.bams}"

rule remove_duplicates:
    input:
        bam = config['datadirs']['bams'] + "/{sample}.sorted.merged.bam",
        java = config['tools']['java']
    output:
        bam = config['datadirs']['picard'] + "/{sample}_rmdup.bam"
    params:
        jar = config['jars']['markduplicates'],
        javaopts = config['tools']['javaopts']
    shell:
        # will (and need the permision to) create a tmp directory
        # with the name of login under specified tmp directory
        # Exception in thread "main" net.sf.picard.PicardException: Exception creating temporary directory.
        """
        {input.java} {params.javaopts} -jar {params.jar} \
        INPUT={input.bam} \
        OUTPUT={output.bam} \
        METRICS_FILE={config[datadirs][picard]}/{wildcards.sample}.txt
        """

rule make_index:
    input:
        bam = config['datadirs']['picard'] + "/{sample}.rmdup.bam",
        java = config['tools']['java']
    output:
        bai = config['datadirs']['picard'] + "/{sample}.rmdup.bai"
    params:
        javaopts = config['tools']['javaopts'],
        jar = config['jars']['buildbamindex']
    shell:
        """
        {input.java} {params.javaopts} -jar {params.jar} \
        INPUT={input.bam}
        """

rule add_readgroup:
    input:
        bam = config['datadirs']['picard'] + "/{sample}.rmdup.bam",
        bai = config['datadirs']['picard'] + "/{sample}.rmdup.bai",
        java = config['tools']['java']
    output:
        bam = config['datadirs']['picard'] + "/{sample}.group.bam"
    params:
        javaopts = config['tools']['javaopts'],
        jar = config['jars']['addorreplacereadgroups']
    shell:
        """
        {input.java} {params.javaopts} -jar {params.jar} \
        I={input.bam} \
        O={output.bam} \
        PL=illumina \
        LB={wildcards.sample} \
        PU={wildcards.sample} \
        SM={wildcards.sample}
        """

#### Variant Calling ####
rule make_gvcf:
    input:
        bam = config['datadirs']['recalibrated'] + "/{sample}.bam",
        java = config['tools']['java']
    output:
        gvcf = config['datadirs']['gvcfs'] + "/{sample}.gvcf"
    params:
        jar = config['jars']['gatk'],
        javaopts = config['tools']['javaopts'],
        refseq = config['refseq']
    shell:
        """
        {input.java} {params.javaopts} -jar {params.jar} \
        -T HaplotypeCaller \
        -R {params.refseq} \
        -I {input.bam} \
        --emitRefConfidence GVCF \
        --variant_index_type LINEAR \
        --variant_index_parameter 128000 \
        --genotyping_mode DISCOVERY \
        -stand_emit_conf 10 \
        -stand_call_conf 30 \
        -o {output.gvcf}
        """

rule combine_gvcfs:
    input: GVCFS,
           java = config['tools']['java']
    output:
        gvcf = config['datadirs']['gvcfs'] + "/joint.vcf"
    params:
        jar = config['jars']['gatk'],
        refseq = config['refseq'],
        list = GVCFSLIST,
        javaopts = config['tools']['javaopts'],
        db = config['siv']
    shell:
        """
        {input.java} {params.javaopts} -jar {params.jar} \
        -T GenotypeGVCFs \
        --dbsnp {params.db} \
        -nt 8 \
        -R {params.refseq} \
        {params.list} \
        -o {output.gvcf}
        """

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
    input: config['sample_table']
    output: config['pedfile']
    run:
        st = pandas.read_table("{0}".format(input))
        st['Sex']=st['Sex'].replace(['M','F'],[1,2])
        st['Affected_status']=st['Affected_status'].replace('unaffected',0)
        st['Affected_status']=st['Affected_status'].replace('[^0].+', 1, regex=True)
        ped = st[[0,1,3,2,4,5]]
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
        momdf['Affected_status']=0
        
        dadfams=ped[ped['Father'].isin(moms)]['FamilyID']
        daddf = pandas.DataFrame(dadfams)
        daddf['Subject']=dads
        daddf['Father']=0
        daddf['Mother']=0
        daddf['Sex']=1
        daddf['Affected_status']=0
        
        ped.append([momdf,daddf])
        
        ped.to_csv("{0}".format(output), sep='\t',index=False)

rule run_phase_by_transmission:
    input:
        vcf = config['datadirs']['gvcfs'] + "/joint.vcf",
        snpeff = config['datadirs']['gvcfs'] + "/snpeff.vcf",
        ped = config['ped'],
        java = config['tools']['java']
    output:
        vcf = config['datadirs']['gvcfs'] + "/phased.vcf",
        mvf = config['datadirs']['gvcfs'] + "/mendelian_violations.txt"
    params:
        jar  = config['jars']['gatk'],
        refseq = config['refseq'],
        javaopts = config['tools']['javaopts']
    log: 
        config['datadirs']['log'] + "/phase_by_transmission.log" 
    shell:
        """
        {input.java} {params.javaopts} -jar {params.jar} \
        -T PhaseByTransmission \
        -R {params.refseq} \
        -V {input.vcf} \
        -ped {input.ped} \
        -mvf {ouput.mvf} \
        -o {output.vcf} >& {log}
        """

#### Annotation ####
# ud - upstream downstream interval length (in bases)
rule run_snpeff:
    input:
        vcf = config['datadirs']['gvcfs'] + "/joint.vcf",
        java = config['tools']['java']
    output:
        vcf = config['datadirs']['gvcfs'] + "/snpeff.vcf"
    params:
        jar  = config['jars']['snpeff']['path'],
        conf = config['jars']['snpeff']['cnf'],
        javaopts = config['tools']['javaopts'],
        database = config['jars']['snpeff']['db'],
        updown = config['jars']['snpeff']['ud']
    shell:
        """
        {input.java} {params.javaopts} -jar {params.jar} \
        -c {params.conf} \
        -t {params.database} \
        -ud {params.updown} \
         {input.vcf} > {output.vcf}
        """

#### Report ####
# create YAML file used in meta-FastQC report
rule makeyaml:
    output: 
        yaml = "fastqc.yaml",
    params:
        projdir = config['projdir'],
        fastqc = config['datadirs']['fastqc']
    run:
        with open(output.yaml, "w") as out:
           idx = 1
           out.write("paired: yes\n") 
           out.write("output: {0}\n".format(params.projdir)) 
           out.write("fastqc:\n") 
           for name in BASENAMES:
               out.write("  {0}\n".format(name)) 
               out.write("  - {0}/{1}_R1_fastqc.zip\n".format(params.fastqc,name))
               out.write("  - {0}/{1}_R2_fastqc.zip\n".format(params.fastqc,name))

#### Create Markdown index of FastQC report files
rule makemd:
    output: 
        md = "fastqc.md",
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

#### Internal
onsuccess:
    print("Workflow finished, no error")
    shell("mail -s 'workflow finished' "+config['admins']+" < {log}")

onerror:
    print("An error occurred")
    shell("mail -s 'an error occurred' "+config['admins']+" < {log}")
