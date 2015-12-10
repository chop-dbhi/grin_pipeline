import glob
import re
from functools import cmp_to_key
"""
run on respublica
source activate snakeenv
snakemake -j -c "qsub -l h_vmem=40G -l mem_free=40G" 
"""

# ISIDIR = "/nas/is1"
ISIDIR = "/mnt/isilon/cbmi/variome"
TOOLDIR = ISIDIR + "/bin"
SCRIPTDIR = "zhangs3/scripts"
JAVA = TOOLDIR + "/java"
FASTQDIR = "fastq"
SEQDIR = "samples"
QCDIR = "fastqc"
SAMDIR = "sams"
BAMDIR = "bams"
GRCH38DIR = ISIDIR + "/reference/human/GRCh38"
SLINK = "GRIN/fastqc/"
PICARDDIR = "picard"
REALIGNDIR = "realigned"
REALIGNDIR2 = "realigned2"
REALIGNDIR3 = "realigned3"
RECALDIR = "recalibrated"
RECALDIR3 = "recalibrated3"
POSTRECALDIR = "postrecalibrated"
POSTRECALDIR3 = "postrecalibrated3"
PDFDIR = "pdfs"
PDFDIR3 = "pdfs3"
VCFDIR = "vcfs"

SAMTOOLS = TOOLDIR + "/Samtools/samtools-1.2/samtools"
PICARDJAR = TOOLDIR + "/picard-1.57.jar"
MarkDuplicates = TOOLDIR + "/picard/MarkDuplicates.jar"
BuildBamIndex = TOOLDIR + "/picard/BuildBamIndex.jar"
AddOrReplaceReadGroups = TOOLDIR + "/picard/AddOrReplaceReadGroups.jar"
GATKJAR = TOOLDIR + "/GATK/GenomeAnalysisTK.jar"
REFSEQFA = GRCH38DIR + "/genome.fa" # novoindex -k 14 -s 1
REFSEQIDX = GRCH38DIR + "/genome.idx" # novoindex -k 14 -s 1
GIV = GRCH38DIR + "/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz"
RTLIST = REALIGNDIR2 + "/realignment_targets.list"
BAMLIST = REALIGNDIR2 + "/bam.list"
SEQLIST = "files.lst"


MEMSET = "-Xmx20g -Xms10g"
JAVAOPTS = "-Djava.io.tmpdir=/mnt/lustre/users/leipzigj/"

TESTDIR = PICARDDIR

FASTQS = glob.glob(FASTQDIR + "/*.gz")
SAMPLES = ['E01612','E01613','E01614']
RESULTS = [re.sub("\.fastq\.gz$","_fastqc.zip", QCDIR + "/" + os.path.basename(name)) for name in FASTQS]
# HTMLS = [re.sub("\.fastq\.gz$","_fastqc.html", QCDIR + "/" + os.path.basename(name)) for name in FASTQS]

BASENAMES = [re.sub("\_R1\.fastq\.gz$", "", os.path.basename(name)) for name in FASTQS if re.search("\_R1\.fastq\.gz$", name)]

BASENAMES.sort()


SAMS = [SAMDIR + "/" + name + ".sam" for name in BASENAMES]
BAMS = [BAMDIR + "/" + name + ".bam" for name in BASENAMES]
DBAIS = [PICARDDIR + "/" + name + ".rmdup.bai" for name in BASENAMES]
DBAMS = [PICARDDIR + "/" + name + ".rmdup.bam" for name in BASENAMES]
GBAIS = [PICARDDIR + "/" + name + ".group.bai" for name in BASENAMES]
GBAMS = [PICARDDIR + "/" + name + ".group.bam" for name in BASENAMES]
RBAMS = [REALIGNDIR + "/" + name + ".bam" for name in BASENAMES]
LISTS = [REALIGNDIR + "/" + name + ".list" for name in BASENAMES]
RBAMS2 = [REALIGNDIR2 + "/" + name + ".bam" for name in BASENAMES]
RBAMS3 = [REALIGNDIR3 + "/" + name + ".bam" for name in BASENAMES]
TABLES = [RECALDIR + "/" + name + ".table" for name in BASENAMES]
RECBAMS = [RECALDIR + "/" + name + ".bam" for name in BASENAMES]
TABLES3 = [RECALDIR3 + "/" + name + ".table" for name in BASENAMES]
RECBAMS3 = [RECALDIR3 + "/" + name + ".bam" for name in BASENAMES]
POSTTABLES = [POSTRECALDIR + "/" + name + ".table" for name in BASENAMES]
POSTTABLES3 = [POSTRECALDIR3 + "/" + name + ".table" for name in BASENAMES]
PDFS = [PDFDIR + "/" + name + ".pdf" for name in BASENAMES]
PDFS3 = [PDFDIR3 + "/" + name + ".pdf" for name in BASENAMES]
VCFS = [VCFDIR + "/" + name + ".vcf" for name in BASENAMES]

INDELS = REALIGNDIR3 + "/indels.list"

SBAMS = [re.sub("\.bam$", ".sorted.bam", name) for name in BAMS]


dlocs = dict()

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
    "combine indels in the lists into a single list"
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

rule all:
    input: RESULTS, "fastqc.md", SBAMS, DBAIS, RBAMS

rule basehead:
    run:
        print(FASTQS)

rule dummy:    # just to test the python codes above
    input: "Snakefile"

rule print_reads3:  # with one combined indel list
    input: RECBAMS3

rule print_reads:
    input: RECBAMS

rule vcfs:
    input: VCFS

rule pdfs3:
    input: PDFS3

rule pdfs:
    input: PDFS

rule recalibrate3:
    input: POSTTABLES3

rule recalibrate:
    input: POSTTABLES

rule combine_indels:
    input: INDELS

rule realign:      # use individual indel list
    input: RBAMS   # realign bams

rule add_group:    # must be run before create_target
    input: GBAIS   # realign bais

rule make_bais:
    input: DBAIS

rule sortbams:
    input: DBAMS

#### run novoalign ####
rule novoalign:
    input:
        pair1 = FASTQDIR + "/{sample}_R1.fastq.gz",
        pair2 = FASTQDIR + "/{sample}_R2.fastq.gz"
    output:
        BAMDIR + "/{sample}.sam"
    threads:
        12
    shell:
        """
        {TOOLDIR}/novoalign -a -k -d {REFSEQIDX} -o SAM -f {input.pair1} {input.pair2} > {output}
        """

rule samtobam:
    input:
        sam = "{sample_path}.sam", samtools = SAMTOOLS
    output:
        temp("{sample_path}.temp.bam")
    threads:
        12
    shell:
        "{input.samtools} view -@ 12 -bS {input.sam} > {output}"

	
#### create yaml file ####
rule makeyaml:
    params:
        absfastqc = ISIDIR + "/" + QCDIR
    output: 
        ymal = "fastqc.yaml",
    run:
        with open(output.ymal, "w") as out:
           idx = 1
           out.write("paired: yes\n") 
           out.write("output: {0}\n".format(absfastqc)) 
           out.write("fastqc:\n") 
           for name in BASENAMES:
               out.write("  " + name + "\n") 
               out.write("  - " + absfastqc + "/" + name + "_R1_fastqc.zip\n")
               out.write("  - " + absfastqc + "/" + name + "_R2_fastqc.zip\n")

#### create md file ####
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
               out.write(" [R1]({{SLINK}}/fastqc/" + name + "_R1_fastqc.html)")
               out.write(" [R2]({{SLINK}}/fastqc/" + name + "_R2_fastqc.html)")
	       out.write("\n\n")
               idx += 1

#### QC ##### 
rule fastqc: 
    input: 
        pair1 = FASTQDIR + "/{sample}1.fastq.gz",
        pair2 = FASTQDIR + "/{sample}2.fastq.gz" 
    log: 
        "logs/{sample}.log" 
    output: 
        pair1 = QCDIR + "/{sample}1_fastqc.zip",
        pair2 = QCDIR + "/{sample}2_fastqc.zip",
        html1 = QCDIR + "/{sample}1_fastqc.html",
        html2 = QCDIR + "/{sample}2_fastqc.html"
    # how to run qsub?
    shell: "{TOOLDIR}/fastqc -o {QCDIR} {input.pair1} {input.pair2} 2> {log}" 



#### combine indel lists  ####

rule combine_lists:
    input: LISTS
    output: INDELS
    run:
         combine()

#### run GATK  ####

rule make_vcfs:
    input:
        bam = RECALDIR3 + "/{sample}.bam"
    output:
        vcf = VCFDIR + "/{sample}.vcf"
    shell:
        """
        {JAVA} {JAVAOPTS} {MEMSET} -jar {GATKJAR} \
        -T HaplotypeCaller \
        -R {REFSEQFA} \
        -I {input.bam} \
        -genotyping_mode DISCOVERY \
        -stand_emit_conf 10 \
        -stand_call_conf 30 \
        -o {output.vcf}
        """

rule recalibrated_result3:
    input:
        table = RECALDIR3 + "/{sample}.table",
        bam = REALIGNDIR3 + "/{sample}.bam"
    output:
        bam = RECALDIR3 + "/{sample}.bam"
    shell:
        """
        {JAVA} {JAVAOPTS} {MEMSET} -jar {GATKJAR} \
        -T PrintReads \
        -R {REFSEQFA} \
        -I {input.bam} \
        -BQSR {input.table} \
        -o {output.bam}
        """

rule recalibrated_result:
    input:
        table = RECALDIR + "/{sample}.table",
        bam = REALIGNDIR + "/{sample}.bam"
    output:
        bam = RECALDIR + "/{sample}.bam"
    shell:
        """
        {JAVA} {JAVAOPTS} {MEMSET} -jar {GATKJAR} \
        -T PrintReads \
        -R {REFSEQFA} \
        -I {input.bam} \
        -BQSR {input.table} \
        -o {output.bam}
        """


rule make_pdfs3:
    input:
        before = RECALDIR3 + "/{sample}.table",
        after = POSTRECALDIR3 + "/{sample}.table"
    output:
        pdf = PDFDIR3 + "/{sample}.pdf"
    shell:
        """
        source ~/.bashrc
        module load R/3.2.2
        {JAVA} {JAVAOPTS} {MEMSET} -jar {GATKJAR} \
        -T AnalyzeCovariates \
        -R {REFSEQFA} \
        -before {input.before} \
        -after {input.after} \
        -plots {output.pdf}
        """

rule make_pdfs:
    input:
        before = RECALDIR + "/{sample}.table",
        after = POSTRECALDIR + "/{sample}.table"
    output:
        pdf = PDFDIR + "/{sample}.pdf"
    shell:
        """
        module load R/3.2.2
        {JAVA} {JAVAOPTS} {MEMSET} -jar {GATKJAR} \
        -T AnalyzeCovariates \
        -R {REFSEQFA} \
        -before {input.before} \
        -after {input.after} \
        -plots {output.pdf}
        """

rule post_recalibrated_table3:
    input:
        table = RECALDIR3 + "/{sample}.table",
        bam = REALIGNDIR3 + "/{sample}.bam"
    output:
        table = POSTRECALDIR3 + "/{sample}.table"
    shell:
        """
        {JAVA} {JAVAOPTS} {MEMSET} -jar {GATKJAR} \
        -T BaseRecalibrator \
        -R {REFSEQFA} \
        -I {input.bam} \
        -knownSites {GIV} \
        -BQSR {input.table} \
        -o {output.table}
        """

rule post_recalibrated_table:
    input:
        table = RECALDIR + "/{sample}.table",
        bam = REALIGNDIR + "/{sample}.bam"
    output:
        table = POSTRECALDIR + "/{sample}.table"
    shell:
        """
        {JAVA} {JAVAOPTS} {MEMSET} -jar {GATKJAR} \
        -T BaseRecalibrator \
        -R {REFSEQFA} \
        -I {input.bam} \
        -knownSites {GIV} \
        -BQSR {input.table} \
        -o {output.table}
        """

rule recalibrated_table3:
    input:
        bam = REALIGNDIR3 + "/{sample}.bam"
    output:
        table = RECALDIR3 + "/{sample}.table"
    shell:
        """
        {JAVA} {JAVAOPTS} {MEMSET} -jar {GATKJAR} \
        -T BaseRecalibrator \
        -R {REFSEQFA} \
        -I {input.bam} \
        -knownSites {GIV} \
        -o {output.table}
        """

rule recalibrated_table:
    input:
        bam = REALIGNDIR + "/{sample}.bam"
    output:
        table = RECALDIR + "/{sample}.table"
    shell:
        """
        {JAVA} {JAVAOPTS} {MEMSET} -jar {GATKJAR} \
        -T BaseRecalibrator \
        -R {REFSEQFA} \
        -I {input.bam} \
        -knownSites {GIV} \
        -o {output.table}
        """
        # -knownSites dbsnp.vcf \


rule bam_list:       # create realign target list
    input:  # deduced bams
        #dbam = PICARDDIR + "/{sample}.rmdup.bam"
        dbam = GBAMS
    output:
        list = BAMLIST
    run:
        with open(output.list, "w") as out:
            for file in input.dbam:
                 out.write(file + "\n")

rule combined_target_list:       # create realign target list
    input:  # deduced bams
        bamlist = BAMLIST
    output:
        rtlist = RTLIST
    shell:
        """
        {JAVA} {JAVAOPTS} {MEMSET} -jar {GATKJAR} \
             -T RealignerTargetCreator \
             -R {REFSEQFA} \
             -I {input.bamlist} \
             -nt 16 \
             -known {GIV} \
             -o {output.rtlist}
        """

rule target_list:       # create realign target list
    input:  # deduced bams
        bam = PICARDDIR + "/{sample}.group.bam"
    output:
        list = REALIGNDIR + "/{sample}.list"
    shell:
        """
        {JAVA} {JAVAOPTS} {MEMSET} -jar {GATKJAR} \
             -T RealignerTargetCreator \
             -R {REFSEQFA} \
             -I {input.bam} \
             -known {GIV} \
             -o {output.list}
        """

rule realign_target3:   # with one combined list file
    input:  # deduced bams
        list = INDELS,
        dbam = PICARDDIR + "/{sample}.group.bai"
    output:
        rbam = REALIGNDIR3 + "/{sample}.bam"
    shell:
        """
        {JAVA} {JAVAOPTS} {MEMSET} -jar {GATKJAR} \
            -T IndelRealigner \
            -R {REFSEQFA} \
            -I {PICARDDIR}/{wildcards.sample}.group.bam \
            -targetIntervals {input.list} \
            -known {GIV} \
            -o {output.rbam}
        """

rule realign_target2:   # with one list from combined bams
    input:  # deduced bams
        list = RTLIST,
        dbam = PICARDDIR + "/{sample}.group.bai"
    output:
        rbam = REALIGNDIR2 + "/{sample}.bam"
    shell:
        """
        {JAVA} {JAVAOPTS} {MEMSET} -jar {GATKJAR} \
            -T IndelRealigner \
            -R {REFSEQFA} \
            -I {PICARDDIR}/{wildcards.sample}.group.bam \
            -targetIntervals {input.list} \
            -known {GIV} \
            -o {output.rbam}
        """
        # touch {REALIGNDIR}/{wildcards.sample}.bam 

rule realign_target:
    input:  # deduced bams
        list = REALIGNDIR + "/{sample}.list",
        dbam = PICARDDIR + "/{sample}.group.bai"
    output:
        rbam = REALIGNDIR + "/{sample}.bam"
    shell:
        """
        {JAVA} {JAVAOPTS} {MEMSET} -jar {GATKJAR} \
            -T IndelRealigner \
            -R {REFSEQFA} \
            -I {PICARDDIR}/{wildcards.sample}.group.bam \
            -targetIntervals {input.list} \
            -known {GIV} \
            -o {output.rbam}
        """
        # touch {REALIGNDIR}/{wildcards.sample}.bam 

#### run picard ####

rule make_group_index:
    input:
        bam = PICARDDIR + "/{sample}.group.bam"
    output:
        bai = PICARDDIR + "/{sample}.group.bai"
    shell:
        """
        {JAVA} {JAVAOPTS} {MEMSET} -jar {BuildBamIndex} INPUT={input.bam}
        """

rule add_readgroup:
    input:
        bam = PICARDDIR + "/{sample}.rmdup.bam"
    output:
        bam = PICARDDIR + "/{sample}.group.bam"
    shell:
        """
	{JAVA} {JAVAOPTS} {MEMSET} -jar {AddOrReplaceReadGroups} I={input.bam} O={output.bam} PL=illumina LB={wildcards.sample} PU={wildcards.sample} SM={wildcards.sample}
        """

rule make_index:
    input:
        bam = PICARDDIR + "/{sample}.rmdup.bam"
    output:
        bai = PICARDDIR + "/{sample}.rmdup.bai"
    shell:
        """
        {JAVA} {JAVAOPTS} {MEMSET} -jar {BuildBamIndex} INPUT={input.bam}
        """

rule remove_duplicates:
    input:
        bam = BAMDIR + "/{sample}.sorted.bam"
    output:
        bam = PICARDDIR + "/{sample}.rmdup.bam"
    shell:
        """
        {JAVA} {JAVAOPTS} {MEMSET} -jar {MarkDuplicates} INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={PICARDDIR}/{wildcards.sample}.txt
        """

#### run novosort ####
rule sorting:
    input:
        bam = BAMDIR + "/{sample}.temp.bam"
    output:
        sorted = BAMDIR + "/{sample}.sorted.bam",
    shell:
        """
        {TOOLDIR}/novosort -m 14g -t . --removeduplicates --keeptags -i -o {output.sorted} {input.bam}
        """



onsuccess:
    print("Workflow finished, no error")
    shell("mail -s 'Workflow finished!' leipzig@gmail.com < {log}")

