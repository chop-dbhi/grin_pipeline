
configfile: "config.json"

BIN_ROOT = "/nas/is1/bin"
ANNOVAR_DIR = BIN_ROOT + "/annovar"
ANNOVAR_PROG = ANNOVAR_DIR + "/annotate_variation.pl"
ANNODB = ANNOVAR_DIR + "/humandb"

ANNO_TARGETS = ['{0}/{1}_{2}.txt'.format(ANNODB,config["buildver"],anno_type) for anno_type in ['refGene','cytoBand','genomicSuperDups','esp6500siv2_all','esp6500siv2_ea','esp6500siv2_aa','exac03','kaviar_20150923','exac03nontcga','exac03nonpsych','hrcr1','1000g2014oct','snp138','dbnsfp30a']]

rule all:
  input: ANNO_TARGETS

rule download:
  input: annovar = ANNOVAR_PROG
  output: ANNODB + "/" + config["buildver"] + "_{type}.txt"
  params: buildver = config["buildver"], annodb = ANNODB
  shell:
    """
    {input.annovar} -buildver {params.buildver} -downdb -webfrom annovar {wildcards.type} {params.annodb}/ 
    """
