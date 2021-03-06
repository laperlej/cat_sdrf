# -*- coding: utf-8 -*-

import os.path
import operator
import re
from collections import OrderedDict
import sys
sys.path.insert(0, '../')
from cat_sdrf import saccer, pombe, celegans

# Assigne la valeur de row à l'index 'filename' comme étant 'csvfile.name'
PREPROCESS = lambda csvfile, row: operator.setitem(row, "filename", os.path.basename(csvfile.name))

def idf_extract(csvfile):
	if "seq.sdrf.txt" in csvfile.name:
		idf_file_path = csvfile.name.replace("seq.sdrf.txt", "idf.txt")
	else:
		idf_file_path = csvfile.name.replace("sdrf.txt", "idf.txt")
	idf_file = open(idf_file_path)
	for line in idf_file.readlines():
		if "Experiment Description" in line:
			expdescription = line.replace("Experiment Description	", "")
			return expdescription
		#else:
			#return "No description available"

PREPROCESS2 = lambda csvfile, row: operator.setitem(row, "expdescription", idf_extract(csvfile))
 
#itère sur chaque ligne et retourne vrai si remplit les 4 conditions.
def split_condition_aux(row, species):
	#les types d'essai qui nous intéressent
	#assays=["chip", "mnase", "dnase","other"]

	#Les types d'essais qui ne nous intéressent pas
	discard_assays=["rip-seq","rna-seq", "unwanted"]
	#dictionnaire des noms d'espèces
	species_dict={
		"saccer": "Saccharomyces cerevisiae",
		"pombe": "Schizosaccharomyces pombe",
		"celegans":"Caenorhabditis elegans"}
	
	return (species_dict[species] in row["3)organism"] and 
		   #[True for assay in assays if assay in row["4)assaytype"].lower()] and 
		   #"fastq" in row["12)fastq"] and
	 	   not [False for discard_assay in discard_assays if discard_assay in row["clean_assay"].lower()])

#selon l'espèce demandée (sys.argv[3]), appelle split_condition avec cette espèce
split_condition = {
	"saccer": lambda row: split_condition_aux(row, "saccer"),
	"pombe": lambda row: split_condition_aux(row, "pombe"),
	"celegans": lambda row: split_condition_aux(row, "celegans")}

#dictionnaire dans lequel les clés sont les espèces et les valeurs sont leur dictionnaire de gène
GENE_DICT = {
	"saccer": saccer.GENE_DICT,
	"pombe": pombe.GENE_DICT,
	"celegans": celegans.GENE_DICT}

GENE_DESCRIP_DICT = {
	"saccer": saccer.GENE_DESCRIP_DICT,
	"pombe": pombe.GENE_DESCRIP_DICT,
	"celegans": celegans.GENE_DESCRIP_DICT}

FIELDNAMES=OrderedDict([
	('1)identifier', lambda title, row: re.search('sourcename', title)),
	('2)filename', lambda title, row: re.search('filename', title)),
	('3)organism', lambda title, row: re.search("organism", title)),
	('clean_assay',lambda title, row: re.search('$a', title)),
	('clean_target', lambda title, row: re.search('$a', title)),
	('reliability', lambda title, row: re.search('$a', title)),
	('4)assaytype', lambda title, row: re.search('(characteristics\[sampledescription\]|\[iporinput\]|\[experimenttype\]|\[test\]|\[type\])',title)),
	('5)antibody', lambda title, row: re.search("(antibody|milliporecatno|vendor|label(?!ed)|antibodies)", title)),
	('6)target', lambda title, row: re.search('(epitopetag|tagged|taptag|protein|h2b|histone|immunoprecipitate|target|\[tag\]|\[pol\sgenotype\]|\[ip\])', title)),
	('7)treatment', lambda title, row: re.search('(phosphate|concentration|medium|media|condition|cycle|culturetype|transformedwith|treatment|temperature|percentage|compound|variable|spikedna|\[time\])', title)),
	('Material_type', lambda title, row: re.search('(cell\stype|organism\spart|sampletype|organelle|\[samplesubtype\]|materialtype|samplecomposition|\[fraction\]|tissue)', title)),
	('clean_celltype', lambda title, row: re.search('$a', title)),
	('cell_type',lambda title, row: re.search('(comment\[sample_source_name\]|\[age\]|cellline|growth|stage|developmental)', title)),
	('8)strain', lambda title, row: re.search('(strain|\[variant\])', title)),
	('9)genotype', lambda title, row: re.search('(genotype|genedeletion|variation\]|genetic|\[yrr1alleletransformed\]|background|\[pgal1\])', title)),
	('10)Technology', lambda title, row: re.search('(platform|instrument_model|technology)', title)),
	('11)description', lambda title, row: re.search('(comment\[sample_description\]|sample_characteristics|\[individual\]|comment\[sample_title\]|comment\[ena_alias\]|\[control\]|\[construct\])', title)),
	('12)fastq', lambda title, row: re.search('fastq_uri', title)),
	('Experiment description', lambda title, row: re.search('expdescription', title)),
	('file_description', lambda title, row: re.search('(arraydatafile|arrayexpress|\[submitted_file_name\]|normalization)', title)),
	('13)other', lambda title, row: re.search('.*', title))])

#dictionnaire des sortes d'essai
ASSAY_DICO = OrderedDict([
	('BrdU ChIP', 'brdu'),
	("MNase",'(mnase|monococcal\snuclease|chec\scleavage|chec\sexperiment|mononucleosomal\sdna)'),
	("ChIP",'(chip|chromatin.immunoprecipitation)'),
	("DNase",'dnase'),
	("rip-seq", 'rip-seq'),
	("IP", "(^ip|ip\s|ip-|\(ip\)|\sip)"),
	("RNA-seq", "(rna$|rna\s)"),
	('ssDNA', 'ssdna'),
	("other",'(other|dna|protein|whole-cell\sextract)'),
	("unwanted", '.*')])

#dictionnaire de marques d'histones
HISTONES_MARKS_DICO = OrderedDict([
	('RNAPII_ser5P', 's5\sphosphorylated\srna\spolii'),
	("DNAPIII",'(pol3|pol\s?iii)'),
	("DNAPII",'(pol2|pol\sɛ|pol\s?ii)'),
	("DNAPI",'(pol1|pol\sα|pol\s?i)'),
	("DNAP31",'(pol31|pol\sδ)'),
	('H3K4me1', '(h3k4me1|monomethylated\sh3k4|h3\s.?mono\smethyl\sK4)'),
	('H3K4me2', '(h3k4me2|dimethylated\sh3k4|h3\s.?di\smethyl\sK4)'),
	('H3K4me3', '(h3k4me3|trimethylated\sh3k4|h3\s.?tri\smethyl\sK4)'),
	('H3K14ac', '(h3k14ac)'),
	('H3K14', '(h3k14)'),
	('H3K36me3','(h3k36me3)'),
	('H3K36me2','(h3k36me2)'),
	('H3K36me','(h3k36me)'),
	('H3K36','(h3k36)'),
	('H3K56ac','(h3k56ac)'),
	('H3K56','(h3k56)'),
	('H3K79me3', 'h3k79me3'),
	('H3K79', 'h3k79'),
	('H3K9ac','(h3k9ac)'),
	('H3K9me3','(h3k9me3|h3.?k9.?me3)'),
	('H3K9me2','(h3k9me2|h3.?k9.?me2)'),
	('H3K9','h3k9'),
	('H3K4me3', '(h3.?k4.?me3)'),
	('H3K4','h3k4'),
	('H3K27me3', 'h3k27me3'),
	('H3R2me2', 'h3r2me2'),
	('H3', 'h3'),
	('H4K16ac','(h4k16ac)'),
	('H4K16','(h4k16)'),
	('H4K12ac', 'h4k12ac'),
	('H4K12', 'h4k12'),
	('H4K44ac','h4k44ac'),
	('H4K44','h4k44'),
	('H4K5ac','(h4k5ac)'),
	('H4K5','(h4k5)'),
	('H4K20me1','h4k20me1'),
	('H4','h4'),
	('H2A.Z','(htz1|h2a\.?z)'),
	('H2A', '(h2a)'),
	('H2B', '(htb1|htb2|h2b)')])
#dictionnaire des cibles des anticorps
TARGET_DICO=OrderedDict([
	('input','(input|whole\scell\sextract)'),
	('Negative control', '(negative\scontrol)'),
	('Mock', '(mock|no.?antibody)'),
	('H3K4me1', '(h3k4me1|monomethylated\sh3k4|h3\s.?mono\smethyl\sk4|ab8895)'),
	('H3K4me2', '(h3k4me2|dimethylated\sh3k4|h3\s.?di\smethyl\sk4)'),
	('H3K4me3', '(h3k4me3|trimethylated\sh3k4|h3\s.?tri\smethyl\sk4|ab8580|39159|305-34819|h3.?k4.?me3|ab8678)'),
	('H3K14ac', '(h3k14ac|07-353)'),
	('H3K14', '(h3k14)'),
	('H3K36me3','(h3k36me3|ab9050|300-95289)'),
	('H3K36me2','(h3k36me2)'),
	('H3K36me','(h3k36me)'),
	('H3K36','(h3k36)'),
	('H3K56ac','(h3k56ac|07-677)'),
	('H3K56','(h3k56)'),
	('H3K79me3', '(h3k79me3)'),
	('H3K79', '(h3k79)'),
	('H3K9ac','(h3k9ac|06-942)'),
	('H3K9me3','(h3k9me3|h3.?k9.?me3|ab8898)'),
	('H3K9me2','(h3k9me2|h3.?k9.?me2)'),
	('H3K9-14ac','(ach3\sk9,14|h3k9-14ac)'),
	('H3K9','(h3k9)'),
	('H4K16ac','(h4k16ac|07-329)'),
	('H4K16','(h4k16)'),
	('H4K12ac', 'h4k12ac'),
	('H4K12', 'h4k12'),
	('H4K44ac','h4k44ac'),
	('H4K44','(h4k44)'),
	('H4K5ac','(h4k5ac|07-327)'),
	('H4K5','(h4k5)'),
	('H4K20me1','(h4k20me1|ab9051)'),
	('H4ac','(h4ac|39177)'),
	('H3K4','h3k4'),
	('H3K27me3', 'h3k27me3'),
	('H3','(h3|ab12079|05-928|07-690|ab1791)'),
	('H4','(h4|ab7311)'),
	('H2A.Z','(htz1|h2a\.?z|ab4626)'),
	('H2A', '(h2a|ab13923)'),
	('H2B', '(htb1|htb2|h2b|ab1790)'),
	('ESA1','(esa1)'),
	('Q105me1','(q105me1)'),
	('AAR','(aar|o-acetyl-adp-ribose)'),
	('DNMT3b','(dnmt3b|ab2851)'),
	('RNAPII_tyr1P','(tyr1|61383|mabe350)'),
	('RNAPII_ser2P','(ser2p|ser2-po4|ab24758|ab193468|ab5095|s2\sphosphorylated\srna\spolii)'),
	('RNAPII_ser5P','(4h8|ab5408|ser5p|ab55208|ab140748|ab5401|ab193467|ab5131|s5\sphosphorylated\srna\spolii)'),
	('RNAPII_ser7P','(ser7p|ab126537)'),
	('RNAPII_CTD','(ctd|8wg16|ab817|mms-126r-200)'),
	('RNAPII_RPB3','(wp012|1Y26|ab202893|rpb3)'),
	('POL30','pcna'),
	('RAD51','rad51'),
	('ORC','orc'),
	('TAF7','ptr6'),
	('SIR2', '(sir2|dam1514081|07131|sirt1)'),
	('Mcm2-7','(mcm2-7|um185)'),
	('RNAPIII','(rpc1|53330002|rnapiii|pol.?3|pol.?i{3})'),
	('RNAPII','(rnapii|pol.?2|pol.?i{2}|rna\spoly?m?e?r?a?s?e?\si?i?2?)'),
	('tag_myc','(myc|05-419|9e10|9e11|ab56|dam1724025)'),
	('tag_HA','(^ha|ha$|anti.ha|ha11|12ca5|ab16918)'),
	('tag_PK','(v5|sv5-pk1|pk|mca1360)'),
	('tag_flag','(flag|f1804)'),
	('tag_T7', 't7'),
	('tag_GFP', '(gfp|egfp|ab290)'),
	('tag_tap','(tap|igg\ssepharose|cab1001)'),
	('tag_brdu', 'brdu'),
	('Mock_IgG','igg'),
	("RNA/DNA hybrid",'(rna/dna\shybrid)'),
	('none','(none|n/a)'),
	('empty', '.*')])

ANTIBODY_DICO = OrderedDict ([
	('H3K4me1', '(ab8895)'),
	('H3K4me3', '(ab8580|39159|305-34819|ab8678)'),
	('H3K14ac', '(07-353)'),
	('H3K36me3','(ab9050|300-95289)'),
	('H3K56ac','(07-677)'),
	('H3K9ac','(06-942)'),
	('H3K9me3', '(ab8898)'),
	('H4K16ac','(07-329)'),
	('H4K5ac','(07-327)'),
	('H4K20me1','(ab9051)'),
	('H4ac','(39177)'),
	('H3','(ab12079|05-928|07-690|ab1791)'),
	('H4','(ab7311)'),
	('H2A', '(ab13923)'),
	('H2A.Z','(ab4626)'),
	('H2B', '(ab1790)'),
	('DNMT3b','(ab2851)'),
	('SIR2', '(dam1514081|07131)'),
	('RNAPIII','(53330002)'),
	('RNAPII_tyr1P','(61383|mabe350)'),
	('RNAPII_ser2P','(ab24758|ab193468|ab5095)'),
	('RNAPII_ser5P','(4h8|ab5408|ab55208|ab140748|ab5401|ab193467|ab5131)'),
	('RNAPII_ser7P','(ab126537)'),
	('RNAPII_CTD','(8wg16|ab817|mms-126r-200)'),
	('RNAPII_RPB3','(wp012|1Y26|ab202893)')])

#dictionnaire des tag et leur regex pour la cible taggée
TAG_DICO=OrderedDict([
	 ('tag_HA','((\w+)\.ha|(\w+)::\d?ha|(\w+)-ha|ha.tagged.(\w+)|ha-(\w+)|\|ha::(\w+)|(\w+)::\S*ha|(\w+-\d+)_ha|(\w+)-\S*ha|ha\S*::(\w+))'),
	 ('tag_GFP','((\w+)\.gfp|(?:anti)(\w+)-gfp|gfp::(\w+)|(\w+-?\d+)_gfp|(\w+\.\d+)_gfp|(\w+)-\S*gfp|gfp\S*::(\w+)|gfp-tagged\s(\w+-?\d+)|(\w+.\d+)::?ty1\se?gfp|(\w+.?\d*)::\S*gfp)'),
	 ('tag_flag','((\w+)\.flag|(\w+)::flag|(\w+)::\S*flag|(\w+)-flag|(\w+-\d+)_flag|(\w+)-\S*flag|flag-(\w+)|flag::(\w+)|flag\S*-(\w+)|flag\S*::(\w+)|(\w+)_flag|flag.tagged\s(\w+))'),
	 ('tag_myc','((\w+)\.myc|(\w+)::myc|(\w+)-myc|myc.tagged\s(\w+)|myc-(\w+)|myc::(\w+)|(\w+)::\w*myc|(\w+-\d+)_myc|(\w+)-\S*myc|myc\S*-(\w+)|myc\S*::(\w+)|(\w+)myc)'),
	 ('tag_PK','((\w+)::\S*pk|pk::(\w+)|(\w+)\s?pk|(\w{2,})-pk|(\w+)-\S*pk|pk\S*-(\w+)|v5-(\w+))'),
	 ('tag_tap','((\w+)\.tap|(\w+)::tap|(\w+)-tap|(\w+)-\w{2,}-tap|(\w+-\d+)_tap|tap::(\w+)|(\w+)::\w*tap|tap\S*-(\w+)|tap\S*::(\w+)|(\w+)\stap|(\w+)tap)'),
	 ('tag_T7','((\w+)::\S*t7|(\w+)-\S*t7|t7\S*-(\w+)|t7\S*::(\w+))'),
	 ('tag_brdu','((\w+)::brdu|(\w+)-brdu|brdu-(\w+)|brdu::(\w+))') ])

#dictionnaire utilisant le mot-clé chip pour trouver la protéine-cible de l'essai
CHIP_DICO = OrderedDict ([
	('protein chip','((\w{2,})\sprotein\schip|(\w+\s\w+)\sprotein\schip)'),
	('BrdU IP','((\w{2,})\sbrdu\sip|(\w+)\s\w*\sbrdu\sip)'),
	('chip','((\w{2,})\s\s?chip|(\w{2,}\.\w+)\schip)'),
	('_chip','(\w{2,})_chip'),
	('chromatin immunoprecipitation', '(\w{2,})\s(?:wildtype|\w+)\snative\schromatin\simmunoprecipi?t?ation'),
	('chromatin ip against', '(chromatin\sip\sagainst\s(\w+\W\d+)|chromatin\sip\sagainst\s(\w{2,}))'),
	('something chromatin input', '(\w+)\schromatin\sinput'),
	('IP','(\w{2,}_(\w+)_ip|(\w{2,})\s?ip|(\w+)_ip)'),
	('IP of', '(ip\sof\s(\w+))'),
	('something IP','((\w{2,})_ip|(\w{2,})\s\w*\sip)'),#vague
	('something chip','((\w{2,})\s\w*\schip|(\w+-\d+)\schip)') #peut être problématique car vague et parfois 2 options
	])
	 
#dictionnaire des les types cellulaires utilisés pour le ChIP-Seq
CELL_TYPE = OrderedDict ([
	('L1 arrest', 'l1\sarrest'),
	('L1', '(\s(l1)|\|(l1)|^l1|_(l1))'),
	('L2-L3', 'l2-l3'),
	('L2', '(\s(l2)|\|(l2)|^l2|_(l2))'),
	('L3', '(\s(l3)|\|(l3)|^l3|_(l3))'),
	('L4-young adult', 'l4/young\sadult'),
	('L4', '(\s(l4)|\|(l4)|^l4|_(l4))'),
	('Late embryo', 'late\sembryo'),
	('Early embryo', 'early\sembryo'),
	('Embryo', '(embryo|mxemb)'),
	('Dauer', 'dauer'),
	('Young adult', '(young|young\sadult)'),
	('Adult','(adult|old)'),
	('Total nuclei', 'total\snuclei'),
	('Muscle nuclei', 'muscle\snuclei'),
	('Oocytes', 'oocytes'),
	('DIP DNA', 'dip\sdna'),
	('ssDNA', 'ssdna'),
	('Genomic DNA', 'genomic\sdna'),
	('DNA','dna'),
	('mRNA', 'mrna'),
	('RNA', 'rna'),
	('Protein', 'protein'),
	('Whole-cell extract', 'whole.cell.extract'),
	('Not specified', '.*')])
