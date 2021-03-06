# -*- coding: utf-8 -*-

import os.path
import operator
import re
from collections import OrderedDict
import saccer, pombe, celegans
from antibody_filter import merge_cols

# Assign 'csvfile.name' as the value of 'row' at the index 'filename' 
PREPROCESS_1 = lambda csvfile, row: operator.setitem(row, "filename", os.path.basename(csvfile.name))

#gets information from each associated idf file
def idf_extract(csvfile, keyword_list):
	if "seq.sdrf.txt" in csvfile.name:
		idf_file_path = csvfile.name.replace("seq.sdrf.txt", "idf.txt")
	else:
		idf_file_path = csvfile.name.replace("sdrf.txt", "idf.txt")
	idf_file = open(idf_file_path)
	for line in idf_file.readlines():
		output1 = ""
		#Search first keyword in the lines. Priority on first keyword, most specific
		if keyword_list[0] in line:
			output1 = line.replace(str(keyword_list[0]), "")
			#returns output1 providing that it is not an empty string
			if output1.strip() is not "" and output1.strip() is not None:
				return output1.replace('\t', ' ')
				break
	#if first keyword not in any lines, use second keyword, providing that there is a second one
	else:
		for line in idf_file.readlines():
			output2 = ""
			if len(keyword_list)== 2 and keyword_list[1] in line:	
				output2 = line.replace(str(keyword_list[1]), "")
				return output2.replace('\t', ' ')
		
#Assign the info from the 'experiment description' line in associated idf file as the value of 'row' at the index 'experiment_description_idf'
PREPROCESS_2 = lambda csvfile, row: operator.setitem(row, "experiment_description_idf", idf_extract(csvfile, ['Experiment Description']))
#Assign the info from the 'PubMed ID' line in associated idf file as the value of 'row' at the index 'pubmed_id_idf'
PREPROCESS_3 = lambda csvfile, row: operator.setitem(row, "pubmed_id_idf", idf_extract(csvfile, ['PubMed ID']))
#Assign the info from the 'Publication Author list' line in associated idf file as the value of 'row' at the index 'author_list_idf'
PREPROCESS_4 = lambda csvfile, row: operator.setitem(row, "author_list_idf", idf_extract(csvfile, ['Publication Author List','Person Last Name']))
#Assign the info from the 'submission date' line in associated idf file as the value of 'row' at the index 'submission_date_idf'
PREPROCESS_5 = lambda csvfile, row: operator.setitem(row, "submission_date_idf", idf_extract(csvfile, ['Comment[ArrayExpressSubmissionDate]', 'Comment[SRASubmissionDate]' ]))
#Assign the info from the 'Public release date' line in associated idf file as the value of 'row' at the index 'release_date_idf'
PREPROCESS_6 = lambda csvfile, row: operator.setitem(row, "release_date_idf", idf_extract(csvfile, ['Public Release Date']))
#Assign the info from the 'Protocol description' line in associated idf file as the value of 'row' at the index 'protocol_description_idf'
PREPROCESS_7 = lambda csvfile, row: operator.setitem(row, "protocol_description_idf", idf_extract(csvfile, ["Protocol Description"]))

#iterate on each ine and return True the conditions are met.
def split_condition_aux(row, species):
	#Assays type to discard
	discard_assays=["rip-seq","rna-seq", "unwanted", 'non-genomic', 'wgs']
	#file types needed
	file_types = ['.fastq.gz', '.bam', '.sam', 'SRX']
	#dictionnary with short name (sys.argv[3]) and full name of the species 
	species_dict={
		"saccer": "Saccharomyces cerevisiae",
		"pombe": "Schizosaccharomyces pombe",
		"celegans":"Caenorhabditis elegans"
		}
	
	return (species_dict[species] in row["3)organism"] and 
		   #[True for assay in assays if assay in row["4)assaytype"].lower()] and 
		   any(file_type in row["18)raw_files"] for file_type in file_types) and
			#not "non-genomic" in row["Material_type"] and
	 	   	not [False for discard_assay in discard_assays if discard_assay in row["4)clean_assay"].lower()])

#Depending on species requested (sys.argv[3]), calls split_condition with the right species
split_condition = {
	"saccer": lambda row: split_condition_aux(row, "saccer"),
	"pombe": lambda row: split_condition_aux(row, "pombe"),
	"celegans": lambda row: split_condition_aux(row, "celegans")
	}

GENE_DICT = {
	"saccer": saccer.GENE_DICT,
	"pombe": pombe.GENE_DICT,
	"celegans": celegans.GENE_DICT
	}

GENE_DESCRIP_DICT = {
	"saccer": saccer.GENE_DESCRIP_DICT,
	"pombe": pombe.GENE_DESCRIP_DICT,
	"celegans": celegans.GENE_DESCRIP_DICT
	}

CELL_TYPE_DICT = {
	"saccer": saccer.CELL_TYPE,
	"pombe": pombe.CELL_TYPE,
	"celegans": celegans.CELL_TYPE
	}
#only keys are useful in this dictionnary, the rest won't serve
FIELDNAMES=OrderedDict([
	('1)identifier', lambda title, row: re.search('sourcename', title)),
	('1,1)Sample_title', lambda title, row: re.search('$a', title)),
	('2)filename', lambda title, row: re.search('filename', title)),
	('3)organism', lambda title, row: re.search("\[organism\]", title)),
	('4)clean_assay',lambda title, row: re.search('$a', title)),
	('5)clean_target', lambda title, row: re.search('$a', title)),
	('6)reliability', lambda title, row: re.search('$a', title)),
	('7)assaytype', lambda title, row: re.search('(\[library.?selection\]|\[library.?strategy\]|characteristics\[sampledescription\]|\[iporinput\]|\[experimenttype\]|\[test\]|\[type\]|\[protocoltype\])',title)),
	('8)antibody', lambda title, row: re.search("(antibody|milliporecatno|vendor|\[label\|label$|antibodies)", title)),
	('9)target', lambda title, row: re.search('(epitopetag|tagged|taptag|protein|h2b|histone|immunoprecipitate|target|\[tag\]|\[pol\sgenotype\]|\[ip\])', title)),
	('10)treatment', lambda title, row: re.search('(\[mnasedigestiontime\]|\[mnaseorexoiii\]|phosphate|concentration|medium|media|condition|cycle|culturetype|transformedwith|treatment|temperature|percentage|compound|spikedna|\[time\])', title)),
	('11)Material_type',lambda title, row: re.search('(\[library_source\]|sampletype|materialtype|samplecomposition|\[samplesubtype\]|\[fraction\]|tissue|cell\stype|organismpart|organelle)', title)),
	('12)clean_celltype', lambda title, row: re.search('$a', title)),
	('13)cell_type',lambda title, row: re.search('(comment\[sample_source_name\]|\[age\]|cellline|growth|stage|developmental)', title)),
	('14)strain', lambda title, row: re.search('(strain|\[variant\])', title)),
	('15)genotype', lambda title, row: re.search('(genotype|genedeletion|variation\]|genetic|\[yrr1alleletransformed\]|background)', title)),
	('16)platform', lambda title, row: re.search('(platform|instrument_model)', title)),
	('17)Sample_description', lambda title, row: re.search('(comment\[sample_description\]|sample_characteristics|\[individual\]|comment\[sample_title\]|comment\[ena_alias\]|\[control\]|variable)', title)),
	('18)raw_files', lambda title, row: re.search('$a', title)),
	('19)all_supp_files', lambda title, row: re.search('fastq_uri', title)),
	('20)SRA_files', lambda title, row: re.search('(comment\[ena_run\]|comment\[ena_experiment\])', title)),
	('21)Experiment description', lambda title, row: re.search('experiment_description_idf', title)),
	('22)Protocol', lambda title, row: re.search('protocol_description_idf', title)),
	('23)Author(s)', lambda title, row: re.search('author_list_idf', title)),
	('24)Submission Date', lambda title, row: re.search('submission_date_idf', title)),
	('25)Release Date', lambda title, row: re.search('release_date_idf', title)),
	('26)Pubmed ID', lambda title, row: re.search('pubmed_id_idf', title)),
	('Other', lambda title, row: re.search('.*', title))])

#assay type dictionnary (WGS:Whole Genome Shotgun sequencing)
ASSAY_DICO = OrderedDict([
	('Non-genomic', '(non.genomic|transcriptomic|total\srna|nascent\srna|polya\srna)'),
	('BrdU-ChIP', '(brdu\sip|brdu-ip)'),
	('ATAC-Seq', 'atac-seq'),
	('BrdU', 'brdu'),
	("WGS", 'wgs'),
	('ChIP-exo', 'chip-exo'),
	("ChIP-eSPAN", "chip-espan"),
	('FAIRE-Seq', 'faire'),
	("ChIP-Seq",'(chip|chip-seq|chromatin\simmunoprecipitation|immunoprecipitation\sof\snative\schromatin)'),
	("MNase-Seq",'(mnase|monococcal\snuclease|micrococcal\snuclease|chec\scleavage|chec\sexperiment|mononucleosomal\sdna|micro-c)'),
	("DNase-Seq",'dnase'),
	("Bisulfite-Seq", "bisulfite"),
	("Rip-Seq", 'rip-seq'),
	("RNA-Seq", 'rna-seq'),
	('Okazaki fragment', 'okazaki\sfragment'),
	("other",'(other|genomic\sdna)'),
	("unwanted", '.*')
	])

#Histone marks and polymerase dictionnary
HISTONES_MARKS_DICO = OrderedDict([
	('RNAPII_tyr1P','(61383|mabe350)'),
	('RNAPII_ser2P','(ser2p|ser2-po4|ab24758|ab193468|ab5095|s2\sphosphorylated\srna\spolii)'),
	('RNAPII_ser5P','(4h8|ab5408|ser5p|ab55208|ab140748|ab5401|ab193467|ab5131|s5\sphosphorylated\srna\spolii)'),
	('RNAPII_ser7P','(ser7p|ab126537)'),
	('RNAPII_CTD','(8wg16|ab817|mms-126r-200|rpb1)'),
	('RNAPII_RPB3','(wp012|1Y26|ab202893|rpb3)'),
	("POLIII",'(pol3$|pol3\s|pol3-|pol\s?iii)'),
	("POLII",'(pol2|pol\sɛ|\\u025b|pol\s?ii)'),
	("POLI",'(pol1|pol\sα|pol\s?i|\\u03b1)'),
	("POL31",'(pol31|pol\sδ|\\u03b4)'),
	('POL', '(pol$|pol\s)'),
	('H3K4me1', '(h3k4me1|monomethylated\sh3k4|h3\s.?mono\smethyl\sk4|ab8895)'),
	('H3K4me2', '(h3k4me2|dimethylated\sh3k4|h3\s.?di\smethyl\sk4)'),
	('H3K4me3', "(h3k4me3|trimethylated\sh3k4|h3\s.?tri\smethyl\sk4|ab8580|39159|305-34819|h3.?k4.?me3|ab8678)"),
	('H3K4ac', '(h3k4ac|07-539)'),
	('H3K4','(h3k4)'),
	('H3K9ac','(h3k9ac|06-942)'),
	('H3K9me3','(h3k9me3|h3.?k9\s?me3|ab8898)'),
	('H3K9me2','(h3k9me2|h3.?k9\s?me2)'),
	('H3K9','(h3k9)'),
	('H3K14ac', '(h3k14ac|07-353)'),
	('H3K14', '(h3k14)'),
	('H3K18ac', '(h3k18ac|ab1191)'),
	('H3K23ac', '(h3k23ac|07-355)'),
	('H3K27ac', '(h3k27ac|07-360)'),
	('H3K27me3', '(h3k27me3)'),
	('H3K36me3','(h3k36me3|ab9050|300-95289)'),
	('H3K36me2','(h3k36me2)'),
	('H3K36me','(h3k36me)'),
	('H3K36','(h3k36)'),
	('H3K56ac','(h3k56ac|07-677)'),
	('H3K56','(h3k56)'),
	('H3K79me3', '(h3k79me3|trimethylated\sh3k79|ab2621)'),
	('H3K79me', '(h3k79me|ab2886)'),
	('H3K79', 'h3k79'),
	('H3R2me2', '(h3r2me2)'),
	('H3K9-14ac','(ach3\sk9,14|h3k9-14ac)'),
	('H3','(h3|ab12079|05-928|07-690|ab1791)'),
	('H4K16ac','(h4k16ac|07-329)'),
	('H4K16','(h4k16)'),
	('H4K12ac', '(h4k12ac)'),
	('H4K12', '(h4k12)'),
	('H4K44ac','(h4k44ac)'),
	('H4K44','(h4k44)'),
	('H4K5ac','(h4k5ac|07-327)'),
	('H4K5','(h4k5)'),
	('H4K8ac', '(h4k8ac|07-328)'),
	('H4K20me1','(h4k20me1|ab9051)'),
	('H4R3me2s','(h4r3me2s|ab5823)'),
	('H4R3me','(h4r3me|ab17339)'),
	('H4ac','(h4ac|39177)'),
	('H4','(h4|ab7311)'),
	('H2A.Z','(htz1|h2a\.?z|pht1|ab4626)'),
	('H2AK5ac', '(h2ak5ac|ab45152)'),
	('H2A', '(h2a|ab13923)'),
	('H2B', '(htb1|htb2|h2b|ab1790)')])
#Target and antibody dictionnary (won't be needed anymore)
TARGET_DICO=OrderedDict([
	('input','(input|input_dna|whole\scell\sextract|none|n/a)'),
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
	('H3K9me3','(h3k9me3|h3.?k\s?me3|ab8898)'),
	('H3K9me2','(h3k9me2|h3.?k9\s?me2)'),
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
	('H2A.Z','(htz1|h2a\.?z|pht1|ab4626)'),
	('H2A', '(h2a|ab13923)'),
	('H2B', '(htb1|htb2|h2b|ab1790)'),
	('ESA1','(esa1)'),
	('H2AQ105me1','(q105me1)'),
	('AAR','(aar|o-acetyl-adp-ribose)'),
	('DNMT3b','(dnmt3b|ab2851)'),
	('RNAPII_tyr1P','(tyr1|61383|mabe350)'),
	('RNAPII_ser2P','(ser2p|ser2-po4|ab24758|ab193468|ab5095|s2\sphosphorylated\srna\spolii)'),
	('RNAPII_ser5P','(4h8|ab5408|ser5p|ab55208|ab140748|ab5401|ab193467|ab5131|s5\sphosphorylated\srna\spolii)'),
	('RNAPII_ser7P','(ser7p|ab126537)'),
	('RNAPII_CTD','(ctd|8wg16|ab817|mms-126r-200|rpb1)'),
	('RNAPII_RPB3','(wp012|1Y26|ab202893|rpb3)'),
	('POL30','pcna'),
	('RAD51','rad51'),
	('ORC','orc'),
	('TAF7','ptr6'),
	('SIR2', '(sir2|dam1514081|07131|sirt1)'),
	('Mcm2-7','(mcm2-7|um185)'),
	("POLIII",'(pol3$|pol3\s|pol3-|pol\s?iii)'),
	("POLII",'(pol2|pol\sɛ|pol\s?ii)'),
	("POLI",'(pol1|pol\sα|pol\s?i)'),
	('RNAPIII','(rpc1|53330002|rnapiii|rnap3)'),
	('RNAPII','(rnapii|rnap2|rna\spoly?m?e?r?a?s?e?\si?i?2?)'),
	('tag_myc','(myc|05-419|9e10|9e11|ab56|dam1724025)'),
	('tag_HA','(^ha|ha$|ha\s|anti.ha|ha11|12ca5|ab16918)'),
	('tag_PK','(v5|sv5-pk1|pk|mca1360)'),
	('tag_flag','(flag|f1804)'),
	('tag_T7', 't7'),
	('tag_GFP', '(gfp|egfp|ab290)'),
	('tag_tap','(tap|igg\ssepharose|cab1001)'),
	('Mock_IgG','igg'),
	("RNA/DNA hybrid",'(rna/dna\shybrid)'),
	('empty', '.*') ])
#target dictionnary, not redundant with the histone dictionnary
TARGET2 = OrderedDict ([
	('input','(input|input_dna|whole\scell\sextract|none|n/a)'),
	('Negative control', '(negative\scontrol)'),
	('Mock', '(mock|no.?antibody)'),
	('ESA1','(esa1)'),
	('H2AQ105me1','(q105me1)'),
	('AAR','(aar|o-acetyl-adp-ribose)'),
	('DNMT3b','(dnmt3b|ab2851)'),
	('RNAPII_tyr1P','(tyr1|61383|mabe350)'),
	('RNAPII_ser2P','(ser2p|ser2-po4|ab24758|ab193468|ab5095|s2\sphosphorylated\srna\spolii)'),
	('RNAPII_ser5P','(4h8|ab5408|ser5p|ab55208|ab140748|ab5401|ab193467|ab5131|s5\sphosphorylated\srna\spolii)'),
	('RNAPII_ser7P','(ser7p|ab126537)'),
	('RNAPII_CTD','(ctd|8wg16|ab817|mms-126r-200|rpb1)'),
	('RNAPII_RPB3','(wp012|1Y26|ab202893|rpb3)'),
	('POL30','pcna'),
	('RAD51','rad51'),
	('ORC','orc'),
	('TAF7','ptr6'),
	('SIR2', '(sir2|dam1514081|07131|sirt1)'),
	('Mcm2-7','(mcm2-7|um185)'),
	('RNAPIII','(rpc1|53330002|rnapiii|rnap3)'),
	('RNAPII','(rnapii|rnap2|rna\spoly?m?e?r?a?s?e?\si?i?2?)'),
	('tag_myc','(myc|05-419|9e10|9e11|ab56|dam1724025)'),
	('tag_HA','(^ha|ha$|ha\s|anti.ha|ha11|12ca5|ab16918)'),
	('tag_PK','(v5|sv5-pk1|pk|mca1360)'),
	('tag_flag','(flag|f1804)'),
	('tag_T7', 't7'),
	('tag_GFP', '(gfp|egfp|ab290)'),
	('tag_tap','(tap|igg\ssepharose|cab1001)'),
	('Mock_IgG','igg'),
	("RNA/DNA hybrid",'(rna/dna\shybrid)'),
	('empty', '.*')])
#catalog number dictionnary
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
	('RNAPII_RPB3','(wp012|1Y26|ab202893)') ])

#dictionnaire des tag et leur regex pour la cible taggée
TAG_DICO_old=OrderedDict([
	 ('tag_HA','((\w+)\.ha|(\w+)::ha|(\w+)-ha|ha-(\w+)|\|ha::(\w+)|(\w+)::\S*ha|(\w+-\d+)_ha|(\w+)-\S*ha|ha\S*::(\w+))'),
	 ('tag_GFP','((\w+)\.gfp|(?:anti)(\w+)-gfp|gfp::(?:3xflag)(\w+)|(\w+-?\d+)_gfp|(\w+\.\d+)_gfp|(?:anti)(\w+)-\S*gfp|gfp\S*::(?:3xflag)(\w+)|gfp-tagged\s(\w+-?\d+)|(\w+.\d+)::?ty1\se?gfp|(\w+.?\d*)::\S*gfp)'),
	 ('tag_flag','((\w+)\.flag|(\w+)::flag|(\w+)::\S*flag|(\w+)-flag|(\w+-\d+)_flag|(\w+)-\S*flag|flag-(\w+)|flag::(\w+)|flag\S*-(\w+)|flag\S*::(\w+)|(\w+)_flag|flag-tagged\s(\w+))'),
	 ('tag_myc','((\w+)\.myc|(\w+)::myc|(\w+)-myc|myc-(\w+)|myc::(\w+)|(\w+)::\w*myc|(\w+-\d+)_myc|(\w+)-\S*myc|myc\S*-(\w+)|myc\S*::(\w+)|(\w+)\smyc|(\w+)myc)'),
	 ('tag_PK','((\w+)::\S*pk|pk::(\w+)|(\w+)\s?pk|(\w{2,})-pk|(\w+)-\S*pk|pk\S*-(\w+)|v5-(\w+))'),
	 ('tag_tap','((\w+)\.tap|(\w+)::tap|(\w+)-tap|(\w+)-\w{2,}-tap|(\w+-\d+)_tap|tap::(\w+)|(\w+)::\w*tap|tap\S*-(\w+)|tap\S*::(\w+)|(\w+)\stap|(\w+)tap)'),
	 ('tag_T7','((\w+)::\S*t7|(\w+)-\S*t7|t7\S*-(\w+)|t7\S*::(\w+))') ])

#dictionnaire des tag et leur regex pour la cible taggée
TAG_DICO=OrderedDict([
	 ('tag_HA','((\w+)\.ha|(\w+)::\d?ha|(\w+)-ha|ha.tagged.(\w+)|ha-(\w+)|\|ha::(\w+)|(\w+)::\S{1,5}ha|(\w+-\d+)_ha|(\w+)-\S{1,5}ha|ha\S{1,5}::(\w+))'),
	 ('tag_GFP','((\w+)\.gfp|(?:anti)(\w+)-gfp|gfp::(\w+)|(\w+-?\d+)_gfp|(\w+\.\d+)_gfp|(\w+)-\S*gfp|gfp\S{1,5}::(\w+)|gfp-tagged\s(\w+-?\d+)|(\w+.\d+)::?ty1\se?gfp|(\w+.?\d*)::\S{1,5}gfp)'),
	 ('tag_flag','((\w+)\.flag|(\w+)::flag|(\w+)::\S*flag|(\w+)-flag|(\w+-\d+)_flag|(\w+)-\S{1,5}flag|flag-(\w+)|flag::(\w+)|flag\S*-(\w+)|flag\S*::(\w+)|(\w+)_flag|flag.tagged\s(\w+))'),
	 ('tag_myc','((\w+)\.myc|(\w+)::myc|(\w+)-myc|myc.tagged\s(\w+)|myc-(\w+)|myc::(\w+)|(\w+)::\S{1,5}myc|(\w+-\d+)_myc|(\w+)-\S{1,5}myc|myc\S*-(\w+)|myc\S{1,5}::(\w+)|(\w+)myc)'),
	 ('tag_PK','((\w+)::\S*pk|pk::(\w+)|(\w+)\s?pk|(\w{2,})-pk|(\w+)-\S*pk|pk\S*-(\w+)|v5-(\w+))'),
	 ('tag_tap','((\w+)\.tap|(\w+)::tap|(\w+)-tap|(\w+)-\w{2,}-tap|(\w+-\d+)_tap|tap::(\w+)|(\w+)::\w*tap|tap\S*-(\w+)|tap\S*::(\w+)|(\w+)\stap|(\w+)tap)'),
	 ('tag_T7','((\w+)::\S*t7|(\w+)-\S*t7|t7\S*-(\w+)|t7\S*::(\w+))') ])

#dictionnaire utilisant le mot-clé chip pour trouver la protéine-cible de l'essai
CHIP_DICO_old = OrderedDict ([
	('protein chip','((\w{2,})\sprotein\schip|(\w+\s\w+)\sprotein\schip)'),
	('BrdU IP','((\w{2,})\sbrdu\sip|(\w+)\s\w*\sbrdu\sip)'),
	('chip','(\w{2,})\s\s?chip'),
	('_chip','(\w{2,})_chip'),
	('chromatin immunoprecipitation', '(\w{2,})\s(?:wildtype|\w+)\snative\schromatin\simmunoprecipi?t?ation'),
	('chromatin ip against', '(chromatin\sip\sagainst\s(\w+\W\d+)|chromatin\sip\sagainst\s(\w{2,}))'),
	('IP','(\w{2,}_(\w+)_ip|(\w{2,})\s?ip|(\w+)_ip)'),
	('something IP','((\w{2,})_ip|(\w{2,})\s.*\sip)'),#vague
	('something chip','((\w{2,})\s\w*\schip|(\w+-\d+)\schip)') #peut être problématique car vague et parfois 2 options
	])

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
	 
#dictionnaire des les types cellulaires utilisés pour le ChIP-Seq avec celegans (not used here anymore)
CELL_TYPE_CELEGANS = OrderedDict ([
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
	('Young adult', '(ya|young|young\sadult)'),
	('Adult','(adult|old)'),
	('Total nuclei', 'total\snuclei'),
	('Muscle nuclei', 'muscle\snuclei'),
	('Oocytes', 'oocytes'),
	('DNA','dna'),
	('N/A', '.*')
	])

#dictionnaire des les types cellulaires utilisés pour le ChIP-Seq avec saccer et pombe (not used here anymore)
CELL_TYPE_SACCER = OrderedDict ([
	('N/A', '.*')
	])

