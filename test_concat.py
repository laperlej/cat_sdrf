# -*- coding: utf-8 -*-
from collections import OrderedDict
#Histone marks and polymerase dictionnary
HISTONES_MARKS_DICO = OrderedDict([
	('RNAPII_ser5P', 's5\sphosphorylated\srna\spolii'),
	("POLIII",'(pol3$|pol3\s|pol3-|pol\s?iii)'),
	("POLII",'(pol2|pol\sɛ|pol\s?ii)'),
	("POLI",'(pol1|pol\sα|pol\s?i)'),
	("POL31",'(pol31|pol\sδ)'),
	('H3K4me1', '(h3k4me1|monomethylated\sh3k4|h3\s.?mono\smethyl\sK4|ab8895)'),
	('H3K4me2', '(h3k4me2|dimethylated\sh3k4|h3\s.?di\smethyl\sK4)'),
	('H3K4me3', '(h3k4me3|trimethylated\sh3k4|h3\s.?tri\smethyl\sK4|ab8580|39159|305-34819|h3.?k4.?me3|ab8678)'),
	('H3K14ac', '(h3k14ac|07-353)'),
	('H3K14', '(h3k14)'),
	('H3K36me3','(h3k36me3|ab9050|300-95289)'),
	('H3K36me2','(h3k36me2)'),
	('H3K36me','(h3k36me)'),
	('H3K36','(h3k36)'),
	('H3K56ac','(h3k56ac|07-677)'),
	('H3K56','(h3k56)'),
	('H3K79me3', 'h3k79me3'),
	('H3K79', 'h3k79'),
	('H3K9ac','(h3k9ac|06-942)'),
	('H3K9me3','(h3k9me3|h3.?k9\s?me3|ab8898)'),
	('H3K9me2','(h3k9me2|h3.?k9\s?me2)'),
	('H3K9','h3k9'),
	('H3K4me3', '(h3.?k4.?me3)'),
	('H3K4','h3k4'),
	('H3K27me3', 'h3k27me3'),
	('H3R2me2', 'h3r2me2'),
	('H3K9-14ac','(ach3\sk9,14|h3k9-14ac)'),
	('H3','(h3|ab12079|05-928|07-690|ab1791)'),
	('H4K16ac','(h4k16ac|07-329)'),
	('H4K16','(h4k16)'),
	('H4K12ac', 'h4k12ac'),
	('H4K12', 'h4k12'),
	('H4K44ac','h4k44ac'),
	('H4K44','h4k44'),
	('H4K5ac','(h4k5ac|07-327)'),
	('H4K5','(h4k5)'),
	('H4K20me1','(h4k20me1|ab9051)'),
	('H4ac','(h4ac|39177)'),
	('H4','(h4|ab7311)'),
	('H2A.Z','(htz1|h2a\.?z|pht1|ab4626)'),
	('H2A', '(h2a|ab13923)'),
	('H2B', '(htb1|htb2|h2b|ab1790)')])

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
	('empty', '.*')
	])


all_targets = OrderedDict ([])
all_targets.update(HISTONES_MARKS_DICO)
#this dict comes last since it ends with a regex catching anything (when not using an empty dict)
all_targets.update(TARGET2)
for key in all_targets.iterkeys():
	print key, all_targets[key]


