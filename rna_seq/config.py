# -*- coding: utf-8 -*-

import os.path
import operator
import re
from collections import OrderedDict
import saccer, pombe, celegans
from antibody_filter import merge_cols

#iterate on each ine and return True the conditions are met.
def split_condition_aux(row, species):
	#Assays type to discard
	discard_assays=["rip-seq","chip-seq", "chip-chip", "unwanted", "mnase", "dnase", 'wgs', "atac", "brdu",'bisulfite-seq']
	#file types needed (updated for chip-chip)
	file_types = ['.sra', 'txt']
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
FIELDNAMES=OrderedDict([
	('1)identifier', ''),
	('1,1)Sample_title', ''),
	('2)filename', ''),
	('3)organism', ''),
	('4)clean_assay',''),
	('7)assaytype', ''),
	('10)treatment', ''),
	('11)Material_type',''),
	('13)cell_type',''),
	('14)strain', ''),
	('15)genotype', ''),
	('16)platform', ''),
	('Manufacturer', ''),
	('17)Sample_description', ''),
	('18)raw_files', ''),
	('19)all_supp_files', ''),
	('20)SRA_accessions', ''),
	('21)Experiment description', ''),
	('22)Protocol', ''),
	('23)Author(s)', ''),
	('Releasing group', ''),
	('24)Submission Date', ''),
	('25)Release Date', ''),
	('26)Pubmed ID', ''),
	('label', ''),
	('Selection', ''),
	('Other', '')])
#assay type dictionnary (WGS:Whole Genome Shotgun sequencing)
ASSAY_DICO = OrderedDict([
	('BrdU-ChIP', '(brdu\sip|brdu-ip)'),
	('ATAC-Seq', 'atac-seq'),
	('BrdU', 'brdu'),
	("WGS", 'wgs'),
	('FAIRE', 'faire'),
	("MNase-chip",'(mnase.treated|monococcal\snuclease|micrococcal\snuclease|chec\scleavage|chec\sexperiment|nucleosomal\sdna|micro-c)'),
	("ChIP-chip",'(chip|chromatin\simmunoprecipitation|immunoprecipitation\sof\snative\schromatin|genome\sbinding.occupancy\sprofiling\sby\sgenome\stiling\sarray)'),
	("ChIP-Seq", "chip-seq"),
	("DNase-Seq",'dnase'),
	("Bisulfite-Seq", "bisulfite"),
	("Rip-Seq", 'rip-seq'),
	("RNA-Seq", '(rna.seq|transcriptomic|total\srna|nascent\srna|polya\srna)'),
	('Non-genomic', '(non.genomic)'),
	('Okazaki fragment', 'okazaki\sfragment'),
	("other",'(other|genomic\sdna)'),
	("unwanted", '.*')
	])
