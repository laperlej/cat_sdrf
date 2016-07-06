# -*- coding: utf-8 -*-
"""
input:
	file_names: list of sdrf.txt files
output:
	tab delimited file, contains first line of each sdrf file	
"""

import config
from csv_manager import CsvManager
import sys
from antibody_filter import bam_sam_filter_rows, filter_rows, assign_tag_multiple 

def main():
	if len(sys.argv) < 5:
		print "Usage: python main.py output_clean output_discard species file_names*"
		exit()

	output_clean= sys.argv[1]
	output_discard= sys.argv[2]
	species= sys.argv[3]
	file_names= sys.argv[4:]

	#species_name = config.species_dict[species]
	#gene_dict is assigned to the species's gene_dict
	gene_dict = config.GENE_DICT[species]
	gene_descrip_dict = config.GENE_DESCRIP_DICT[species]
	CELL_TYPE = config.CELL_TYPE_DICT[species]

	fieldnames = config.FIELDNAMES
	# read_csv section
	csv_manager = CsvManager(fieldnames.keys(), fieldnames)
	for file_name in file_names:
		# Preprocess get information from idf files (see config)
		csv_manager.read_csv(open(file_name, 'r'), preprocess_1=config.PREPROCESS_1, preprocess_2=config.PREPROCESS_2, preprocess_3=config.PREPROCESS_3, preprocess_4=config.PREPROCESS_4, preprocess_5=config.PREPROCESS_5, preprocess_6=config.PREPROCESS_6, preprocess_7=config.PREPROCESS_7 )

	#Filter section for antibody. Take in argument the dictionnary, the columns in which we search and the column to modify
	csv_manager.rows = filter_rows(csv_manager.rows, config.TARGET_DICO, ["4)assaytype","5)antibody", "6)target"], "clean_target")
	#Filter section for the assay type. Take in argument the dictionnary, the columns in which we search and the column to modify
	csv_manager.rows = filter_rows(csv_manager.rows, config.ASSAY_DICO, ["4)assaytype", "Material_type", 'cell_type', "11)description"], "clean_assay")
	#Filter section for the cell type. Take in argument the dictionnary, the columns in which we search and the column to modify
	csv_manager.rows = filter_rows(csv_manager.rows, CELL_TYPE, ["cell_type"], "clean_celltype")
	#Filter section for the Bam-Sam files. Take in argument the dictionnary, the columns in which we search and the column to modify
	csv_manager.rows = bam_sam_filter_rows(csv_manager.rows, config.FILETYPES ,["13)other", "Protocol"], "BAM-SAM")
	#section filtre pour les tags et vides. Prend en argument les dictionnaires utilisés
	csv_manager.rows = assign_tag_multiple(csv_manager.rows, config.TAG_DICO, config.HISTONES_MARKS_DICO, gene_dict, gene_descrip_dict, config.CHIP_DICO, config.ANTIBODY_DICO)
	
	csv_manager.fix_dup_gsm(['1)identifier', '2)filename', '3)organism', 'clean_assay', 'clean_target', 'reliability', '4)assaytype', '5)antibody', '6)target', '7)treatment', '8)strain', '9)genotype', '10)platform', '11)description'])
	#section pour dupliquer la liste en un fichier clean et les lignes indésirables dans un autre fichier; prend en argument le dict d'espèces pour savoir laquelle utiliser
	csv_managers=csv_manager.split(config.split_condition[species])	
	
	#section output
	csv_managers[0].write_csv(open(output_clean,'w'))
	csv_managers[1].write_csv(open(output_discard, 'w'))

if __name__ == "__main__":
	main()