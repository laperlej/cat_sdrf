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
from antibody_filter import filter_rows, assign_tag_multiple 

def main():
	file_names=sys.argv[1:]
	#remplacer par quelque chose comme
	#output_file= sys.argv[1]
	#discard_file= sys.argv[2]
	#input_files= sys.argv[3:]

	
	fieldnames = config.FIELDNAMES
	#section read_csv
	csv_manager = CsvManager(fieldnames.keys(), fieldnames)
	for file_name in file_names:
		csv_manager.read_csv(open(file_name, 'r'), preprocess=config.PREPROCESS)

	#ajouter fonction discard_line (csv_manager, où on aurait aussi la commande pour écrire dans discard) 
	#
	#section filtre pour anticorps
	csv_manager.rows = filter_rows(csv_manager.rows, config.GENE_scerevisiae, config.TARGET_DICO, ["5)antibody", "6)target"], "clean_target")
	#section filtre pour le type d'essai
	csv_manager.rows = filter_rows(csv_manager.rows, config.ASSAY_DICO, ["4)assaytype"], "clean_assay")
	#section filtre pour les tags
	csv_manager.rows = assign_tag_multiple(csv_manager.rows, config.TAG_DICO, config.GENE_scerevisiae)
		
	#section output
	csv_manager.write_csv(sys.stdout)
	#remplacer par 
	#csv_manager.write_csv(output_file)

if __name__ == "__main__":
	main()