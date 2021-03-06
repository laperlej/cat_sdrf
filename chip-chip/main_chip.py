# -*- coding: utf-8 -*-
"""
input:
	file_names: list of sdrf.txt files
output:
	tab delimited file, contains first line of each sdrf file	
"""

import config_chip
from cat_sdrf.csv_manager import CsvManager
import sys
from antibody_filter_chip import filter_rows, assign_tag_multiple 

def main():
	if len(sys.argv) < 5:
		print "Usage: python main.py output_clean output_discard species file_names*"
		exit()

	output_clean= sys.argv[1]
	output_discard= sys.argv[2]
	species= sys.argv[3]
	file_names= sys.argv[4:]

	#species_name = config.species_dict[species]
	#gene_dict est assigné au gene_dict de l'espèce appellée (sys.argv[3]) dans la commande
	gene_dict = config_chip.GENE_DICT[species]
	gene_descrip_dict = config_chip.GENE_DESCRIP_DICT[species]

	fieldnames = config_chip.FIELDNAMES
	#section read_csv
	csv_manager = CsvManager(fieldnames.keys(), fieldnames)
	for file_name in file_names:
		csv_manager.read_csv(open(file_name, 'r'), preprocess=config_chip.PREPROCESS, preprocess2=config_chip.PREPROCESS2)

	#section filtre pour anticorps. Prend en argument le dictionnaire utilisé, les colonnes dans lesquelles on cherche, et la colone qu'on veux ensuite modifier
	csv_manager.rows = filter_rows(csv_manager.rows, config_chip.TARGET_DICO, ["4)assaytype","5)antibody", "6)target"], "clean_target")
	#section filtre pour le type d'essai. Prend en argument le dictionnaire utilisé, les colonnes dans lesquelles on cherche, et la colone qu'on veux ensuite modifier
	csv_manager.rows = filter_rows(csv_manager.rows, config_chip.ASSAY_DICO, ["4)assaytype", "5)antibody", "cell_type", "11)description", "13)other"], "clean_assay")
	#section filtre pour le type de cellules utilisées pour le ChIP-Seq. Prend en argument le dict utilisé, les col dans lesquelles on cherche, et la colone qu'on veux ensuite modifier
	csv_manager.rows = filter_rows(csv_manager.rows, config_chip.CELL_TYPE, ["cell_type"], "clean_celltype")
	#section filtre pour les tags et vides. Prend en argument les dictionnaires utilisés
	csv_manager.rows = assign_tag_multiple(csv_manager.rows, config_chip.TAG_DICO, config_chip.HISTONES_MARKS_DICO, gene_dict, gene_descrip_dict, config_chip.CHIP_DICO, config_chip.ANTIBODY_DICO)
	#section pour dupliquer la liste en un fichier clean et les lignes indésirables dans un autre fichier; prend en argument le dict d'espèces pour savoir laquelle utiliser
	csv_managers=csv_manager.split(config_chip.split_condition[species])	
	
	#section output
	csv_managers[0].write_csv(open(output_clean,'w'))
	csv_managers[1].write_csv(open(output_discard, 'w'))

if __name__ == "__main__":
	main()