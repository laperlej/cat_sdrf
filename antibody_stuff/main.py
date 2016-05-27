"""
input:
	file_names: list of sdrf.txt files
output:
	tab delimited file, contains first line of each sdrf file	
"""

import config
from csv_manager import CsvManager
import sys
from antibody_filter import filter_rows

def main():
	file_names=sys.argv[1:]

	
	fieldnames = config.FIELDNAMES
	#section read_csv
	csv_manager = CsvManager(fieldnames.keys(), fieldnames)
	for file_name in file_names:
		csv_manager.read_csv(open(file_name, 'r'), preprocess=config.PREPROCESS)

	#section filtre pour anticorps
	csv_manager.rows = filter_rows(csv_manager.rows, config.TARGET_DICO, ["5)antibody", "6)target"], "clean_target")
	#section filtre pour le type d'essai
	csv_manager.rows = filter_rows(csv_manager.rows, config.ASSAY_DICO, ["4)assaytype"], "clean_assay")
	#section filtre pour les tags
	csv_manager.rows = tagfilter_rows(csv_manager.rows, config.TAG_DICO, ["8)strain", "9)genotype"], "clean_target")
		
	#section output
	csv_manager.write_csv(sys.stdout)

if __name__ == "__main__":
	main()