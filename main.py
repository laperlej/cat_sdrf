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

	#section filtre
	csv_manager.rows = filter_rows(csv_manager.rows, config.TARGET_DICO, "5)antibody", "clean_target")
		
	#section output
	csv_manager.write_csv(sys.stdout)

if __name__ == "__main__":
	main()