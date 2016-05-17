"""
input:
	file_names: list of sdrf.txt files
output:
	tab delimited file, contains first line of each sdrf file	
"""

import config
from csv_manager import CsvManager
import sys

def main():
	file_names=sys.argv[1:-1]
	out_name=sys.argv[-1]
	
	fieldnames = config.FIELDNAMES

	csv_manager = CsvManager(fieldnames)
	for file_name in file_names:
		csv_manager.read_csv(open(file_name, 'r'), preprocess=config.PREPROCESS)
	csv_manager.write_csv(open(out_name, 'w'))

if __name__ == "__main__":
	main()