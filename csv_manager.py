"""
input:
	file_names: list of sdrf.txt files
output:
	tab delimited file, contains first line of each sdrf file	
"""

import csv
import os.path
import utils
import re

class CsvManager(object):
	def __init__(self, column_dict, sep=" | "):
		self.column_dict = column_dict
		self.sep = sep
		self.rows = []

	def read_csv(self, csvfile, preprocess=lambda csvfile, row: None):
		reader = csv.DictReader(csvfile, dialect='excel-tab')
		for row in reader:
			preprocess(csvfile, row)
			self.rows.append(self.translate_row(row))

	def translate_row(self, row):
		norm_row = utils.norm_keys(row)
		new_row = {}
		for key1, regex in self.column_dict.iteritems():
			content = set()
			for key2 in norm_row.keys():
				if re.search(regex, key2):
					info = norm_row.pop(key2, "").strip()
					if info:
						content.add(info)
			new_row[key1] = self.sep.join(content)
		return new_row

	def write_csv(self, outfile):
		writer = csv.DictWriter(outfile,
			                    fieldnames=self.column_dict.keys(),
			                    dialect='excel-tab')
		writer.writeheader()
		for row in self.rows:
			writer.writerow(row)
