"""
input:
	file_names: list of sdrf.txt files
output:
	tab delimited file, contains first line of each sdrf file	
"""

import csv
import os.path
import utils

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
		for key, value in self.column_dict.iteritems():
			content = set()
			for title in value:
				info = norm_row.get(title, "").strip()
				if info:
					content.add(info)
			new_row[key] = self.sep.join(content)
		return new_row

	def write_csv(self, outfile);
		writer = csv.DictWriter(outfile,
			                    fieldnames=self.column_dict.keys(),
			                    dialect='excel-tab')
		for row in self.rows:
			writer.writerow(row)
