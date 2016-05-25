# -*- coding: utf-8 -*-

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
	def __init__(self, column_names, lambda_dict, sep=" | "):
		self.column_names = column_names
		self.lambda_dict = lambda_dict
		self.sep = sep
		self.rows = []

	def read_csv(self, csvfile, preprocess=lambda csvfile, row: None):
		reader = csv.DictReader(csvfile, dialect='excel-tab')
		for row in reader:
			preprocess(csvfile, row)
			self.rows.append(self.translate_row(row))

	def empty_row(self):
		return {column_name:set() for column_name in self.column_names}

	def translate_row(self, row):
		new_row = self.empty_row()
		norm_row = utils.norm_keys(row)
		for title in norm_row.iterkeys():
			for column_name in self.column_names:
				if self.lambda_dict[column_name](title, norm_row):
					info = str(norm_row.get(title, " ")).strip()
					if info not in ["", "None"]:
						new_row[column_name].add(info)
					break
		return {key:self.sep.join(content) for key, content in new_row.iteritems()}

	

	def write_csv(self, outfile):
		writer = csv.DictWriter(outfile,
			                    fieldnames=self.column_names,
			                    dialect='excel-tab')
		writer.writeheader()
		for row in self.rows:
			writer.writerow(row)

	
		
