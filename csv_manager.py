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
import copy
import StringIO
import json
import itertools
from xml.etree import ElementTree

class CsvManager(object):
	def __init__(self, column_names, lambda_dict, sep=" | "):
		self.column_names = column_names
		self.lambda_dict = lambda_dict
		self.sep = sep
		self.rows = []

	def read_csv(self, csvfile, preprocess_1=lambda csvfile, row: None, preprocess_2=lambda csvfile, row: None, preprocess_3=lambda csvfile, row: None, preprocess_4=lambda csvfile, row: None, preprocess_5=lambda csvfile, row: None, preprocess_6=lambda csvfile, row: None, preprocess_7=lambda csvfile, row: None ):
		csvcontent = csvfile.read()
		csvcontent = self.fixduplicates(csvcontent)
		#StringIO permits to read a string as a file
		reader = csv.DictReader(StringIO.StringIO(csvcontent), dialect='excel-tab')
		#reader = csv.DictReader(csvfile, dialect='excel-tab')
		for row in reader:
			#preprocess gets information from idf files (see config)
			preprocess_1(csvfile, row)
			preprocess_2(csvfile, row)
			preprocess_3(csvfile, row)
			preprocess_4(csvfile, row)
			preprocess_5(csvfile, row)
			preprocess_6(csvfile, row)
			preprocess_7(csvfile, row)
			self.rows.append(self.translate_row(row))

	def fixduplicates(self,csvcontent):
		"""Checks all the headers and numerates the headers that have the same name"""
		#Splits the first line from the rest (separated by a newline) and then in a list of string (tab-separated)
		headers = csvcontent.split('\n')[0].split("\t")
		#Splits the other lines (which are newline-separated)
		content = csvcontent.split('\n')[1:]
		new_headers = []
		count = 1
		#Iterate on the headers
		for title in headers:
			if title not in new_headers:
				new_headers.append(title)
			else:
				#Modifies the header if it is a duplicate (adds a number equivalent to the count)
				new_title = title + str(count)
				new_headers.append(new_title)
				count += 1
		#Joins the headers (tab-separated) and then joins the content of the file (newline-separated) and returns the result	
		return "\n".join(["\t".join(new_headers)]+ content)   

	def fix_dup_gsm(self, uniq_titles):
		""" Reunites lines that are identical for the information in the columns listed (uniq titles) and concatenate their informations if it is not the same"""
		uniq_lines = {}
		for row in self.rows:
			#if the key (information) is in uniq_titles, then it is added  to the dictionnary uniq_cols
			uniq_cols = {key:value for key, value in row.iteritems() if key in uniq_titles}
			non_uniq_cols = {key:value for key, value in row.iteritems() if key not in uniq_titles}
			#json.dumps returns a string for uniq_cols
			string_dict = json.dumps(uniq_cols, sort_keys=True, ensure_ascii=False)
			#if string_dict is not in uniq_lines
			if uniq_lines.get(string_dict, False):
				for key in non_uniq_cols.iterkeys():
					cell1 = uniq_lines[string_dict][key]
					cell2 = non_uniq_cols[key]
					#splits the content of cell 1 and 2 
					cell1 =  cell1.split(' | ')
					cell1 = [x for x in cell1 if x not in ["", " "]]
					cell2 =  cell2.split(' | ')
					cell2 = [x for x in cell2 if x not in ["", " "]]
					#Joins the content of cell 1 and 2 without repetition of content
					new_cell = set(cell1)|set(cell2)
					non_uniq_cols[key]= ' | '.join(new_cell)
				uniq_lines[string_dict]= non_uniq_cols	
			else:
				uniq_lines[string_dict]= non_uniq_cols
		self.rows = [dict(itertools.chain({key.encode('utf-8'):value.encode('utf-8') for key, value in json.loads(key).iteritems()}.iteritems(), value.iteritems())) for key,value in uniq_lines.iteritems()]


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

	#Create a copy of rows; result[0] gets the lines returning True to the condition, the other lines go into result[1]
	def split(self,condition=lambda row:True):
		result=[self, copy.deepcopy(self)]
		result[1].rows=[row for row in self.rows if not condition(row)]
		result[0].rows=[row for row in self.rows if condition(row)]
		return result

	def write_csv(self, outfile):
		writer = csv.DictWriter(outfile,
			                    fieldnames=self.column_names,
			                    dialect='excel-tab')
		writer.writeheader()
		for row in self.rows:
			writer.writerow(row)

	
		
