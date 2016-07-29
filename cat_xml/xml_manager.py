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
import xmltodict
from collections import OrderedDict

class XmlManager(object):
	def __init__(self, column_names, lambda_dict, sep=" | "):
		self.column_names = column_names
		self.lambda_dict = lambda_dict
		self.sep = sep
		#self.rows will be a list of OrderedDict (row) which will contain pairs of (row title:info)
		self.rows = []

	def read_xml(self, opened_file, sep=' | '):
		#Makes a list of orderedDict containing key:value pairs and lists and orderedDict
		mon_dict = xmltodict.parse(opened_file.read())
		self.contributor_dict = {}
		if 'Series' in mon_dict['MINiML']:
			self.series_dict = {'GSE':'', 'Title':'', 'Summary':'', 'Pubmed':''}
			self.series_to_gsm(mon_dict['MINiML']['Series'])
		if 'Contributor' in mon_dict['MINiML']:
			self.author_list(mon_dict['MINiML']['Contributor'])
		if 'Sample' not in mon_dict['MINiML']:
			pass
		else:	
			#If the series contains only one sample
			if type(mon_dict['MINiML']['Sample']) is not list:
				row = OrderedDict ([('1)identifier', ''), ('2)filename', ''),('3)organism', ''), ('4)clean_assay', ''), ('5)clean_target',''),('6)reliability', ''), ('7)assaytype', ''), ('8)antibody', ''), ('9)target', ''), ('10)treatment', ''),('11)Material_type', ''), ('12)clean_celltype',''), ('13)cell_type', ''), ('14)strain',''), ('15)genotype', ''), ('16)platform', ''), ('17)description', ''), ('18)raw_files', ''), ('19)all_supp_files', ''), ('20)SRA_files', ''), ('21)Experiment description', ''), ('22)Protocol', ''), ('23)Author(s)', ''), ('24)Submission Date', ''), ('25)Release Date', ''), ('26)Pubmed ID', ''), ('Other', '') ])
				self.id_list = []
				self.id_list.append(mon_dict['MINiML']['Sample']['@iid'])
				row['1)identifier'] = sep.join(self.id_list)
				self.rows.append(row)
			else:
				#Iteration on the list of all samples
				for x in range(len(mon_dict['MINiML']['Sample'])):
					row = OrderedDict ([('1)identifier', ''), ('2)filename', ''),('3)organism', ''), ('4)clean_assay', ''), ('5)clean_target',''),('6)reliability', ''), ('7)assaytype', ''), ('8)antibody', ''), ('9)target', ''), ('10)treatment', ''),('11)Material_type', ''), ('12)clean_celltype',''), ('13)cell_type', ''), ('14)strain',''), ('15)genotype', ''), ('16)platform', ''), ('17)description', ''), ('18)raw_files', ''), ('19)all_supp_files', ''), ('20)SRA_files', ''), ('21)Experiment description', ''), ('22)Protocol', ''), ('23)Author(s)', ''), ('24)Submission Date', ''), ('25)Release Date', ''), ('26)Pubmed ID', ''), ('Other', '') ])
					#Resets the lists after each sample
					self.id_list = []
					self.org_list = []
					self.treatment_list = []
					self.antibody_list = []
					self.target_list = []
					self.assay_list = []
					self.material_list = []
					self.gene_list = []
					self.strain_list = []
					self.descrip_list = []
					self.protocol_list = []
					self.supp_data = []
					self.other_list = []
					#since some GSM don't even have supplementary files
					if 'Supplementary-Data' not in mon_dict['MINiML']['Sample'][x]:
						self.id_list = []
						self.id_list.append(mon_dict['MINiML']['Sample'][x]['@iid'])
						row['1)identifier'] = sep.join(self.id_list)
						self.rows.append(row)
					else:	
						#Iteration on one Sample dictionnary
						for section in mon_dict['MINiML']['Sample'][x]:
							if section == '@iid':
								self.general_sample(self.id_list, mon_dict['MINiML']['Sample'][x]['@iid'])
							elif section == 'Title':
								self.descrip_sample(mon_dict['MINiML']['Sample'][x]['Title'])
							elif section == 'Organism':
								self.organism_sample(self.org_list, mon_dict['MINiML']['Sample'][x]['Organism'])
							elif section == 'Library-Source':
								self.general_sample(self.material_list, mon_dict['MINiML']['Sample'][x]['Library-Source'])
							elif section == 'Library-Strategy':
								self.general_sample(self.assay_list, mon_dict['MINiML']['Sample'][x]['Library-Strategy'])
							elif section == 'Library-Selection':
								self.general_sample(self.assay_list, mon_dict['MINiML']['Sample'][x]['Library-Selection'])
							elif section == 'Channel':
								if type(mon_dict['MINiML']['Sample'][x]['Channel']) is OrderedDict:
									#Iteration on the ordereddict that is Channel
									for key in mon_dict['MINiML']['Sample'][x]['Channel']:
										if 'Organism' in key:
											self.organism_sample(mon_dict['MINiML']['Sample'][x]['Channel']['Organism'])
										elif 'Source' in key:
											self.gene_list.append(mon_dict['MINiML']['Sample'][x]['Channel']['Source'])
										elif 'Molecule' in key:
											self.general_sample(self.material_list, mon_dict['MINiML']['Sample'][x]['Channel']['Molecule'])	
										elif 'Characteristics' in key:
											self.characteristics_sample(mon_dict['MINiML']['Sample'][x]['Channel']['Characteristics'])
										elif 'Treatment-Protocol' in key:
											self.general_sample(self.protocol_list, mon_dict['MINiML']['Sample'][x]['Channel']['Treatment-Protocol'])
										elif 'Growth-Protocol' in key:
											self.general_sample(self.protocol_list, mon_dict['MINiML']['Sample'][x]['Channel']['Growth-Protocol'])
										elif 'Extract-Protocol' in key:
											self.general_sample(self.protocol_list, mon_dict['MINiML']['Sample'][x]['Channel']['Extract-Protocol'])
										else:
											self.other_stuff_sample(mon_dict['MINiML']['Sample'][x]['Channel'][key])
								elif type(mon_dict['MINiML']['Sample'][x]['Channel']) is list:
									for num in range(len(mon_dict['MINiML']['Sample'][x]['Channel'])):
										for key in mon_dict['MINiML']['Sample'][x]['Channel'][num]:
											if 'Organism' in key:
												self.organism_sample(mon_dict['MINiML']['Sample'][x]['Channel'][num]['Organism'])
											elif 'Source' in key:
												self.gene_list.append(mon_dict['MINiML']['Sample'][x]['Channel'][num]['Source'])
											elif 'Molecule' in key:
												self.general_sample(self.material_list, mon_dict['MINiML']['Sample'][x]['Channel'][num]['Molecule'])	
											elif 'Characteristics' in key:
												self.characteristics_sample(mon_dict['MINiML']['Sample'][x]['Channel'][num]['Characteristics'])
											elif 'Treatment-Protocol' in key:
												self.general_sample(self.protocol_list, mon_dict['MINiML']['Sample'][x]['Channel'][num]['Treatment-Protocol'])
											elif 'Growth-Protocol' in key:
												self.general_sample(self.protocol_list, mon_dict['MINiML']['Sample'][x]['Channel'][num]['Growth-Protocol'])
											elif 'Extract-Protocol' in key:
												self.general_sample(self.protocol_list, mon_dict['MINiML']['Sample'][x]['Channel'][num]['Extract-Protocol'])
											else:
												self.other_stuff_sample(mon_dict['MINiML']['Sample'][x]['Channel'][num][key])		

							elif section == 'Characteristics':
								self.descrip_sample(mon_dict['MINiML']['Sample'][x]['Characteristics'])
							elif section == 'Instrument-Model':
								row['16)platform'] = mon_dict['MINiML']['Sample'][x]['Instrument-Model']['Predefined']
							elif section == 'Description':
								self.descrip_sample(mon_dict['MINiML']['Sample'][x]['Description'])
							elif section == 'Supplementary-Data':
								self.supp_data_sample(mon_dict['MINiML']['Sample'][x]['Supplementary-Data'])
							elif section == 'Contact-Ref':
								contact = mon_dict['MINiML']['Sample'][x]['Contact-Ref']['@ref']
								#row['23)Author(s)'] = mon_dict['MINiML']['Sample'][x]['Contact-Ref']['@ref']
							elif section == 'Status':
								for key in mon_dict['MINiML']['Sample'][x]['Status']:
									if 'Submission-Date' in key:
										row['24)Submission Date'] = mon_dict['MINiML']['Sample'][x]['Status']['Submission-Date']
									elif 'Release-Date' in key:
										row['25)Release Date'] = mon_dict['MINiML']['Sample'][x]['Status']['Release-Date']
									else:
										self.other_stuff_sample(mon_dict['MINiML']['Sample'][x]['Status'][key])
							else:
								self.other_stuff_sample(mon_dict['MINiML']['Sample'][x][section])

						row['1)identifier'] = sep.join(self.id_list)
						row['2)filename'] = self.series_dict['GSE']
						row['3)organism'] = sep.join(self.org_list)
						row['7)assaytype'] = sep.join(self.assay_list)
						row['8)antibody'] = sep.join(self.antibody_list)
						row['9)target'] = sep.join(self.target_list)
						row['10)treatment'] = sep.join(self.treatment_list)
						row['11)Material_type'] = sep.join(self.material_list)
						#row 12)clean_celltype and 13)cell_type not very useful now
						row['14)strain'] = sep.join(self.strain_list)
						row['15)genotype'] = sep.join(self.gene_list)
						#row['16)platform'] done earlier
						row['17)description'] = sep.join(self.descrip_list)
						row['19)all_supp_files'] = sep.join(self.supp_data)
						Exp_descrip = self.series_dict['Title'] + ' : ' + self.series_dict['Summary']
						#row['20)SRA_files'] not very useful now
						row['21)Experiment description'] = Exp_descrip
						row['22)Protocol'] = sep.join(self.protocol_list)
						row['23)Author(s)'] = self.contributor_dict[contact]
						#row['24)Submission Date'] and row['25)Release Date'] done earlier
						row['26)Pubmed ID'] = self.series_dict['Pubmed']
						row['Other'] = sep.join(self.other_list)
						self.rows.append(row)
	

	def general_sample(self, my_list, section):
		my_list.append(section)

	def organism_sample(self, section):
		if type(section) is list:
			for list_index in range(len(section)):
				if type(section[list_index]) is OrderedDict:
					self.org_list.append(section[list_index]['#text'])
				elif type(section[list_index]) is not OrderedDict:
					self.org_list.append(section[list_index])
		else:
			self.org_list.append(section['#text'])		
	def characteristics_sample(self, section):
		if type(section) is list:
			conditions = ['treatment', 'condition', 'growth', 'time', 'cycle', 'cell', 'temperature', 'fragmentation', 'synchronized']
			strain = ['strain', 'variant']
			gene = ['genetic', 'genotype', 'background']
			target = ['protein', 'epitope', 'target']
			#Iteration on the list, which is often composed of orderedDict ([('@tag', '...') , ('#text', '...')])
			for list_index in range(len(section)):
				#when some part of the list is just text
				if type(section[list_index]) is not OrderedDict:
					self.other_stuff_sample(section[list_index])
				#When the components of the list are OrderedDict
				elif any(condition in section[list_index]['@tag'] for condition in conditions):
					self.treatment_list.append(section[list_index]['#text'])
				elif any(item in section[list_index]['@tag'] for item in target):
					self.target_list.append(section[list_index]['#text'])	
				elif 'antibody' in section[list_index]['#text'] or 'antibody' in section[list_index]['@tag']:
					self.antibody_list.append(section[list_index]['#text'])
				elif any(item in section[list_index]['@tag'] for item in strain):
					self.strain_list.append(section[list_index]['#text'])
				elif any(item in section[list_index]['@tag'] for item in gene):	
					self.gene_list.append(section[list_index]['#text'])
				elif 'molecule' in section[list_index]['@tag']:	
					self.material_list.append(section[list_index]['#text'])
				elif 'library selection' in section[list_index]['@tag']:
					self.assay_list.append(section[list_index]['#text'])
				else:
					self.other_stuff_sample(section[list_index]['#text'])
		#probably a useless section 			
		elif type(section) is OrderedDict:
			for key in section.keys():
				if 'treatment' in section[key]:
					self.treatment_list.append(section['#text'])
				elif 'antibody' in key:
					self.antibody_list.append(section['#text'])
				elif 'strain' in key:
					self.strain_list.append(section['#text'])
				elif 'genotype' in key:	
					self.gene_list.append(section['#text'])
				else:
					self.other_stuff_sample(section)	
		else:
			if type(section) is OrderedDict:
				self.descrip_sample(section['#text'])
			else:
				if 'strain' in section:
					self.strain_list.append(section)
				else:	
					self.descrip_sample(section)

	def descrip_sample(self, section):
		if type(section) is OrderedDict:
			for key in section.keys():
				self.descrip_list.append(section['#text'])
		else:
			self.descrip_list.append(section)

	def supp_data_sample(self, section):
		if type(section) is list:
			for list_index in range(len(section)):
				self.supp_data.append(section[list_index]['#text'])
		elif type(section) is OrderedDict:
			self.supp_data.append(section['#text'])
		else:	
			self.supp_data.append(section)
	def other_stuff_sample(self, section):
		if type(section) is list:
			for list_index in range(len(section)):
				if type(section[list_index]) is list:
					for num in range(len(section[list_index])):
						self.other_list.append(section[list_index][item])
				elif type(section[list_index]) is OrderedDict:
					for key in section[list_index]:
						self.other_list.append(section[list_index][key])
				else:
					self.other_list.append(section[list_index])
		elif type(section) is OrderedDict:
			for key in section:
				if type(key) is OrderedDict:
					self.other_list.append(section[key])
		else:	
			self.other_list.append(section)
	
	def author_list(self, section):
		if type(section) is list:
			#Iteration on the list of contributors
			for num in range(len(section)):
				#number is associated to the contributor's number (ex: contrib1)
				number = section[num]['@iid']
				if 'Person' not in section[num]:
					names = []
					if 'Organization' in section[num]:
						if type(section[num]['Organization']) is list:
							org_name = " , ".join(item for item in section[num]['Organization'] if item is not None)
							#names.append(org_name)
						else:
							names.append(section[num]['Organization'])
				elif 'Person' in section[num]:
					names = []
					#Iteration on the First, middle and Last name of each contributor
					for person_name in section[num]['Person']:
						names.append(section[num]['Person'][person_name])		
				self.contributor_dict[number] = " ".join(names)	
							

	def series_to_gsm(self, section):
		if type(section) is list:
			for num in range(len(section)):
				if '@iid' in section[num]:
					self.series_dict['GSE'] = section[num]['@iid']
				elif 'Title' in section[num]:
					self.series_dict['Title'] = section[num]['Title']
				elif 'Pubmed' in section[num]:
					self.series_dict['Pubmed'] = section[num]['Pubmed-ID']
				elif 'Summary' in section[num]:
					self.series_dict['Summary'] = section[num]['Summary']
		elif type(section) is OrderedDict:
			for key in section:
				if '@iid' in key:
					self.series_dict['GSE'] = section['@iid']
				elif 'Title' in key:
					self.series_dict['Title'] = section[key]	
				elif 'Pubmed' in key:
					self.series_dict['Pubmed'] = section[key]
				elif 'Summary' in key:
					self.series_dict['Summary'] = section[key]									


	#Not needed with the xml files
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
#		print self.rows #in the shell, the special characters are written \u2103 and \xb5

	#not needed with xml files
	def empty_row(self):
		return {column_name:set() for column_name in self.column_names}
	#not needed with xml files
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
		#message d'erreur ici
		self.rows = [{k.encode('unicode-escape'):v.encode('unicode-escape') for k, v in row.iteritems()} for row in self.rows] 
		writer = csv.DictWriter(outfile,
			                    fieldnames=self.column_names,
			                    dialect='excel-tab')
		writer.writeheader()
		for row in self.rows:
			writer.writerow(row)

	
		
