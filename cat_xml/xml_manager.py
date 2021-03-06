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
		#mon_dict = xmltodict.parse(opened_file.read(), encoding='utf-8')
		mon_dict = xmltodict.parse(opened_file.read())
		self.contributor_dict = {}
		special_characters = {'∆': 'Delta', 'ɛ': 'Epsilon', 'δ': 'Delta', 'α':'Alpha'}
		if 'Series' in mon_dict['MINiML']:
			# information left from GSE section of file: 'Status database', 'Submission-Date', 'Release-Date', 'Last-Update-Date', 'Accession database', 'Type', 'Contributor-Ref', 'Sample-Ref', 'Contact-Ref', 'Supplementary-Data' and 'Relation-Type'
			#Information left from GPL section of file: 'Platform iid', 'Status database', 'Submission-Date', 'Release-Date', 'Last-Update-Date', 'Title', 'Accession database', 'Technology', 'Distribution', 'Organism', 'Description', 'Manufacturer', 'Manufacture-Protocol'
			self.series_dict = {'GSE':'', 'Title':'', 'Summary':'', 'Pubmed':'' , 'Overall-Design':''}
			self.series_to_gsm(mon_dict['MINiML']['Series'])
		if 'Contributor' in mon_dict['MINiML']:
			self.author_list(mon_dict['MINiML']['Contributor'])
		if 'Sample' not in mon_dict['MINiML']:
			pass
		else:	
			#If the series contains only one sample
			if type(mon_dict['MINiML']['Sample']) is not list:
				row = OrderedDict ([('1)identifier', ''), ('1,1)Sample_title', ''), ('2)filename', ''),('3)organism', ''), ('4)clean_assay', ''), ('5)clean_target',''),('6)reliability', ''), ('7)assaytype', ''), ('8)antibody', ''), ('9)target', ''), ('10)treatment', ''),('11)Material_type', ''), ('12)clean_celltype',''), ('13)cell_type', ''), ('14)strain',''), ('15)genotype', ''), ('16)platform', ''), ('17)Sample_description', ''), ('18)raw_files', ''), ('19)all_supp_files', ''), ('20)SRA_accessions', ''), ('21)Experiment description', ''), ('22)Protocol', ''), ('23)Author(s)', ''), ('24)Submission Date', ''), ('25)Release Date', ''), ('26)Pubmed ID', ''), ('Other', '') ])
				self.id_list = []
				self.id_list.append(mon_dict['MINiML']['Sample']['@iid'])
				row['1)identifier'] = sep.join(self.id_list)
				self.rows.append(row)
			else:
				#Iteration on the list of all samples
				for x in range(len(mon_dict['MINiML']['Sample'])):
					row = OrderedDict ([('1)identifier', ''), ('1,1)Sample_title', ''), ('2)filename', ''),('3)organism', ''), ('4)clean_assay', ''), ('5)clean_target',''),('6)reliability', ''), ('7)assaytype', ''), ('8)antibody', ''), ('9)target', ''), ('10)treatment', ''),('11)Material_type', ''), ('12)clean_celltype',''), ('13)cell_type', ''), ('14)strain',''), ('15)genotype', ''), ('16)platform', ''), ('17)Sample_description', ''), ('18)raw_files', ''), ('19)all_supp_files', ''), ('20)SRA_accessions', ''), ('21)Experiment description', ''), ('22)Protocol', ''), ('23)Author(s)', ''), ('24)Submission Date', ''), ('25)Release Date', ''), ('26)Pubmed ID', ''), ('Other', '') ])
					#Resets the lists after each sample
					self.id_list = []
					self.org_list = []
					self.treatment_list = []
					self.antibody_list = []
					self.target_list = []
					self.assay_list = []
					self.material_list = []
					self.cell_list = []
					self.gene_list = []
					self.strain_list = []
					self.platform_list = []
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
							all_protocols = ['Growth-Protocol', 'Extract-Protocol', 'Treatment-Protocol', 'Label-Protocol', 'Scan-Protocol', 'Hybridization-Protocol']
							if section == '@iid':
								self.general_sample(self.id_list, mon_dict['MINiML']['Sample'][x]['@iid'])
							elif section == 'Title':
								row['1,1)Sample_title'] = mon_dict['MINiML']['Sample'][x]['Title']
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
											self.descrip_list.append(mon_dict['MINiML']['Sample'][x]['Channel']['Source'])
										elif 'Molecule' in key:
											self.general_sample(self.material_list, mon_dict['MINiML']['Sample'][x]['Channel']['Molecule'])	
										elif 'Characteristics' in key:
											self.characteristics_sample(mon_dict['MINiML']['Sample'][x]['Channel']['Characteristics'])
										elif any(protocol in key for protocol in all_protocols):
											self.general_sample(self.protocol_list, mon_dict['MINiML']['Sample'][x]['Channel'][key])
										else:
											self.characteristics_sample(mon_dict['MINiML']['Sample'][x]['Channel'][key])
								elif type(mon_dict['MINiML']['Sample'][x]['Channel']) is list:
									for num in range(len(mon_dict['MINiML']['Sample'][x]['Channel'])):
										for key in mon_dict['MINiML']['Sample'][x]['Channel'][num]:
											if 'Organism' in key:
												self.organism_sample(mon_dict['MINiML']['Sample'][x]['Channel'][num]['Organism'])
											elif 'Source' in key:
												self.descrip_list.append(mon_dict['MINiML']['Sample'][x]['Channel'][num]['Source'])
											elif 'Molecule' in key:
												self.general_sample(self.material_list, mon_dict['MINiML']['Sample'][x]['Channel'][num][key])	
											elif 'Characteristics' in key:
												self.characteristics_sample(mon_dict['MINiML']['Sample'][x]['Channel'][num]['Characteristics'])
											elif any(protocol in key for protocol in all_protocols):
												self.general_sample(self.protocol_list, mon_dict['MINiML']['Sample'][x]['Channel'][num][key])
											else:
												self.characteristics_sample(mon_dict['MINiML']['Sample'][x]['Channel'][num][key])		

							elif any(protocol in section for protocol in all_protocols):
								self.general_sample(self.protocol_list, mon_dict['MINiML']['Sample'][x][section])
							elif section == 'Characteristics':
								#anything goes here?
								self.descrip_sample(mon_dict['MINiML']['Sample'][x]['Characteristics'])
							elif 'Platform-Ref' in section:
								self.platform_list.append(mon_dict['MINiML']['Sample'][x]['Platform-Ref']['@ref'])
							elif section == 'Instrument-Model':
								# Used tag: 'Instrument-Model'
								self.platform_list.append(mon_dict['MINiML']['Sample'][x]['Instrument-Model']['Predefined'])
							elif section == 'Description':
								#verify what goes here
								self.descrip_sample(mon_dict['MINiML']['Sample'][x]['Description'])
							elif section == 'Data-Processing':
								self.general_sample(self.protocol_list, mon_dict['MINiML']['Sample'][x]['Data-Processing'])	
							elif section == 'Supplementary-Data':
								self.supp_data_sample(mon_dict['MINiML']['Sample'][x]['Supplementary-Data'])
							# Assign to the variable contact the contributor number
							elif section == 'Contact-Ref':
								contact = mon_dict['MINiML']['Sample'][x]['Contact-Ref']['@ref']	
							elif section == 'Status':
								for key in mon_dict['MINiML']['Sample'][x]['Status']:
									if 'Submission-Date' in key:
										row['24)Submission Date'] = mon_dict['MINiML']['Sample'][x]['Status']['Submission-Date']
									elif 'Release-Date' in key:
										row['25)Release Date'] = mon_dict['MINiML']['Sample'][x]['Status']['Release-Date']
									elif 'Last-Update-Date' in key:
										self.other_list.append(mon_dict['MINiML']['Sample'][x]['Status']['Last-Update-Date'])
									else:
										self.other_stuff_sample(mon_dict['MINiML']['Sample'][x]['Status'][key])
							else:
								self.other_stuff_sample(mon_dict['MINiML']['Sample'][x][section])

						# Used tag: 'Sample iid' in sample part of file
						row['1)identifier'] = sep.join(self.id_list)
						# Used tag for row['1,1)Sample_title'] : 'Title' in the sample part of file 
						#Used tag: 'Series iid' in the GSE part of the file
						row['2)filename'] = self.series_dict['GSE']
						#Used tag: 'Organism' from 'Channel' section of sample part of file
						row['3)organism'] = sep.join(self.org_list)
						#Used tag: 'Library-Strategy', 'Library-Selection'
						row['7)assaytype'] = sep.join(self.assay_list)
						#Used tag:
						row['8)antibody'] = sep.join(self.antibody_list)
						#Used tag:
						row['9)target'] = sep.join(self.target_list)
						#Used tag: 'Growth-Protocol', 'Treatment-Protocol' and 'Extract-Protocol' (from 'Channel','Characteristics' section) and 'Data-Processing' 
						row['10)treatment'] = sep.join(self.treatment_list)
						#Used tag: 'Library-Source', 'Molecule' in 'Channel' section
						row['11)Material_type'] = sep.join(self.material_list)
						row ['12)clean_celltype'] = sep.join(self.cell_list)
						#13)cell_type not very useful now
						#Used tag:
						row['14)strain'] = sep.join(self.strain_list)
						#Used tag:
						row['15)genotype'] = sep.join(self.gene_list)
						# Used tag : 'Platform-Ref' and 'Instrument-Model'
						row['16)platform'] = sep.join(self.platform_list)
						#Used tag: 'Description'
						row['17)Sample_description'] = sep.join(self.descrip_list)
						#Used tag:
						row['19)all_supp_files'] = sep.join(self.supp_data)
						#row['20)SRA_files'] not very useful now
						Exp_descrip = self.series_dict['Title'] + ' | ' + self.series_dict['Summary'] + ' | ' + self.series_dict['Overall-Design']
						#Used tag: concatenation of 'Title', 'Summary' and 'Overall-Design' from the GSE part
						row['21)Experiment description'] = Exp_descrip
						row['22)Protocol'] = sep.join(self.protocol_list)
						#Consist of the name associated to a contributor number mentionned in the GSM part; the contributor number and name are taken from a list or contributors described in the GSE part
						row['23)Author(s)'] = self.contributor_dict[contact]
						#row['24)Submission Date'] and row['25)Release Date'] done earlier
						#Used tag: Pubmed ID in the GSE part
						row['26)Pubmed ID'] = self.series_dict['Pubmed']
						#Used tag:
						row['Other'] = sep.join(self.other_list)
						#replace the special characters (ɛ, δ, α, ∆)
						for key in special_characters:
							#iteration on the dictionnary row
							for section in row:
								row[section] = row[section].replace(key,special_characters[key])
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
		key_value = ''
		cell_type = ['cell type', 'cell line']
		#carefull that 'tag' is not specific but it is used as a tag
		target = ['protein', 'epitope', 'target', 'tag', 'flag', 'ChIP', 'h2b', 'histone', 'IP against', 'target of ip', 'tagged protein']
		conditions = ['treatment', 'condition', 'growth', 'time', 'timing', 'cycle', 'cell', 'temperature', 'fragmentation', 'synchronized', 'media', 'medium', 'buffer', 'culture', 'stage', 'status', 'carbon', 'glucose', 'selection', 'plasmid', 'vector', 'drug', 'dmso', 'stress', 'concentration', 'mnase', 'agent', 'mononucleosome', 'spike-in', 'enzyme', 'ploid', 'environnement', 'treated', 'ymc', 'digested with', 'digestion', 'addition', 'transformation', 'depleted factor', 'sucrose',  'sex', 'h2o2', 'hours at 37','triton', 'immunodepletion', 'knock', 'equivalents of ercc spike', 'transformed with','break induction', 'rna purification', 'fluorescence','transfection', 'facs-sorted population', 'construct', 'transposon','resistance', 'transcription', 'factor', 'fluor', 'cyanine dye', 'cy3', 'cy5', 'sirna', 'rna deletion', 'rnai deletion', 'crispri guiderna', 'mixed percentage', 'incubation', 'harvest', 'passage']
		#removed : 'rna', 'dna'
		material = ['molecule', 'tissue', 'organelle', 'cell part', 'mrna type', 'shrna', 'rna subtype', 'material', 'genomic dna', 'nucleosomal DNA', 'Input', 'input', 'chromatin']
		strain = ['strain', 'Strain', 'background', 'variant', 'mutant', 'yeast', 'parents', 'wild type', 'MAT-a', 'direct rna sequence from']
		gene = ['genetic', 'genotype', 'allele', 'phenotype', 'gene deletion', 'rnai deletion', 'modification', 'bearing', 'genome', 'variation']
		
		junk = ['hotspot',  'batch', 'repetition', 'replicate', 'repeat', 'experiment', 'isolate number', 'index pair', 'grna libraries',  'matched wild type sample', 'barcode', 'sample identifier', 'tandem repeat', 'chd1-ume6 fusion', 'primer', 'index', 'strategy',  'sort', 'capture', 'lentivirally',  'sequencing chip', 'crosslink',  'cmc use', 'ID', 'fragment size', 'application', 'paired-end',  'vendor', 'oligonucleotide', 'processed data', 'Biotin', 'biotin', 'Sample', 'SAMPLE', 'replication', 'Affymetrix', 'sequenced with', 'matched wild type sample', 'average']
		if type(section) is list:
			#Iteration on the list, which is often composed of orderedDict ([('@tag', '...') , ('#text', '...')])
			for list_index in range(len(section)):
				#when some part of the list is just text
				if type(section[list_index]) is not OrderedDict:
					self.other_stuff_sample(section[list_index])
				#When the components of the list are OrderedDict 
				#Used tags: 'cell type', 'cell line', 'tissue' (catches 'tissue/cell line')
				elif any(item in section[list_index]['@tag'].lower() for item in cell_type):
					print (section[list_index]['@tag'])
					print (section[list_index]['#text'])
					self.cell_list.append(section[list_index]['#text'])
				# Used tag: 'protocol', but does not contain much info; catches 'protocol', 'growth protocol' ,'treatment protocol', 'culture protocol', 'harvest method' and also 'growt protocol';
				elif 'protocol' in section[list_index]['@tag'].lower() or 'harvest method' in section[list_index]['@tag'].lower():	
					self.protocol_list.append(section[list_index]['#text'])	
				#Used tag here catches 'Treatment', 'culture condition', 'growth condition' and 'growth protocol'; valid info
				elif 'mg/l' in section[list_index]['#text'] or 'uM' in section[list_index]['#text']:
					key_value = section[list_index]['@tag'] + ' : ' +  section[list_index]['#text']
					self.treatment_list.append(section[list_index]['#text'])	
				#Some info; valid; catches tag 'ip'
				elif 'IP against' in section[list_index]['#text']:
					self.target_list.append(section[list_index]['#text'])		
				# Used tag 'antibody' catches 'antibody', 'chip-antibody', 'chip antibody', 'chip antibody lot #', 'chip antibody cat. #', 'chip antibody vendor', 'chip antibody reference'
				elif 'antibody' in section[list_index]['#text'] or 'antibody' in section[list_index]['@tag']:
					self.antibody_list.append(section[list_index]['#text'])
				#Some info; valid
				elif 'catalog' in section[list_index]['@tag']:
					self.antibody_list.append(section[list_index]['#text'])	
				# Used tag: 'sample type'; Maybe should go in Material; some lines are descriptive
				elif 'sample type' in section[list_index]['@tag']:
					self.descrip_list.append(section[list_index]['#text'])
				#Used tag: 'genotype', '' ; lots of info; valid
				elif any(item in section[list_index]['@tag'] for item in gene):
					self.gene_list.append(section[list_index]['#text'])
				#Lots of info; valid
				elif any(item in section[list_index]['@tag'].lower() for item in strain):
					self.strain_list.append(section[list_index]['#text'])
				# Used tag : 'assay' catches 'library selection', 'assayed molecue' and 'assay'
				elif 'assay' in section[list_index]['@tag'] or 'library selection' in section[list_index]['@tag']:
					self.assay_list.append(section[list_index]['#text'])
				# Used tag: 'library' catches 'library strategy' and 'library type'
				elif 'library' in section[list_index]['@tag'].lower():
					self.assay_list.append(section[list_index]['#text'])
				#Lots of info; valid; many used tags: 'material type', 'molecule subtype', 'tissue', 'cell part', 'organelle', 'rna subtype', 'mrna type', 'stable expression of shrna'
				elif any(item in section[list_index]['@tag'] for item in material):
					#nothing here
					if any(item in section[list_index]['@tag'] for item in junk):
						key_value =  section[list_index]['@tag'] + ' : ' +  section[list_index]['#text']
						self.other_stuff_sample(key_value)
					#nothing here
					elif 'rna subset' in section[list_index]['@tag'].lower():
						self.descrip_list.append(section[list_index]['#text'])
					else:
						self.material_list.append(section[list_index]['#text'])
				#Lots of info here; to many used tags for treatment
				elif any(condition in section[list_index]['@tag'].lower() for condition in conditions):
					# Assign to the string key_value both the key and the value (like 'Mnase concentration : 10 mM')
					key_value =  section[list_index]['@tag'] + ' : ' +  section[list_index]['#text']
					self.treatment_list.append(key_value)
				# Valid info, but for 'tag': 'MATa ade2-1 can1-100 HIS3 leu2-3,112 trp1-1 ura3-1 RAD5+ ISW1-FL3-KanMX snf2-delta::URA3'; 'tag' not specific
				elif any(item in section[list_index]['@tag'] for item in target):
					key_value =  section[list_index]['@tag'] + ' : ' +  section[list_index]['#text']
					self.target_list.append(key_value)		
				elif any(item in section[list_index]['@tag'] for item in junk):
					key_value =  section[list_index]['@tag'] + ' : ' +  section[list_index]['#text']
					self.other_stuff_sample(key_value)
				#Not specific tags, but valid info
				elif 'od' in section[list_index]['@tag'] or 'age' in section[list_index]['@tag']:
					key_value =  section[list_index]['@tag'] + ' : ' +  section[list_index]['#text']
					self.treatment_list.append(key_value)
				elif 'rna' in section[list_index]['@tag'].lower():
					self.material_list.append(section[list_index]['#text'])
				else:
					# The leftover goes in the 'Other' section
					key_value =  section[list_index]['@tag'] + ' : ' +  section[list_index]['#text']
					self.other_stuff_sample(key_value)
				
		elif type(section) is OrderedDict:
			for key in section.keys():
				#Some info; valid
				if 'antibody' in section['@tag']:
					self.antibody_list.append(section['#text'])
				# Some info here; valid
				elif 'stage' in section['@tag']:
					key_value = section['@tag']+ ': ' + section['#text']
					self.treatment_list.append(key_value)
				#Nothing here
				elif any(item in section['@tag'] for item in target):
					self.target_list.append(section['#text'])	
		#		elif section['@tag']:	
		#			self.descrip_list.append(section['#text'])
				# Lots of info; valid
				elif any(item in section['@tag']for item in gene):
					self.gene_list.append(section['#text'])
				#Info going here; valid
				elif any(item in section['@tag'] for item in strain):
					#Not sure if should go in genotype, strain or description
					if 'Strain Background is ' in section['@tag']:
						key_value = section['@tag'] + ': ' + section['#text']
						self.strain_list.append(key_value)
					else:
						#lots of info goes here
						self.strain_list.append(section['#text'])
				# Some info; valid
				elif any(item in section['@tag'] for item in cell_type):
					self.cell_list.append(section['#text'])
				#lots of info; valid
				elif any(condition in section['@tag'] for condition in conditions):
					key_value = section['@tag']+ ': ' + section['#text']
					self.treatment_list.append(key_value)
				#Some info here
				elif any(item in section['@tag'] for item in material):	
					self.material_list.append(section['#text'])
				# Some info goes here; valid
				elif 'protocol' in section['@tag'] or 'harvest' in section['@tag']:	
					self.protocol_list.append(section['#text'])	
				#Not specific tags
				elif 'od' in section['@tag'] or 'age' in section['@tag']:
					key_value =  section['@tag'] + ' : ' +  section['#text']
					self.treatment_list.append(key_value)
				else:
					# The leftover goes in the 'Other' section
					key_value = section['@tag']+ ': ' + section['#text']
					self.other_list.append(key_value)
	
		else:
			#nothing here
			if 'library strategy' in section:
				self.protocol_list.append(section)
			#Some info 
			elif any(condition in section for condition in conditions):
				self.treatment_list.append(section)
			#nothing here
			elif 'antibody' in section:
				self.antibody_list.append(section)
			#nothing here
			elif any(item in section for item in target):
				# catch publication if it contains 'tag' or 'flag'
				if 'et al.' in section:
					self.other_list.append(section)
				else:	
					self.target_list.append(section)	
			#noting here
			elif any(item in section for item in strain):
				self.strain_list.append(section)
			#noting here
			elif any(item in section for item in gene):
				self.gene_list.append(section)
			#Nothing here
			elif any(item in section for item in material):	
				self.material_list.append(section)
			#nothing here
			elif 'BrdU' in section or 'brdu' in section:
				self.assay_list.append(section)	
			else:	
				# The leftover (there is some info) goes in the 'Other' section
				self.other_list.append(section)
				#self.other_stuff_sample(section)

	def descrip_sample(self, section):
		if type(section) is list:
			for list_index in range(len(section)):
				if type(section[list_index]) is list:
					for num in range(len(section[list_index])):
						self.descrip_list.append(section[list_index][num])
				elif type(section[list_index]) is OrderedDict:
					for key in section[list_index]:
						self.descrip_list.append(section[list_index]['#text'])
				else:
					self.descrip_list.append(section[list_index])
		elif type(section) is OrderedDict:
			for key in section.keys():
				if type(key) is OrderedDict:
					for key2 in section[key]:
						self.descrip_list.append(section[key]['#text'])
				else:
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
				if type(section[key]) is OrderedDict:
					for key2 in section[key]:
						self.other_list.append(section[key][key2])
				elif type(section[key]) is list:
					for num in range(len(section[key])):
						if type(section[key][num]) is OrderedDict:
							for key2 in section[key][num]:
								self.other_list.append(section[key][num][key2])
						elif type(section[key]) is not OrderedDict:
							self.other_list.append(section[key][num])		
				else:
					self.other_list.append(section[key])
		else:	
			self.other_list.append(section)
	
	def author_list(self, section):
		if type(section) is list:
			#Iteration on the list of contributors at the begining of the GSE file
			for num in range(len(section)):
				#number is associated to the contributor's number (ex: contrib1)
				number = section[num]['@iid']
				if 'Person' not in section[num]:
					names = []
					if 'Organization' in section[num]:
						if type(section[num]['Organization']) is list:
							org_name = " , ".join(item for item in section[num]['Organization'] if item is not None)
							names.append(org_name)
						else:
							names.append(section[num]['Organization'])	
				elif 'Person' in section[num]:
					names = []
					#Iteration on the First, middle and Last name of each contributor
					for person_name in section[num]['Person']:
						#names is a list of all the names for one contributor (first name, middle name and last name)
						names.append(section[num]['Person'][person_name])		
				# Assign a complete name as the value to the key that is the contributor number (ex contrib1 : John Doe)
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
					if type(section[key]) is list:
						for i in range(len(section[key])):
							Ids = []
							Ids.append(section[key][i])
						self.series_dict['Pubmed'] = "|".join(Ids)
					else:
						self.series_dict['Pubmed'] = section[key]
				elif 'Summary' in key:
					self.series_dict['Summary'] = section[key]
				elif 'Overall-Design' in key:
					self.series_dict['Overall-Design'] = section[key]									


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
			uniq_cols = {key:value for key, value in row.items() if key in uniq_titles}
			non_uniq_cols = {key:value for key, value in row.items() if key not in uniq_titles}
			#json.dumps returns a string for uniq_cols
			string_dict = json.dumps(uniq_cols, sort_keys=True, ensure_ascii=False)
			#if string_dict is not in uniq_lines
			if uniq_lines.get(string_dict, False):
				for key in non_uniq_cols:
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
		self.rows = [dict(itertools.chain({key.encode('utf-8'):value.encode('utf-8') for key, value in json.loads(key).items()}.items(), value.items())) for key,value in uniq_lines.items()]
#		print (self.rows) #in the shell, the special characters are written \u2103 and \xb5

	#not needed with xml files
	def empty_row(self):
		return {column_name:set() for column_name in self.column_names}
	#not needed with xml files
	def translate_row(self, row):
		new_row = self.empty_row()
		norm_row = utils.norm_keys(row)
		for title in norm_row:
			for column_name in self.column_names:
				if self.lambda_dict[column_name](title, norm_row):
					info = str(norm_row.get(title, " ")).strip()
					if info not in ["", "None"]:
						new_row[column_name].add(info)
					break
		return {key:self.sep.join(content) for key, content in new_row.items()}

	#Create a copy of rows; result[0] gets the lines returning True to the condition, the other lines go into result[1]
	def split(self,condition=lambda row:True):
		result=[self, copy.deepcopy(self)]
		result[1].rows=[row for row in self.rows if not condition(row)]
		result[0].rows=[row for row in self.rows if condition(row)]
		return result

	def write_csv(self, outfile):
		tmp_rows = []
		for row in self.rows:
			tmp_row = {}
			for k, v in row.items():
				if isinstance(k, bytes):
					tmp_row[k.decode("utf-8")]=v.decode("utf-8")
				else:
					tmp_row[k]=v
			tmp_rows.append(tmp_row)
		self.rows = tmp_rows
		#self.rows = [{k.decode("utf-8"):v.decode("utf-8") if isinstance(k, bytes) else k:v for k, v in row.items()} for row in self.rows]
		#self.rows = [{k.encode('utf-8'):v.encode('utf-8') for k, v in row.items()} for row in self.rows] 
		#self.rows = [{k.decode('utf-8'):v for k, v in row.items()} for row in self.rows] 
		#self.rows = [{k:v for k, v in row.items()} for row in self.rows] 
		writer = csv.DictWriter(outfile,
			                    fieldnames=self.column_names,
			                    dialect='excel-tab')
		writer.writeheader()
		for row in self.rows:
			writer.writerow(row)

	
		
