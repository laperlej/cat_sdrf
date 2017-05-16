# -*- coding: utf-8 -*-
"""
input:
	file_names: list of xml files
output:
	tab delimited file, contains lines from each xml file	
"""
import codecs
import config
from xml_manager import XmlManager
import sys
from antibody_filter import filter_rows, assign_tag_multiple, condition_rows 

def main():
	if len(sys.argv) < 5:
		print ("Usage: python main.py output_clean output_discard species file_names*")
		exit()

	output_clean= sys.argv[1]
	output_discard= sys.argv[2]
	species= sys.argv[3]
	file_names= sys.argv[4:]

	#species_name = config.species_dict[species]
	#gene_dict is assigned to the species's gene_dict
	gene_dict = config.GENE_DICT[species]
	gene_descrip_dict = config.GENE_DESCRIP_DICT[species]
	CELL_TYPE = config.CELL_TYPE_DICT[species]

	fieldnames = config.FIELDNAMES
	# assings to a variable the name of the class and its arguments
	xml_manager = XmlManager(fieldnames.keys(), fieldnames)
	# read_xml section	
	for file_name in file_names:
		# Preprocesses get information from idf files (see config)
		xml_manager.read_xml(codecs.open(file_name, 'r', encoding='utf-8'))
	xml_manager.fix_dup_gsm(['1)identifier'])
	
	# Filter section for antibody. Take in argument the dictionnary, the columns in which we search and the column to modify
	xml_manager.rows = filter_rows(xml_manager.rows, config.TARGET2, config.HISTONES_MARKS_DICO, ["7)assaytype","8)antibody", "9)target"], "5)clean_target")
	# Filter section for the assay type. Take in argument the dictionnary, the columns in which we search and the column to modify
	xml_manager.rows = filter_rows(xml_manager.rows, config.ASSAY_DICO, {}, ["7)assaytype", "1,1)Sample_title", "11)Material_type", '13)cell_type', "17)Sample_description"], "4)clean_assay")
	# Filter section for tags and 'empty' lines in 'clean_target'. Takes in argument all the dictionnaries used
	xml_manager.rows = assign_tag_multiple(xml_manager.rows, config.TAG_DICO, config.HISTONES_MARKS_DICO, gene_dict, gene_descrip_dict, config.CHIP_DICO, config.ANTIBODY_DICO)
	#Fills a column 'Selection' which explains why the sample is in the 'clean file' or in the 'discard file'
	xml_manager.rows = condition_rows(xml_manager.rows, species)
	# Duplicate the list of dictionnaries (rows) in a 'clean' file and puts the unwanted lines in a 'discard' file; the conditions depends on the species
	xml_managers=xml_manager.split(config.split_condition[species])	
	
	# Output section
	xml_managers[0].write_csv(codecs.open(output_clean,'w', encoding='utf-8'))
	xml_managers[1].write_csv(codecs.open(output_discard,'w', encoding='utf-8')) 
	

if __name__ == "__main__":
	main()
