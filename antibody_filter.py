# -*- coding: utf-8 -*-
"""
analyser la colonne 5)antibody et 6)taget pour ajouter la cible ou le tag à la colonne clean_target
analyser la colonne 8)strain, 9)genotype, 11)description pour déterminer la protéine-cible avec les flag (regex +match avec dictionnaire de genes)
"""

import re
from collections import OrderedDict

# Itère sur chaque ligne du dictionnaire rows (contient les info des fichiers sdrf)
def filter_rows(rows, target_dico, input_cols, output_col):
	for row in rows:
		row = filter_row(row, target_dico, input_cols, output_col) 
	return rows

#Joint le contenu des colonnes (input_cols) sans séparateur
def merge_cols(row, input_cols):
	return "|".join(row[input_col].lower() for input_col in input_cols)

def filter_row(row, target_dico, input_cols, output_col):
	"""
	itère sur les regex de target-dico et les compare avec l'information dans input_col jusqu'à ce qu'il y ait un match.
	Fonction multi-task, peut être utilisée aussi pour filtrer les essais 
	input: 
		row: dictionnaire, la clé est le titre de la colonne et la valeur est le contenu de la colonne
		target_dico: le dictionnaire de regex pour les cibles
		input_cols: liste des colonnes dans lesquelles on cherche (concaténées).
		output_col: colonne changée lorsqu'il y a un match
	output:
		row: colonne output_col = clé du dictionnaire (info) s'il y a eu un match
	"""
	
	#section qui itère sur target_dico	
	new_value = ""
	for info in target_dico.keys():
		searchtarget = merge_cols(row, input_cols)
		if re.search(target_dico[info], searchtarget):
			new_value = info
			break
	 
	row[output_col] = new_value
	return row


def assign_tag_multiple(rows, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico):
	"""
	input: 
		tag_dico: dictionnaire de regex pour matcher la cible des tag (..... - tag)
		histone_dico: dictionnaire des marques d'histones avec regex
		gene_dico: dictionnaire de gènes avec regex
		gene_descrip_dico: dictionnaire de noms standards et des alias pour les gènes (non redondant)
		chip_dico: dictionaire de regex pour matcher la cible des chip (.... chip)
	"""
	#fonction qui permet d'utiliser tous les cpu disponibles (accélère)
	import multiprocessing
	multiprocessing.cpu_count()
	num_cpu = multiprocessing.cpu_count()
	try:	
		p = multiprocessing.Pool(num_cpu)
		arguments = ((row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico) for row in rows)
		results = p.map(assign_tag_parallel, arguments)
		for row, result in zip(rows, results):		
			row["clean_target"], row["reliability"] = result
		return rows	 		
	finally:
		p.close()
		p.join()
	exit(1)


def assign_tag_parallel(data):
	return assign_tag(*data)


#section qui trie entre les input, mock, tag et types d'essai
#appelle compare_tag, compare_chip et compare_directly (différents niveaux de comparaisons) 
def assign_tag(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico):
	#Assigne 'N/A' à la colonne clean_target si l'essai est mnase ou dnase
	if "MNase" in row["clean_assay"] or "DNase" in row["clean_assay"] or "FAIRE" in row["13)other"] or "FAIRE" in row["11)description"]:
		return "N/A", "assay type (1)"
	#Assigne 'input' à la colonne clean_target si le mot-clé input est trouvé
	elif "input" in merge_cols(row, ["4)assaytype", "11)description","13)other"]).lower():
		return "input", "keyword (2)"
	#Assigne 'mock' à la colonne clean_target si le mot-clé 'mock' est trouvé
	elif "mock" in merge_cols(row, ["4)assaytype", "11)description","13)other"]).lower():
		return "Mock", "keyword (2)"
	
	elif "empty" in row["clean_target"]:
		#trouve les mock dans les essais qui n'ont pas été identifiées comme portant un tag
		if "notag" in merge_cols(row, ["13)other"]).lower():
			return "Mock", "keyword (2)"
		else:
			#appelle compare_chip; 
			return compare_chip(row, histones_dico, gene_dico, gene_descrip_dico, chip_dico)
			
	# section qui appelle compare_tag pour les lignes qui contiennent 'tag' et qui ne sont pas des 'Mock'
	elif "tag" in row["clean_target"]:
	#	si untagged dans les col 4, 11 et 13
		if "untagged" in merge_cols(row, ["4)assaytype", "11)description","13)other"]).lower() or "no tag" in merge_cols(row, ["4)assaytype", "11)description","13)other"]).lower():
	#		ajouter Mock à clean target:
			return "Mock","keyword (1)"
		else:	
	#		si tagged, retourne fonction compare_tag	
			return compare_tag(row,tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico)	
	else:
		return row["clean_target"], "target_dico (1)" 

def compare_tag(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico): 
	#ajouter filtre avec histones_dico ici
	tagged = row["clean_target"]
	#compare le regex du tag de clean_target à ce qu'il y a dans les colonnes 4-6-8-9-11
	#idéalement, enlever ici la colonne strain, amène des faux positifs
	match = re.search(tag_dico[tagged],merge_cols(row,["4)assaytype", "6)target", "8)strain", "9)genotype", "11)description"]).lower())
	if match:
		#compare match du tag avec histones_dico
		for hist in histones_dico.keys():
			if re.search(histones_dico[hist], match.group(1)):	
				return hist, "tag to histone (2)"
		#compare match du tag avec gene_dico
		for gene in gene_dico.keys():
			if re.search(gene_dico[gene], match.group(1)):
				return gene, "tag to gene (3)"
		
		#compare match du tag avec gene_descrip_dico
		for gene in gene_descrip_dico.keys():
			if re.search(gene_descrip_dico[gene], match.group(1)):	
				return gene, "tag to gene descr (4)"
	
	return compare_chip(row, histones_dico, gene_dico, gene_descrip_dico, chip_dico)

#cherche le mot-clé "ChIP" (entre autres) et compare ce qui vient avec chip avec les dictionnaires d'histones de gènes et d'alias
def compare_chip(row, histones_dico, gene_dico, gene_descrip_dico, chip_dico): 
	#itère sur les regex du chip_dico (pour trouver la cible des ChIP)
	for regex in chip_dico.keys():
		match = re.search(chip_dico[regex], row["11)description"].lower())
		if match:
			for hist in histones_dico.keys():
				#compare match du regex chip avec le dico d'histones
				if re.search(histones_dico[hist], match.group(1)):
					return hist, "chip to histone (3)"
				
			for gene in gene_dico.keys():
				#compare match du regex avec gene_dico
				if re.search(gene_dico[gene], match.group(1)):
					return gene, "chip to gene (3)"	
				
			for gene in gene_descrip_dico.keys():
				#compare match du regex avec gene_descrip_dico
				if re.search(gene_descrip_dico[gene], match.group(1)):
					return gene, "chip to gene descr (4)"
	return compare_directly(row, histones_dico, gene_dico, gene_descrip_dico)		

#compare le contenu des colonnes (contenu plus spécific au plus large) avec les dictionnaires d'histones, de gènes et d'alias
def compare_directly(row, histones_dico, gene_dico, gene_descrip_dico):
	#itère sur le dictionnaire d'histones 
	for hist in histones_dico.keys():
		if re.search(histones_dico[hist], merge_cols(row, ["1)identifier", "5)antibody", "6)target", "11)description"]).lower()):
			return hist, "histone mark (3)"
		elif re.search(histones_dico[hist], row["9)genotype"].lower()):
			return hist, "histone mark (4)"
	#itère sur le dict de gènes et compare à plusieurs niveaux (moins au plus large en principe)
	for gene in gene_dico.keys():
		if re.search(gene_dico[gene], merge_cols(row, ["5)antibody","6)target", "11)description"]).lower()):
			return gene, "gene (4)"
		elif re.search(gene_dico[gene], row["8)strain"].lower()):
			return gene, "gene (5)"
		elif re.search(gene_dico[gene], row["9)genotype"].lower()):
			return gene, "gene (5)"

	#itère sur le dict d'alias et compare sur plusieurs niveaux (moins au plus large en principe)
	for gene in gene_descrip_dico.keys():
		if re.search(gene_descrip_dico[gene], merge_cols(row, ["5)antibody","6)target", "11)description"]).lower()):
			return gene, "gene descr (5)"
		elif re.search(gene_descrip_dico[gene], row["8)strain"].lower()):
			return gene, "gene descr (5)"
		elif re.search(gene_descrip_dico[gene], row["9)genotype"].lower()):
			return gene, "gene descr (5)"
	
	for hist in histones_dico.keys():
		elif re.search(histones_dico[hist], row["13)other"].lower()):
			return hist, "histone mark (4)"
	#itère sur le dictionnaire de gènes et compare sur les lignes les moins spécifiques
	for gene in gene_dico.keys():
		if re.search(gene_dico[gene], merge_cols(row, ["1)identifier", "13)other"]).lower()):
			return gene, "gene (6)"
	#itère sur le dictionnaire d'alias et compare sur les lignes les moins spécifiques
	for gene in gene_descrip_dico.keys():
		if re.search(gene_descrip_dico[gene], merge_cols(row, ["1)identifier", "13)other"]).lower()):
			return gene, "gene descr (6)"		


	return row["clean_target"], "target_dico (1)"	






	"""
		elif empty in row["clean_target]:
		#section qui itère sur histones_dico pour trouver un match sur deux niveaux, moins au plus large
		for hist in histones_dico.keys():
			if re.search(histones_dico[hist], merge_cols(row, ["1)identifier", "5)antibody", "6)target", "11)description"]).lower()):
				return hist, "histone mark (1)"
			elif re.search(histones_dico[hist], merge_cols(row, ["9)genotype", "13)other"]).lower()):
				return hist, "histone mark (2)"		
		#section qui itère sur les dictionnaire de gènes, d'abord le dict de gènes, puis le dict de gènes contenant les alias
		for gen in gene_dico.keys():
			if re.search(gene_dico[gen], merge_cols(row, ["5)antibody","6)target"]).lower()):
				return gen, "gene (1)"
			elif re.search(gene_dico[gen], row["11)description"].lower()):
				return gen, "gene (2-3)"
		for gen_descr in gene_descrip_dico.keys():
			if re.search(gene_descrip_dico[gen_descr], merge_cols(row, ["5)antibody","6)target"]).lower()):
				return gen_descr, "gene descr (3)"
			elif re.search(gene_descrip_dico[gen_descr], row["11)description"].lower()):
				return gen_descr, "gene descr (4)"		
		#section qui itère sur le dictionnaire de gènes, à deux niveaux (moins au plus large)
		for gen in gene_dico.keys():
			if re.search(gene_dico[gen], row["9)genotype"].lower()):
				return gen, "gene large (4)"
			elif re.search(gene_dico[gen], merge_cols(row, ["1)identifier", "13)other"]).lower()):
				return gen, "gene large (5)"	

		#section qui itère sur le dictionnaire de nom standard et description de gènes, à deux niveaux (moins au plus large)
		for gen_descr in gene_descrip_dico.keys():
			if re.search(gene_descrip_dico[gen_descr], row["9)genotype"].lower()):
				return gen_descr, "gene descr large (4)"
			elif re.search(gene_descrip_dico[gen_descr], merge_cols(row, ["1)identifier", "13)other"]).lower()):
				return gen_descr, "gene descr large (5)"
		return row["clean_target"], "target_dico (1)"
	"""	
		
