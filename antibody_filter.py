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


def assign_tag_multiple(rows, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico, antibody_dico):
	"""
	input: 
		tag_dico: dictionnaire de regex pour matcher la cible des tag (..... - tag)
		histone_dico: dictionnaire des marques d'histones avec regex
		gene_dico: dictionnaire de gènes avec regex
		gene_descrip_dico: dictionnaire de noms standards et des alias pour les gènes (non redondant)
		chip_dico: dictionaire de regex pour matcher la cible des chip (.... chip)
		antibody_dico: dictionnaire de numéros d'anticorps
	"""
	#fonction qui permet d'utiliser tous les cpu disponibles (accélère)
	import multiprocessing
	multiprocessing.cpu_count()
	num_cpu = multiprocessing.cpu_count()
	try:	
		p = multiprocessing.Pool(num_cpu)
		arguments = ((row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico, antibody_dico) for row in rows)
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
def assign_tag(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico, antibody_dico):
	#Assigne 'N/A' à la colonne clean_target si l'essai est mnase, dnase ou FAIRE-Seq
	if "MNase" in row["clean_assay"] or "DNase" in row["clean_assay"]:
		return "N/A", "assay type (1)"
	elif "faire" in merge_cols(row,["13)other", "11)description"]):
		return "N/A", "assay type (1)"	
	#Assigne 'input' à la colonne clean_target si le mot-clé input est trouvé dans la ligne
	elif "input" in merge_cols(row, ["4)assaytype", "11)description","13)other"]).lower():
		return "input", "keyword (2)"
	#Assigne 'mock' à la colonne clean_target si le mot-clé 'mock' est trouvé dans la ligne
	elif "mock" in merge_cols(row, ["4)assaytype", "11)description","13)other"]).lower():
		return "Mock", "keyword (2)"
	#Assigne 'mock' à la colonne clean_target si le mot-clé 'non antibody control' est trouvé
	elif "non antibody control" in merge_cols(row, ["4)assaytype", "11)description","13)other"]).lower():
		return "Mock", "keyword (1)"

	elif "empty" in row["clean_target"]:
		#trouve les mock dans les essais qui n'ont pas été identifiées comme portant un tag
		if "notag" in merge_cols(row, ["13)other"]).lower():
			return "Mock", "keyword (2)"
		else:
			#appelle différentes fonctions pour trouver une cible
			var_search_target = search_target(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico)
			if var_search_target is not None:
				return var_search_target
			var_search_antibody = search_antibody(row, antibody_dico)
			if var_search_antibody is not None:
				return var_search_antibody	
			var_compare_tag2 = 	compare_tag2(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico)
			if var_compare_tag2 is not None:
				return var_compare_tag2
			var_compare_chip1 = compare_chip1(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico)
			if var_compare_chip1 is not None:
				return var_compare_chip1
			var_compare_tag_larger2 = 	compare_tag_larger2(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico)
			if var_compare_tag_larger2 is not None:
				return var_compare_tag_larger2
			var_compare_chip2 = compare_chip2(row, histones_dico, gene_dico, gene_descrip_dico, chip_dico)
			if var_compare_chip2 is not None:
				return var_compare_chip2
			var_compare_directly =	compare_directly(row, histones_dico, gene_dico, gene_descrip_dico)
			if var_compare_directly is not None:
				return var_compare_directly	
				
	# section qui appelle compare_tag pour les lignes qui contiennent 'tag' et qui ne sont pas des 'Mock'
	elif "tag" in row["clean_target"]:
	#	si untagged dans les col 4, 11 et 13
		if "untagged" in merge_cols(row, ["4)assaytype", "11)description","13)other"]).lower() or "no tag" in merge_cols(row, ["4)assaytype", "11)description","13)other"]).lower():
	#		ajouter Mock à clean target:
			return "Mock","keyword (1)"
		else:	
			#appelle différentes fonctions pour trouver une cible
			var_search_target = search_target(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico)
			if var_search_target is not None:
				return var_search_target
			var_search_antibody = search_antibody(row, antibody_dico)
			if var_search_antibody is not None:
				return var_search_antibody	
			var_compare_tag1 = 	compare_tag1(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico)
			if var_compare_tag1 is not None:
				return var_compare_tag1
			var_compare_chip1 = compare_chip1(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico)
			if var_compare_chip1 is not None:
				return var_compare_chip1
			var_compare_tag_larger1 = compare_tag_larger1(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico)
			if var_compare_tag_larger1 is not None:
				return var_compare_tag_larger1
			var_compare_chip2 = compare_chip2(row, histones_dico, gene_dico, gene_descrip_dico, chip_dico)
			if var_compare_chip2 is not None:
				return var_compare_chip2
			var_compare_directly =	compare_directly(row, histones_dico, gene_dico, gene_descrip_dico)
			if var_compare_directly is not None:
				return var_compare_directly	
	else:
		return row["clean_target"], "target_dico (1)"


def search_target(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico):
	#cherche la cible dans les col 4-5-6 en comparant directement dans le dict de gènes
	for gene in gene_dico.keys():
		if re.search(gene_dico[gene],merge_cols(row,["4)assaytype", "5)antibody", "6)target"])):
			return gene, "target (2)"
	return None	

def search_antibody(row, antibody_dico):
	#cherche pour des numéros d'anticorps dans la col 11)description
	for antibody in antibody_dico.keys():
		if re.search(antibody_dico[antibody],merge_cols(row,["11)description"])):
			return antibody, "antibody no (2)"		

def compare_tag1(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico): 
	#pour "tag" dans clean_target
	#cherche la cible d'un tag dans les colonnes 4-9 (plus spécifique normalement)
	tagged = row["clean_target"]
	#compare le regex (valeur) dont le tag de clean_target est la clé à ce qu'il y a dans les colonnes 4-6-11 (contiennent des occurences de cible-tag)
	match = re.search(tag_dico[tagged],merge_cols(row,["4)assaytype", "11)description"]))
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
	return None

#cherche le mot-clé "ChIP" (entre autres) et compare le match (ce qui vient avec chip) avec les dictionnaires d'histones, de gènes et d'alias
def compare_chip1(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico): 
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
	return None

def compare_tag_larger1(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico): 
	#pour "tag" dans clean_target
	#cherche la cible d'un tag dans les colonnes 8-9 (plus large)
	tagged = row["clean_target"]
	#compare le regex du tag de clean_target à ce qu'il y a dans les colonnes 8-9
	#la colonne 8)strain et 9)genotype amènent parfois des faux positifs (beaucoup de gènes)
	match = re.search(tag_dico[tagged],merge_cols(row,["cell_type", "8)strain", "9)genotype"]).lower())
	if match:
		#compare match du tag avec histones_dico
		for hist in histones_dico.keys():
			if re.search(histones_dico[hist], match.group(1)):	
				return hist, "tag to histone (3)"
		#compare match du tag avec gene_dico
		for gene in gene_dico.keys():
			if re.search(gene_dico[gene], match.group(1)):
				return gene, "tag to gene (4)"
		
		#compare match du tag avec gene_descrip_dico
		for gene in gene_descrip_dico.keys():
			if re.search(gene_descrip_dico[gene], match.group(1)):	
				return gene, "tag to gene descr (4)"
	return None


def compare_tag2(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico):
	#pour "empty" dans clean_target
	#cherche pour des tags qui n'ont pas été indiqués dans les colonnes 5-6; recherche assez spécifique
	for tag in tag_dico:
		#compare le regex (valeur) dont le tag de clean_target est la clé à ce qu'il y a dans les colonnes 4-6-11
		match = re.search(tag_dico[tag],merge_cols(row,["cell_type", "11)description"]).lower())
		if match:
			#compare match du tag avec histones_dico
			for hist in histones_dico.keys():
				if re.search(histones_dico[hist], match.group(1)):	
					return hist, "search tag to histone (2)"
			#compare match du tag avec gene_dico
			for gene in gene_dico.keys():
				if re.search(gene_dico[gene], match.group(1)):
					return gene, "search tag to gene (4)"
			#compare match du tag avec gene_descrip_dico
			for gene in gene_descrip_dico.keys():
				if re.search(gene_descrip_dico[gene], match.group(1)):	
					return gene, "search tag to gene descr (4)"
	return None

def compare_tag_larger2(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico):
	#pour "empty" dans clean_target
	#cherche pour des tags qui n'ont pas été indiqués dans les colonnes 5-6; recherche plus large
	for tag in tag_dico:
		#compare le regex (valeur) dont le tag de clean_target est la clé à ce qu'il y a dans les colonnes 8-9
		match = re.search(tag_dico[tag],merge_cols(row,["cell_type", "8)strain", "9)genotype"]).lower())
		if match:
			#compare match du tag avec histones_dico
			for hist in histones_dico.keys():
				if re.search(histones_dico[hist], match.group(1)):	
					return hist, "search tag to histone (2)"
			#compare match du tag avec gene_dico
			for gene in gene_dico.keys():
				if re.search(gene_dico[gene], match.group(1)):
					return gene, "search tag to gene (4)"
			#compare match du tag avec gene_descrip_dico
			for gene in gene_descrip_dico.keys():
				if re.search(gene_descrip_dico[gene], match.group(1)):	
					return gene, "search tag to gene descr (4)"				
	return None		

#cherche le mot-clé "ChIP" (entre autres) et compare ce qui vient avec chip avec les dictionnaires d'histones de gènes et d'alias
def compare_chip2(row, histones_dico, gene_dico, gene_descrip_dico, chip_dico): 
	#itère sur les regex du chip_dico (pour trouver la cible des ChIP)
	for regex in chip_dico.keys():
		match = re.search(chip_dico[regex], merge_cols(row,["1)identifier", "cell_type", "8)strain", "9)genotype"]).lower())
		if match:
			for hist in histones_dico.keys():
				#compare match du regex chip avec le dico d'histones
				if re.search(histones_dico[hist], match.group(1)):
					return hist, "chip to histone (4)"
				
			for gene in gene_dico.keys():
				#compare match du regex avec gene_dico
				if re.search(gene_dico[gene], match.group(1)):
					return gene, "chip to gene (4)"	
				
			for gene in gene_descrip_dico.keys():
				#compare match du regex avec gene_descrip_dico
				if re.search(gene_descrip_dico[gene], match.group(1)):
					return gene, "chip to gene descr (4)"
	return None			

#compare le contenu des colonnes (contenu plus spécific au plus large) avec les dictionnaires d'histones, de gènes et d'alias
def compare_directly(row, histones_dico, gene_dico, gene_descrip_dico):
	#Première passe, plus spécifique
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
	#itère sur le dict d'alias et compare sur plusieurs niveaux (moins au plus large en principe)
	for gene in gene_descrip_dico.keys():
		if re.search(gene_descrip_dico[gene], merge_cols(row, ["5)antibody","6)target", "11)description"]).lower()):
			return gene, "gene descr (5)"
	
	#Deuxième passe, plus large
	#itère sur le dictionnaire de gènes et compare sur les lignes les moins spécifiques
	for gene in gene_dico.keys():
		if re.search(gene_dico[gene], merge_cols(row,["cell_type", "8)strain"]).lower()):
			return gene, "gene (5)"
		elif re.search(gene_dico[gene], row["9)genotype"].lower()):
			return gene, "gene (5)"		
	#itère sur le dictionnaire d'alias et compare sur les lignes les moins spécifiques
	for gene in gene_descrip_dico.keys():
		if re.search(gene_descrip_dico[gene], merge_cols(row,["cell_type", "8)strain"]).lower()):
			return gene, "gene descr (5)"
		elif re.search(gene_descrip_dico[gene], row["9)genotype"].lower()):
			return gene, "gene descr (5)"

	#Troisième passe, moins spécifique
	#itère sur le dictionnaire d'histones
	for hist in histones_dico.keys():
		if re.search(histones_dico[hist], row["13)other"].lower()):
			return hist, "histone mark (5)"
	#itère sur le dictionnaire de gènes et compare sur les lignes les moins spécifiques
	for gene in gene_dico.keys():
		if re.search(gene_dico[gene], merge_cols(row, ["1)identifier", "13)other"]).lower()):
			return gene, "gene (6)"
	#itère sur le dictionnaire d'alias et compare sur les lignes les moins spécifiques
	for gene in gene_descrip_dico.keys():
		if re.search(gene_descrip_dico[gene], merge_cols(row, ["1)identifier", "13)other"]).lower()):
			return gene, "gene descr (6)"		

	return row["clean_target"], "target_dico (1)"


	
