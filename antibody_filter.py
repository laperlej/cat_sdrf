# -*- coding: utf-8 -*-
"""
analyser la colonne 5)antibody pour trier entre les protéines ciblées, les flag et le reste (numéro de catalogue)
analyser la colonne 6)target pour avoir la protéine cible (match avec dictionnaire de gènes de Scerevisiae ou avec dictionnaire de cibles)
analyser la colonne 8)strain et 9)genotype pour déterminer la protéine-cible avec les flag (regex +match avec dictionnaire de genes)
"""

	#si l'information est dans les 23 premières listes du dico, écrire clé dans la rangée clean_target
	
	#si info = tag_HA, 
		#si la ligne contient 'untagged',
			
		#sinon,
			#chercher regex (.+::\d?HA|\d?HA::.+|.+-HA) dans colonne 9
			#ensuite comparer match avec dictionnaire de  de S cerevisiae.

	#si info = flag, 
		
		#sinon,	
			#chercher regex (.+::flag|.+-flag) dans colonne 9
			#ensuite comparer match avec dictionnaire de  de S cerevisiae.

	#si info = tag
		#sinon,		
			#chercher regex (.+[::|-]myc|myc[::|-].+) dans colonne 9
			#ensuite comparer match avec dictionnaire de  de S cerevisiae.

	#si info = tag_V5, 
		
		#sinon,	
			#chercher regex (.+[::|-]V5|V5[::|-].+) dans colonne 9
			#ensuite comparer match avec dictionnaire de  de S cerevisiae.

	#si info = tag_T7,
	
		#sinon,	
			#chercher regex (.+[::|-]T7|T7[::|-].+) dans colonne 9
			#ensuite comparer match avec dictionnaire de  de S cerevisiae.
import re
def filter_rows(rows, target_dico, input_col, output_col):
	for row in rows:
		row = filter_row(row, target_dico, input_col, output_col) 
	return rows

def filter_row(row, target_dico, input_col, output_col):
	"""
	itère sur les regex de target-dico et les compare avec l'information dans input_col jusqu'à ce qu'il y ait un match.
	input: 
		row: dictionnaire, la clé est le titre de la colonne et la valeur est le contenu de la colonne
		target_dico: le dictionnaire de regex
		input_col: colonne dans laquelle on cherche.
		output_col: colonne changée lorsqu'il y a un match
	output:
		row: colonne output_col = clé du target_dico s'il y a eu un match
	"""
	new_value = ""
	for info in target_dico.keys():
		if re.search(target_dico[info], row[input_col].lower()):
			new_value = info
			break	
	 
	row[output_col] = new_value
	return row

	