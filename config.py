import os.path
import operator
from collections import OrderedDict

PREPROCESS = lambda csvfile, row: operator.setitem(row, "filename", os.path.basename(csvfile.name)) #row['filename']=os.path.basename(csvfile.name)

FIELDNAMES=OrderedDict([
	('1)identifiant',
		'sourcename'),
	('2)filename',
		'filename'),
	('3)organism',
		"organism"),
	('4)assaytype',
		'(comment\[library_selection\]|comment\[library_strategy\]|characteristics\[sampledescription\])'),
	('5)antibody',
		"(antibody|immunoprecipitate)"),
	('6)target',
		'(epitopetag|taggedprotein|taptag|protein)'),
	('7)treatment',
		'(phosphate|concentration|growth|medium|media|condition|cyclestage|developmentalstage|culturetype|transformedwith)'),
	('8)strain',
		'strain'),
	('9)genotype',
		'(genotype)'),
	('10)description',
		'(comment\[sample_source_name\]|comment\[sample_description\]|characteristics\[individual\]|comment\[sample_title\])'),
	('11)other',
		'.*'),
	])
#4,5 , 10cd ..
