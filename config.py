import os.path
import operator
from collections import OrderedDict

PREPROCESS = lambda csvfile, row: operator.setitem(row, "filename", os.path.basename(csvfile.name)) #row['filename']=os.path.basename(csvfile.name)

FIELDNAMES=OrderedDict([
	('1)identifier',
		'sourcename'),
	('2)filename',
		'filename'),
	('3)organism',
		"organism"),
	('4)assaytype',
		'(comment\[library_selection\]|comment\[library_strategy\]|characteristics\[sampledescription\])'),
	('5)antibody',
		"(antibody|milliporecatno)"),
	('6)target',
		'(epitopetag|tagged|taptag|protein|h2b|immunoprecipitate|target|\[tag\])'),
	('7)treatment',
		'(phosphate|concentration|growth|medium|media|condition|cycle|stage|developmental|stage|culturetype|transformedwith|treatment|temperature|sac_cerandgla_glamixedpercentage|compound)'),
	('8)strain',
		'(strain|cellline)'),
	('9)genotype',
		'(genotype|genedeletion|\[variation\]|genetic|\[yrr1alleletransformed\])'),
	('10)platform',
		'(platform|instrument_model)'),
	('11)description',
		'(comment\[sample_source_name\]|comment\[sample_description\]|characteristics\[individual\]|comment\[sample_title\])'),
	('12)fastq',
		'fastq_uri'),
	('13)other',
		'.*')
	])
#4,5 , 10cd ..
