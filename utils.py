def norm_string(s):
	return s.replace(" ","").replace("-","").lower()

def norm_keys(d):
	return {norm_string(key):d[key] for key in d}
