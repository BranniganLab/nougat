from itertools import product

def prep_config(category_names, lol_of_categories):
	""" Create a text file that the user can manually add paths to.

	Keyword arguments:
	category_names -- a list of strings containing category names
	lol_of_categories -- list of lists (lol) containing the entries for each category
	"""

	#error check
	if len(category_names) != len(lol_of_categories):
		print("different number of categories than lists!")
		return

	with open('comp_config.txt', 'w') as f:
		
		# write comments at top
		f.write("### Please enter paths to nougat-generated output folder after each entry below\n")
		f.write("### If no path exists, write NULL\n")
		f.write("### Entry format is "+';'.join(category_names)+":[path]\n")
		f.write("\n")

		#create entries for user to fill in
		for item in product(*lol_of_categories):
			f.write(";".join(item)+':\n')


def strip_blank_lines(f):
	"""Generator function to return only non-blank lines

	Keyword arguments:
	f -- the list of lines
	"""
	for l in f:
		line = l.strip()
		if line:
			yield line


def read_config(config_file, category_names):
	"""Reads in the config file provided by user and returns a dict

	Keyword arguments:
	config_file -- config filename (must be in working directory)
	category_names -- a list of strings containing category names
	"""

	config_dict = {}
	config_dict['TOC'] = {}
	counter = 0

	with open(config_file, "r+") as f:
		#skip blank lines
		for line in strip_blank_lines(f):

			#ignore lines starting with comment
			if line.startswith("#") is True:
				continue
			
			else:
				# ignore in-line comment if present
				line = line.partition('#')[0]

				#double-check that is probably unnecessary
				if line.rstrip():

					#initialize dict key to contain nested dict
					config_dict[str(counter)] = {}
					
					#split the line into 2 parts
					values = line.split(":")
					path = values[1]
					cat_vals = values[0].split(";")

					#save to nested dict
					for key,val in zip(category_names,cat_vals):
						config_dict[str(counter)][str(key)]= val
						if key in config_dict['TOC']:
							if val not in config_dict['TOC'][str(key)]:
								config_dict['TOC'][str(key)].append(str(val))
						else: 
							config_dict['TOC'][str(key)] = [str(val)]

					config_dict[str(counter)]["path"]= path
					
					counter=counter+1

	return config_dict


def generate_combinations(config_dict):
	"""Generator that reads in a config dict and yields a list of the different paths needed to make each figure

	Keyword arguments:
	config_dict -- the dictionary made by read_config
	"""

	key_list = config_dict.keys()
	category_list = config_dict['TOC'].keys()
	for category in category_list:
		for comparison_category in category_list:
			if comparison_category == category:
				continue
			else:
				
		for value in config_dict['TOC'][category]:
			combo_list = []
			for key in key_list:
				if key == 'TOC':
					continue
				else:
					if config_dict[key][category] == value:
						combo_list.append(config_dict[key]['path'])
			print(value)
			print(combo_list)

if __name__ == "__main__": 
	#prep_config(["Lipid Tail Length", "Saturation", "Structure"], [["2", "3", "4", "5", "6"],["Saturated", "Mono-unsaturated"],["capped","uncapped","protein-less"]])
	config_dict = read_config('comp_config.txt',["length","saturation","structure"])
	generate_combinations(config_dict)