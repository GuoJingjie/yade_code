import os
import sys
import pickle
import json

#-------------------#
# Constants
#-------------------#
bigSep = "\n====================== "
sep = "---------------------- "

#-------------------#
# Import measures functions
#-------------------#
execfile("common/simulationPyRunners.py")
execfile("common/measures.py")
execfile("params_extract.py")


def read_ids(dr):
	print(sep + "Reading ids.")
	f = open(dr + '/.ids', 'r')
	ids = eval(f.read())
	f.close()
	return ids


def read_data(dr):
	print(sep + "Loading data.")
	stime = []
	data = {}
	for key in measures:
		data[key] = []
	for f in os.listdir(dr + "/data"):
		print("Loading file: " + dr + "/data/" + f)
		### Loading data.
		O.load(dr + "/data/" + f)
		### Getting time.
		stime.append(O.time)
		### Measure data.
		for key in measures:
			data[key].append(eval(measures[key]))
	return stime, data


def sort_data(stime, data):
	print(sep + "Sorting data.")
	for key in data:
		sstime, data[key] = zip(*sorted(zip(stime, data[key])))
	return sstime, data


def extract(dr):
	# Update name_value
	name_value = dr.split("_")[-1]

	#ids = read_ids(dr)
	stime, data = read_data(dr)
	stime, data = sort_data(stime, data)
	# Adding time to data
	data["time"] = stime
	# Saving :
	print(sep + dr + '/data.pkl')
	with open(dr + '/data.pkl', 'wb') as f:
		pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
	if json_save:
		print(sep + dr + '/data.json')
		with open(dr + '/data.json', 'w') as f:
			json.dump(data, f, indent=4)


#-------------------#
# Post Processing 1D
#-------------------#
for dr in sys.argv[1:]:
	print(bigSep + dr)
	os.chdir(dr)
	execfile("params.py")
	os.chdir("..")
	extract(dr)
