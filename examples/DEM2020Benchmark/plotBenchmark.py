import matplotlib as mpl
mpl.use('Agg')  # so it runs without a X server (e.g. HPC)

import matplotlib.pyplot as plt
import os

for case in ['Case1', 'Case2', 'Case3']:
	files = [file for file in os.listdir("./outputData") if (file.startswith(case) and file.endswith(".txt"))]
	if len(files) == 0:
		continue  # no data for this case
	for f in files:
		X, Y, titles = [], [], []
		for line in open("./outputData/" + f, 'r'):
			if len(titles) == 0:
				titles = [str(s) for s in line.split()]
				continue
			values = [float(s) for s in line.split()]
			X.append(values[0] if case != 'Case1' else values[-1])
			Y.append(values[1])
		plt.plot(X, Y, label=f)
	plt.xlabel(titles[1 if case != 'Case1' else -1])
	plt.ylabel(titles[2])
	plt.legend()
	print("saving", case, "...")
	plt.savefig(case + '.png')
	plt.clf()
