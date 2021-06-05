import os, sys

TESTS = ""

path = sys.argv[-1]
files = os.listdir(path)

for file in files:
	# Number 7 is the number of letters that we want to skip in current naming convention "testGui".
	# -3 is the ending ".py" we do that because bash script needs only distinguishing names e.g. "Empty", "Simple"
	if(file[len(file) - 3:] == ".py"):
		TESTS += file[7:-3] + " "
print(TESTS)
