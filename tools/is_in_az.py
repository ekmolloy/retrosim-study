import numpy
import sys


logps = []
trees = []

for line in sys.stdin:
	if line[:5] == "Total":
		logp = float(line.split()[3])
		logps.append(logp)

	if line[:14] == "CalGTProb net1":
		tree = line.split()[2].replace('(', '').replace(')', '')
		trees.append(tree)

i = numpy.argmax(logps)

if i == 0:
	sys.stdout.write("Did NOT find gene tree with higher probability\n")
else:
	sys.stdout.write("Found gene tree %s with higher probability\n" % trees[i])
