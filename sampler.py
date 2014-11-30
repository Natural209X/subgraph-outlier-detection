""" This samples the reSampledMotifs taking 1 in 100 uniformly """
import random

fin = open("sampledMotifs", "r")
fout = open("reSampledMotifs", "w")

for line in fin:
	x = random.randint(0, 99)
	if x==0:
		fout.write(line)

fin.close()
fout.close()
