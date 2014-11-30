from cvxopt import matrix, solvers
import itertools
import logging
import random

logger = None

def createLogger(logfile):
	global logger
	# set up logging to file - see previous section for more details
	logging.basicConfig(level=logging.DEBUG,
	                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
	                    datefmt='%m-%d %H:%M',
	                    filename=logfile,
	                    filemode='a')
	# define a Handler which writes INFO messages or higher to the sys.stderr
	console = logging.StreamHandler()
	console.setLevel(logging.INFO)
	# set a format which is simpler for console use
	formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
	# tell the handler to use this format
	console.setFormatter(formatter)
	# add the handler to the root logger
	logging.getLogger('').addHandler(console)

	logger = logging.getLogger('Simplex')
	logger.info("========================Logger created successfully!!===========================\n\n")

C=100.0
def outlierScore(neighbors,belongingnessVectorDict,fields,edges,C,line):
	""" Returns outlier score"""
	pairs=list(itertools.combinations(neighbors,2))
	Hm=[]
	Lm=[]
	w=[]
	zeta=[]
	b=[]
	c=[]

	for f in fields:
		w.append([])
	for rg in range(len(pairs)):
		zeta.append([])

	j=0
	for p in pairs:			#Equation 3.9-3.12
		n1=min(p)
		n2=max(p)

		if n2 in edges[n1]:	#Equation 3.9 and 3.10 (alternative lines for each)
			Hm.append(-1.0)
			Hm.append(0.0)
			Lm.append(0.0)
			Lm.append(0.0)
			for i in range(len(w)):
				w[i].append(1.0*abs(belongingnessVectorDict[n1][i]-belongingnessVectorDict[n2][i]))
				w[i].append(0.0)
			for i in range(len(zeta)):
				if (i==j):
					zeta[i].append(-1.0)
					zeta[i].append(-1.0)
				else:
					zeta[i].append(0.0)
					zeta[i].append(0.0)

		else:				#Equation 3.11 and 3.12 (alternative lines for each)
			Hm.append(0.0)
			Hm.append(0.0)
			Lm.append(1.0)
			Lm.append(0.0)
			for i in range(len(w)):
				w[i].append(-1.0*abs(belongingnessVectorDict[n1][i]-belongingnessVectorDict[n2][i]))
				w[i].append(0.0)
			for i in range(len(zeta)):
				if (i==j):
					zeta[i].append(-1.0)
					zeta[i].append(-1.0)
				else:
					zeta[i].append(0.0)
					zeta[i].append(0.0)

		b.append(0.0)
		b.append(0.0)


		j+=1
	for i in range(len(w)):	#Equation 3.13
		Hm.append(0.0)
		Hm.append(0.0)
		Lm.append(0.0)
		Lm.append(0.0)
		for i in range(len(zeta)):
			zeta[i].append(0.0)
			zeta[i].append(0.0)
		for j in range(len(w)):
			if(i==j):
				w[j].append(1.0)
				w[j].append(-1.0)
			else:
				w[j].append(0.0)
				w[j].append(0.0)
		b.append(1.0)
		b.append(0.0)

							#Equation 3.14
	Hm.append(0.0)
	Hm.append(0.0)
	Lm.append(0.0)
	Lm.append(0.0)
	for i in range(len(zeta)):
		zeta[i].append(0.0)
		zeta[i].append(0.0)
	for i in range(len(w)):
		w[i].append(1.0)
		w[i].append(-1.0)
	b.append(1.0)
	b.append(-1.0)

	# Now populating c, the order of coefficients in c is Hm, Lm, w, zeta
	c.append(1.0)
	c.append(-1.0)
	for i in w:
		c.append(0.0)
	for i in zeta:
		c.append(1.0*C/len(pairs))

	tempList=[]
	tempList.append(Hm)
	tempList.append(Lm)


	for i in w:
		tempList.append(i)
	for i in zeta:
		tempList.append(i)

	A=matrix(tempList)
	b=matrix(b)
	c=matrix(c)
	# Now feeding the input to simplex algo

	# print A
	# print b
	# print c

	try:
		sol=solvers.lp(c,A,b, solver='glpk')
		return sol['x'][0]-sol['x'][1]
	except Exception as e:
		logger.error("Failed to process " + str(e))
		return 0



belongingnessVectorDict={}
fields={}
edges={}

def readData():
	""" It reads graph, mapping and feature vector"""

	tempList=[]

	fin = open("map","r")
	mapping={}
	for line in fin:
		line=line.rstrip()
		tempList=line.split(" ")
		mapping[int(tempList[1])]=int(tempList[0])
	fin.close()

	fin = open("fieldProfile","r")

	i=0
	for line in fin:
		line=line.rstrip()
		tempList=line.split(" ")
		fields[tempList[0]]=i
		i+=1
	fin.close()

	fin = open("graph","r")
	for line in fin:
		line=line.rstrip()
		tempList=line.split(" ")
		if int(tempList[0]) not in edges:
			edges[int(tempList[0])]=set()
		edges[int(tempList[0])].add(int(tempList[1]))

		if int(tempList[1]) not in edges:
			edges[int(tempList[1])]=set()
		edges[int(tempList[1])].add(int(tempList[0]))

		# if int(tempList[1])<int(tempList[0]):
		# 	print "Error"
	fin.close()


	fin = open("belongingnessVectorDict","r")

	tempDict={}
	for line in fin:
		line=line.rstrip()
		tempList=line.split(" ",1)

		if int(tempList[0]) not in mapping:
			continue
		authId=mapping[int(tempList[0])]
		belongingnessVectorDict[authId]=[]
		for f in fields:
			belongingnessVectorDict[authId].append(0.0)
		tempDict=eval(tempList[1])
		for f in tempDict:
			(belongingnessVectorDict[authId])[fields[f]]=tempDict[f]

	fin.close()


def processMotifs():
	""" It processes all the motifs and writes outlier scores"""
	"""
		Input file line format: node1 node2 node3 node4 motifType
	"""
	fin = open('reSampledMotifs','r')
	fout = open('outlierScores','w')
	neighbors=set()

	for i, line in enumerate(fin):
		tempList=(line.rstrip().split(" "))
		tempList=map(int,tempList)
		neighbors.clear()
		motifAuths=tempList[:4]

		for nod in motifAuths:
			for auth in edges[nod]:
				neighbors.add(auth)

		""" 
			Since our dataset is too large and contains motifs having as many as n = 500 neigbhours.
			There as many as O(n^4) constraints. Simplex cannot handle so many constraints.
			We ignore such motifs because it is highly improbable that such motifs are outliers due to there connectedness.
		"""
		if len(neighbors)<=100:
			fout.write(str(i) + " " + str(len(neighbors))+" "+str(outlierScore(neighbors,belongingnessVectorDict,fields,edges,C, i))+"\n")

		# if(i%10000==0):
		# 	print str(i)

	fin.close()
	fout.close()


if __name__ == "__main__":
	createLogger('logger.log')
	print "Reading files........"
	readData()
	print "All files read!!"
	print "Processing motifs......"
	processMotifs()



