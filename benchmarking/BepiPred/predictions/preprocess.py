#This file processes the raw results from the Bepipred method.
import sys
import numpy as np
inputfile = open(sys.argv[1],'r')
lines = inputfile.readlines()
dic={}
current=''
for i in range(len(lines)):
	line = lines[i]
	if line[0:2]=='##':
		current = line.split()[-1].strip()
		#print(current, lines[(i+3)].split())
		dic[line.split()[-1].strip()]=[float(lines[i+3].split()[5].strip())]
		#print(dic)
	elif line[0] != '#':
		#print(current,dic[current])
		dic[current].append(float(lines[i].split()[5].strip()))
	else:
		continue
#print(dic)
inputfile.close()

output = open(sys.argv[1]+".cleaned", 'w')
print("id"+"\t"+"target"+"\t"+"pred"+"\t"+"pred-score")
for i in dic.keys():
	pred_score = np.average(np.array(dic[i]))
	print(str(i)+"\t"+str(i.split('-')[-1])+"\t"+str("1" if pred_score >= 0.35 else "0")+"\t"+str(pred_score), file=output)



