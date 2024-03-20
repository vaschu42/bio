from pathlib import Path
import math


DIR = Path('blocks') # path for your folder with .txt files
aa = {} # Number of amino acids
f = {} # Amino acid frequencies
n = {} #Observed number of substitutions
O = {} #Observed frequency of substitutions
S = {} #Matrix Score
comb = 0

for file in DIR.iterdir():
	fl = open(file)
	columns = {}
	c = 0 # number of lines
	l = 0 # number of columns

	for line in fl: 
		line = line.rstrip()
		l += 1
	
		for idx, letter in enumerate(line):			
			if idx not in columns:
				columns[idx] = {}			
			columns[idx][letter] = columns[idx].get(letter, 0) + 1			
			if letter not in aa:
				aa[letter] = 1
			else:
				aa[letter] += 1
				
	total_aa = sum(aa.values())
	c += len(line)
	comb += c * ((l * (l - 1)) / 2)

	for col in columns.values():
		for k, v in col.items():
			for k2, v2 in col.items():
				if k not in n:
					n[k] = {}
				n[k][k2] = n[k].get(k2, 0) + (v * v2 if k != k2 else (v * (v - 1)) / 2)
				
	for x in n.values():
		for k, v in x.items():
			for k2, v2 in x.items():
				if k not in O:
					O[k] = {}
				O[k][k2] =  n[k].get(k2, 0) / comb
	
for a in aa: 
	f[a] = aa[a] / total_aa

E = {k: {k2: f[k] * f[k2] for k2 in f} for k in f} #Expected frequency of substitutions

for y in O.values():
	for k, v in y.items():
		for k2, v2 in y.items():
			if k not in S:
				S[k] = {}
			S[k][k2] = 2 * math.log2(O[k][k2] / E[k][k2])

print(f' ', end="")
[print(f'  {p}', end="")for p in sorted(S.keys())]
print("\t")
for k, v in sorted(S.items()):
	print(k, end="")
	for i, j in sorted(v.items()):
		print(f'  {int(round(j))}', end="") if int(round(j)) >= 0 else print(f' {int(round(j))}', end="") 
	print("\t")


fl.close()
	
