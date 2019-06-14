import re, glob, os, sys, itertools
import pandas as pd
from collections import Counter

def get_pairs(nwk):
	tt = open(nwk, 'r').read()[:-2]
	taxa = re.sub(r'[()]','',re.sub(':\d+\.\d+','',tt)).split(',')
	taxa.sort()
	# print("taxa: " + ', '.join(taxa))
	pairs = [','.join(p) for p in (itertools.combinations(taxa,2))]
	codes = ['1100', '1010', '1001', '1001', '1010', '1100']
	return dict(zip(pairs, codes))

def split_tree(tree_name):
	tt = open(tree_name, 'r').read()[1:-3]
	split = re.split(r'[()]',re.sub(':\d+\.\d+','',tt))
	taxa = [[i for i in s.split(',') if i] for s in split]
	for t in taxa:
		if len(t) == 2:
			t.sort()
			return ','.join(t)
			
data_folder = sys.argv[1]
os.chdir(data_folder)
data_name =  data_folder.split('/')[1]
nwks = glob.glob("*.nwk")
one_nwk = nwks[0]
print(f"data name: {data_name}")

types = [n.split('-')[0] for n in nwks]
all_b = list(set([int(re.split("[a-z]",t)[1]) for t in types]))
all_s = list(set([int(re.split("[a-z]",t)[2]) for t in types]))

all_s.sort()
all_b.sort()
pairs = get_pairs(one_nwk)
print("\n==== PAIRS ====")
for pair in pairs:
	print(f"{pair}:\t{pairs[pair]}")
	
data = pd.DataFrame(columns = ['1100','1010','1001'])
b_list = []
s_list = []
for b in all_b:
	b_data = [n for n in nwks if f"b{b}s" in n]
	for s in all_s:
		s_data = [n for n in b_data if f"s{s}r" in n]
		s_list.append(s)
		b_list.append(b)
		print(f"reading {b} {s} {len(s_data)}", b_data)
		result = [pairs[split_tree(n)] for n in s_data]
		counter = Counter(result)
		data = data.append(counter, ignore_index=True).fillna(0)
	
data['b'] = b_list
data['s'] = s_list

data_div = data.div(100)
data_div['b'] = b_list
data_div['s'] = s_list
summary = data_div.groupby(['b']).mean().drop("s", axis = 1)
print(summary)

os.chdir("../../")
data_div.to_excel(f"{data_name}_scores.xlsx")