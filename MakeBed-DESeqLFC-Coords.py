#get coordinates for a list of gene names

file1="/path/to/gene/list/DE_genes.txt"
file2="/path/to/promoter/coordinates/UCSC-TSS-mm10_PROMOTER-COORDS.txt"

array={}



with open(file2,'r') as f:
	for line in f:
		f=line.strip().split('\t')
		array[f[3]]=[f[0],f[1],f[2]]

#bb are trying to print columns in file 1. bb[0] is column 1, add more to print step if more columns
with open(file1,'r') as f1:
	for line in f1:
		bb=line.strip().split('\t')
		if(bb[0] in array):
			(a,b,c)=array[bb[0]]

			print (a, b, c, bb[0],sep='\t')


		else:
			pass

