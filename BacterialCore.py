# ETE implements the three most common strategies: preorder, levelorder and postorder.
# 
#     preorder: 1)Visit the root, 2) Traverse the left subtree , 3) Traverse the right subtree.
#     postorder: 1) Traverse the left subtree , 2) Traverse the right subtree, 3) Visit the root
#     levelorder (default): every node on a level before is visited going to a lower level
# 
# Note
# 
#     Preorder traversal sequence: F, B, A, D, C, E, G, I, H (root, left, right)
#     Inorder traversal sequence: A, B, C, D, E, F, G, H, I (left, root, right); note how this produces a sorted sequence
#     Postorder traversal sequence: A, C, E, D, B, H, I, G, F (left, right, root)
#     Level-order traversal sequence: F, B, G, A, D, I, C, E, H
# analtype=3; cutoff=0.01; file='/home/silviatm/micro/bc/tomate_subset.fa'; final_lim=0.74; ilvl=0.99; compute_core=1.0; output='/home/silviatm/micro/bc/tomate_bc_0.01/'; process=3; taxo_p=0.9; threads=16; tree_l=99
# file_no_fa=os.path.basename(str(file.split(".")[0]))+"_otus.txt" ; tree_l_name="99_otus_nodes.tree" ; init_var_tree=0.99

from ete3 import Tree #ete3 is a library needed to transverse a tree object
import numpy as np
import pandas as pd # new dependency...
from scipy import stats
import argparse
import os
import fileinput

def main():  
  ##########################################################################################
	#--Argument
	parser = argparse.ArgumentParser(description="")

        #--Input/output files
	# Must use absolute paths
	parser.add_argument('-f', dest = 'file', type = str, default = './', required = True, help = "Reads .fasta file") #Fasta file with all the sequences from all subjects of a certain condition
	parser.add_argument('-o', dest = 'output', type = str, default = './', required = True, help = "Output file") #Folder name where output files are going to be located
	parser.add_argument('-p', dest = 'process', type = int, required = True, help = "1-Complete process, OTUs and Tree. 2- OTUS. 3- Tree. 4- Venn comparation: ") #Choice which process do you want to be executed
	parser.add_argument('-initial_level', dest = 'ilvl', type = float, default = 0.97, choices=[Range(0.5, 1.0)], required = False, help = "First aggregation level to measure") #Initial identity level to use in pick_otus
	parser.add_argument('-final_level', dest = 'flvl', type = float, default = 0.74, choices=[Range(0.5, 1.0)], required = False, help = "Last aggregation level to measure") #Las identity level to use in pick_otus
	parser.add_argument('-tree_type_analysis', dest = 'analtype', type = int, default = 1, required = False, choices=[1,2,3,4], help = "'1: Core 100. 2: Binomial. 3: Min Percentage. 4: Print all nodes'") #Way to evaluate core levels while trasversing the 97_nodes tree
	parser.add_argument('-tree_level', dest = 'tree_l', type = int, default = 97, required = False, choices=(97, 99), help = "97: 97_otus_nodes.tree, 99: 99_otus_nodes.tree. Default: 97") #Choice which level
	parser.add_argument('-min_core_percentage', dest = 'min_core', type = float, default = 1, required = False, choices=[Range(0.5, 1.0)], help = "Minimum percentage of samples to describe a core") #Here we can change the "1" (100%) value for another percentage
	parser.add_argument('-cutoff', dest = 'cutoff', type = float, default = 0, required = False, choices=[Range(0.0, 1.0)], help = "Minimum relative abundance needed to consider a node to be 'present' in one sample (Tree approach, process==3). Default: 0")
	parser.add_argument('-taxo_p', dest = 'taxo_p', type = float, default = 0.9, required = False, choices=[Range(0.5, 1.0)],	help = "minimum percentage of the same taxonomic group within all OTUs contained into the same Node") #Here we can change the "1" (100%) value for another percentage
	parser.add_argument('-threads', dest = 'threads', type = int, default = 1, required = False, help = "For paralellization of all parallelizable steps")
	
	args = parser.parse_args()
  ##########################################################################################

	cores = []

	init_var=args.ilvl #Saving the arguments as new variables
	last_var=init_var
	final_lim=args.flvl
	output = args.output
	tree_l = args.tree_l
	if tree_l==97:
		tree_l_name="97_otus_nodes.tree"
		init_var_tree=0.97
	elif tree_l==99:
		tree_l_name="99_otus_nodes.tree"
		init_var_tree=0.99
	compute_core=args.min_core
	cutoff=args.cutoff
	taxo_p=args.taxo_p
	threads = str(args.threads)
	
	file_no_fa=os.path.basename(str(args.file.split(".")[0]))+"_otus.txt" #Create a variable that has the variable name summed to the sufix _otus.txt, to permit pick_rep_set finde the file generated from pick_otus

	if args.process==1 or  args.process==2: #The denovo OTUs picking will be performed if we choose the complete analysis (1) or only the OTU picking (2)
		os.system("pick_otus.py -i "+str(args.file)+" -o "+str(output)+"/OTUs/"+str(init_var)+"/ -m usearch61 -z -s "+str(init_var)+" --threads "+threads) #Denovo OTUs picking at the initial identity level, using usearch61. This version is necessary for large datasets, but could be changed for usearch (32 bits version) for smaller ones
		os.system("pick_rep_set.py -i "+str(output)+"/OTUs/"+str(init_var)+"/"+str(file_no_fa)+" -f "+str(args.file)+" -o "+str(output)+"/OTUs/"+str(init_var)+"/rep_set.fasta") #At this point a representative sequence is choosen within each OTU, with default parameters
		os.system("./prinseq-lite.pl -fasta "+str(output)+"/OTUs/"+str(init_var)+"/rep_set.fasta -out_format 1 -out_good "+str(output)+"/OTUs/"+str(init_var)+"/rep_set_"+str(init_var)+".good") #We need prinseq in order to chage the format of the multifasta, to avoid sequence import problems with the next steps
		os.system("assign_taxonomy.py -i "+str(output)+"/OTUs/"+str(init_var)+"/rep_set_"+str(init_var)+".good.fasta -o "+str(output)+"/OTUs/"+str(init_var)+"/rep_set.good."+str(init_var)+".assignedTax -m rdp") #Here we assign taxonomy to each OTU using rdp classifier, with default parameters
#		os.system("parallel_assign_taxonomy_rdp.py -i "+str(output)+"/OTUs/"+str(init_var)+"/rep_set_"+str(init_var)+".good.fasta -o "+str(output)+"/OTUs/"+str(init_var)+"/rep_set.good."+str(init_var)+".assignedTax --jobs_to_start "+threads) # HAS TO BE AN ABSOLUTE PATH
		os.system("make_otu_table.py -i "+str(output)+"/OTUs/"+str(init_var)+"/"+str(file_no_fa)+" -t "+str(output)+"/OTUs/"+str(init_var)+"/rep_set.good."+str(init_var)+".assignedTax/rep_set_"+str(init_var)+".good_tax_assignments.txt -o "+str(output)+"/OTUs/"+str(init_var)+"/All_table_"+str(init_var)+".biom") #Using the taxonomy and representative OTUs file, we create here an OTU table
		os.system("biom convert -i "+str(output)+"/OTUs/"+str(init_var)+"/All_table_"+str(init_var)+".biom -o "+str(output)+"/OTUs/"+str(init_var)+"/table.from_biom_"+str(init_var)+".txt --to-tsv --header-key taxonomy") #This is an extra step to generate an OTU table as .txt, in case we need it for other to check the results at this point
		os.system("compute_core_microbiome.py -i "+str(output)+"/OTUs/"+str(init_var)+"/All_table_"+str(init_var)+".biom -o "+str(output)+"/OTUs/"+str(init_var)+"/All_table_normKO_core"+str(int(compute_core*100))+" --min_fraction_for_core "+str(compute_core)+" --max_fraction_for_core "+str(compute_core)) #Here we compute the core microbiome
		os.system("biom convert -i "+str(output)+"/OTUs/"+str(init_var)+"/All_table_normKO_core"+str(int(compute_core*100))+"/core_table_"+str(int(compute_core*100))+".biom -o "+str(output)+"/OTUs/"+str(init_var)+"/All_table_normKO_core"+str(int(compute_core*100))+"/"+str(init_var)+".biom --to-json") #In this and the next step we convert the biom files to json format, needed for further analysis
		os.system("biom convert -i "+str(output)+"/OTUs/"+str(init_var)+"/All_table_"+str(init_var)+".biom -o "+str(output)+"/OTUs/"+str(init_var)+"/All_table_json_"+str(init_var)+".biom --to-json")

		if(os.path.exists(str(output)+"/OTUs/"+str(init_var)+"/All_table_normKO_core"+str(int(compute_core*100))+"/"+str(init_var)+".biom")): #With this if loop we modify some strings in the biom table to help the visualization of the results once the statistics are done
			os.system("cp "+str(output)+"/OTUs/"+str(init_var)+"/All_table_normKO_core"+str(int(compute_core*100))+"/"+str(init_var)+".biom "+str(output)+"/OTUs/"+str(init_var)+"/All_table_normKO_core"+str(int(compute_core*100))+"/"+str(init_var)+"._changed_name.biom")
			for line in fileinput.input(str(output)+"/OTUs/"+str(init_var)+"/All_table_normKO_core"+str(int(compute_core*100))+"/"+str(init_var)+"._changed_name.biom", inplace=True):
				print (line.replace("]}},{\"id\": \"", "]}},{\"id\": \"Core "+str(init_var)+" "))
			for line in fileinput.input(str(output)+"/OTUs/"+str(init_var)+"/All_table_normKO_core"+str(int(compute_core*100))+"/"+str(init_var)+"._changed_name.biom", inplace=True):
				print (line.replace("\"rows\": [{\"id\": \"","\"rows\": [{\"id\": \"Core "+str(init_var)+" "))

		os.system("mv "+str(output)+"/OTUs/"+str(init_var)+"/"+str(file_no_fa)+" "+str(output)+"/OTUs/"+str(init_var)+"/rep_set_"+str(init_var)+".good_otus.txt") #Here we modified the name of the first otus_file in order to generate and standerd name that the loop could understand in the next levels
		par_a=str(output)+"/OTUs/"+str(init_var)+"/All_table_normKO_core"+str(int(compute_core*100))+"/core_otus_"+str(int(compute_core*100))+".txt"
		par_b=str(output)+"/OTUs/"+str(init_var)+"/rep_set_"+str(init_var)+".good_otus.txt"
		par_c=str(output)+"/OTUs/filtrar"+str(init_var)+".txt"
		CoreFilt_a(par_a,par_b,par_c) #With this module we can extract the reads that are part of the core at this level
		
		output_file=open(str(output)+"/OTUs/"+str(init_var)+"/rep_set_"+str(init_var)+".good_otus.txt.modified.txt", 'w') #This line modifies the \t with ; in the OTUs file, needed to perform the next statisticals analysis with R.  

		for line in open(str(output)+"/OTUs/"+str(init_var)+"/rep_set_"+str(init_var)+".good_otus.txt", 'r'):  #With this loop we write the last file
			line = line.rstrip()
			line_info=line.split('\t')
			output_file.write(str(line_info[0])+'\t')
			for x in range (1,len(line_info)):
				output_file.write(str(line_info[x])+";")
			output_file.write("\n")

		for x in np.arange (init_var-0.01,final_lim,-0.01): #With this for loop we permorf exactly the same as before between the initial and final levels desired, with the particularity that we avoid to evaluate reads that have been already part of a core in a last aggregation level
			x=round(x,2)
			os.system("pick_otus.py -i "+str(output)+"/OTUs/"+str(last_var)+"/rep_set_"+str(last_var)+".good.fasta -o "+str(output)+"/OTUs/"+str(x)+"/ -m usearch61 -z -s "+str(x))
			os.system("merge_otu_maps.py -i "+str(output)+"/OTUs/"+str(last_var)+"/rep_set_"+str(last_var)+".good_otus.txt,"+str(output)+"/OTUs/"+str(x)+"/rep_set_"+str(last_var)+".good_otus.txt -o "+str(output)+"/OTUs/"+str(x)+"/otus"+str(x)+".txt")
			os.system("mv "+str(output)+"/OTUs/"+str(x)+"/otus"+str(x)+".txt "+str(output)+"/OTUs/"+str(x)+"/rep_set_"+str(x)+".good_otus.txt")
			os.system("pick_rep_set.py -i "+str(output)+"/OTUs/"+str(x)+"/rep_set_"+str(x)+".good_otus.txt -f "+str(args.file)+" -o "+str(output)+"/OTUs/"+str(x)+"/rep_set.fasta")
			os.system("./prinseq-lite.pl -fasta "+str(output)+"/OTUs/"+str(x)+"/rep_set.fasta -out_format 1 -out_good "+str(output)+"/OTUs/"+str(x)+"/rep_set_"+str(x)+".good")
			os.system("assign_taxonomy.py -i "+str(output)+"/OTUs/"+str(x)+"/rep_set_"+str(x)+".good.fasta -o "+str(output)+"/OTUs/"+str(x)+"/rep_set.good."+str(x)+".assignedTax -m rdp")
			os.system("make_otu_table.py -i "+str(output)+"/OTUs/"+str(x)+"/rep_set_"+str(x)+".good_otus.txt -t "+str(output)+"/OTUs/"+str(x)+"/rep_set.good."+str(x)+".assignedTax/rep_set_"+str(x)+".good_tax_assignments.txt -o "+str(output)+"/OTUs/"+str(x)+"/All_table_"+str(x)+".biom")
			os.system("biom convert -i "+str(output)+"/OTUs/"+str(x)+"/All_table_"+str(x)+".biom -o "+str(output)+"/OTUs/"+str(x)+"/table.from_biom_"+str(x)+".txt --to-tsv --header-key taxonomy")
			
			last_var=x

			os.system("cat "+str(output)+"/OTUs/filtrar*.txt > "+str(output)+"/OTUs/filt.txt")
			par_a=str(output)+"/OTUs/filt.txt"
			par_b=str(output)+"/OTUs/"+str(x)+"/rep_set_"+str(x)+".good_otus.txt"
			par_c=str(output)+"/OTUs/"+str(x)+"/otus"+str(x)+".Filt.txt"
			CoreFilt_b(par_a,par_b,par_c) #Here we ommit reads that already have been part of a core
			os.system("make_otu_table.py -i "+str(output)+"/OTUs/"+str(x)+"/otus"+str(x)+".Filt.txt -t "+str(output)+"/OTUs/"+str(x)+"/rep_set.good."+str(x)+".assignedTax/rep_set_"+str(x)+".good_tax_assignments.txt -o "+str(output)+"/OTUs/"+str(x)+"/"+str(x)+".Filt_table.biom")
			os.system("compute_core_microbiome.py -i "+str(output)+"/OTUs/"+str(x)+"/"+str(x)+".Filt_table.biom -o "+str(output)+"/OTUs/"+str(x)+"/coreFilt."+str(x)+" --num_fraction_for_core_steps 2 --min_fraction_for_core "+str(compute_core)+" --max_fraction_for_core "+str(compute_core))
			os.system("biom convert -i "+str(output)+"/OTUs/"+str(x)+"/coreFilt."+str(x)+"/core_table_"+str(int(compute_core*100))+".biom -o "+str(output)+"/OTUs/"+str(x)+"/coreFilt."+str(x)+"/"+str(x)+".biom --to-json")
			os.system("biom convert -i "+str(output)+"/OTUs/"+str(x)+"/All_table_"+str(x)+".biom -o "+str(output)+"/OTUs/"+str(x)+"/All_table_json_"+str(x)+".biom --to-json")
			
			if(os.path.exists(str(output)+"/OTUs/"+str(x)+"/coreFilt."+str(x)+"/"+str(x)+".biom")):
				os.system("cp "+str(output)+"/OTUs/"+str(x)+"/coreFilt."+str(x)+"/"+str(x)+".biom "+str(output)+"/OTUs/"+str(x)+"/coreFilt."+str(x)+"/"+str(x)+"._changed_name.biom")
				for line in fileinput.input(str(output)+"/OTUs/"+str(x)+"/coreFilt."+str(x)+"/"+str(x)+"._changed_name.biom", inplace=True):
					print (line.replace("]}},{\"id\": \"", "]}},{\"id\": \"Core "+str(x)+" "))
				for line in fileinput.input(str(output)+"/OTUs/"+str(x)+"/coreFilt."+str(x)+"/"+str(x)+"._changed_name.biom", inplace=True):
					print (line.replace("\"rows\": [{\"id\": \"","\"rows\": [{\"id\": \"Core "+str(x)+" "))
			
			par_a=str(output)+"/OTUs/"+str(x)+"/coreFilt."+str(x)+"/core_otus_"+str(int(compute_core*100))+".txt"
			par_b=str(output)+"/OTUs/"+str(x)+"/rep_set_"+str(x)+".good_otus.txt"
			par_c=str(output)+"/OTUs/filtrar"+str(x)+".txt"
			CoreFilt_a(par_a,par_b,par_c)

			output_file=open(str(output)+"/OTUs/"+str(x)+"/rep_set_"+str(x)+".good_otus.txt.modified.txt", 'w')

			for line in open(str(output)+"/OTUs/"+str(x)+"/rep_set_"+str(x)+".good_otus.txt", 'r'):
				line = line.rstrip()
				line_info=line.split('\t')
				output_file.write(str(line_info[0])+'\t')
				for x in range (1,len(line_info)):
					output_file.write(str(line_info[x])+";")
				output_file.write("\n")		


		os.system("cat "+str(output)+"/OTUs/filtrar*.txt > "+str(output)+"/OTUs/filt.txt") #Core reads filtered after each evaluated levels are concatenated into a single file
		os.system("rm -rf "+str(output)+"/OTUs/filtrar*.txt")
		os.system("Rscript --vanilla OTUs_statistics.R "+str(output)+"/OTUs/ "+str(init_var)) #Finally, the R script "OTUs_statistics" perform the statistical analysis of the cores obtained within all levels
#		print("pick_otus.py -i "+str(args.file)+" -o "+str(output)+"/Tree/"+str(init_var_tree)+"/ -m usearch61_ref -C -z -s "+str(init_var_tree)) #OTUs picking at the initial identity level, using u$
#		print("pick_rep_set.py -i "+str(output)+"/Tree/"+str(init_var_tree)+"/"+str(file_no_fa)+" -f "+str(args.file)+" -o "+str(output)+"/Tree/"+str(init_var_tree)+"/rep_set.fasta")
#		print("./prinseq-lite.pl -fasta "+str(output)+"/Tree/"+str(init_var_tree)+"/rep_set.fasta -out_format 1 -out_good "+str(output)+"/Tree/"+str(init_var_tree)+"/rep_set_"+str(init_var_tree)+".good")
#		print("assign_taxonomy.py -i "+str(output)+"/Tree/"+str(init_var_tree)+"/rep_set_"+str(init_var_tree)+".good.fasta -o "+str(output)+"/Tree/"+str(init_var_tree)+"/rep_set.good."+str(init_var_tree)+".assignedTax -m rdp")
	if args.process==1 or  args.process==3: #The analysis of the tree will be performed if we choose the complete analysis (1) or only the tree analysis (3)
#
		os.system("pick_otus.py -i "+str(args.file)+" -o "+str(output)+"/Tree/"+str(init_var_tree)+"/ -m usearch61_ref -C -z -s "+str(init_var_tree)+" --threads "+threads)  #OTUs picking at the initial identity level, using usearch61 against reference. This version is necessary for large datasets, but could be changed for usearch (32 bits version) for smaller ones
		os.system("pick_rep_set.py -i "+str(output)+"/Tree/"+str(init_var_tree)+"/"+str(file_no_fa)+" -f "+str(args.file)+" -o "+str(output)+"/Tree/"+str(init_var_tree)+"/rep_set.fasta") #At this point a representative sequence is choosen within each OTU, with default parameters
		os.system("./prinseq-lite.pl -fasta "+str(output)+"/Tree/"+str(init_var_tree)+"/rep_set.fasta -out_format 1 -out_good "+str(output)+"/Tree/"+str(init_var_tree)+"/rep_set_"+str(init_var_tree)+".good")
		os.system("assign_taxonomy.py -i "+str(output)+"/Tree/"+str(init_var_tree)+"/rep_set_"+str(init_var_tree)+".good.fasta -o "+str(output)+"/Tree/"+str(init_var_tree)+"/rep_set.good."+str(init_var_tree)+".assignedTax -m rdp")
#		os.system("parallel_assign_taxonomy_rdp.py -i "+str(output)+"/Tree/"+str(init_var_tree)+"/rep_set_"+str(init_var_tree)+".good.fasta -o "+str(output)+"/Tree/"+str(init_var_tree)+"/rep_set.good."+str(init_var_tree)+".assignedTax --poll_directly --jobs_to_start "+threads) # HAS TO BE AN ABSOLUTE PATH
		os.system("make_otu_table.py -i "+str(output)+"/Tree/"+str(init_var_tree)+"/"+str(file_no_fa)+" -t "+str(output)+"/Tree/"+str(init_var_tree)+"/rep_set.good."+str(init_var_tree)+".assignedTax/rep_set_"+str(init_var_tree)+".good_tax_assignments.txt -o "+str(output)+"/Tree/"+str(init_var_tree)+"/All_table_"+str(init_var_tree)+".biom") #Using the taxonomy and representative OTUs file, we create here an OTU table
		os.system("biom convert -i "+str(output)+"/Tree/"+str(init_var_tree)+"/All_table_"+str(init_var_tree)+".biom -o "+str(output)+"/Tree/"+str(init_var_tree)+"/table.from_biom_"+str(init_var_tree)+".txt --to-tsv --header-key taxonomy")
		os.system("sed 's/\"//g' "+str(output)+"/Tree/"+str(init_var_tree)+"/table.from_biom_"+str(init_var_tree)+".txt > "+str(output)+"/Tree/"+str(init_var_tree)+"/table.from_biom_no_com_"+str(init_var_tree)+".txt")

		par_a=str(tree_l_name)
		par_b=str(output)+"/Tree/"+str(init_var_tree)+"/table.from_biom_no_com_"+str(init_var_tree)+".txt"
		par_c=str(output)+"/Tree/results.txt"
		par_d=args.analtype
		par_e=str(output)+"/Tree/abundances.txt"
		Tree_analysis(par_a,par_b,par_c,par_d,par_e,compute_core,cutoff,taxo_p) #The "Tree analysis" module is used to perform th search of the core Nodes thorugh the tree using the OTUs abundances obtained with the otu picking against reference.
	
		output_file=open(str(output)+"/Tree/"+str(init_var_tree)+"/"+str(file_no_fa)+".modified.txt", 'w')

		for line in open(str(output)+"/Tree/"+str(init_var_tree)+"/"+str(file_no_fa), 'r'):
			line = line.rstrip()
			line_info=line.split('\t')
			output_file.write(str(line_info[0])+'\t')
			for x in range (1,len(line_info)):
				output_file.write(str(line_info[x])+";")
			output_file.write("\n")

		os.system("perl Distances_16S.pl "+str(output)+"/Tree/ "+str(tree_l)+" "+threads)
		
	if args.process==1 or  args.process==4: #This option is included to perform the Venn analysis using the results obtained with OTUs and Tree functions
		os.system("Rscript --vanilla Venn_OTUS_Tree.R "+str(output))
	
	
def CoreFilt_a(r,s,o):
	ids=[]
	with open(r,'r') as f:
		for l in f:
			l=l.strip()
			if not l.startswith("#"):
				ids.append(l.split()[0])
  
	with open(s,'r') as f, open(o,'w') as fout:
		for l in f:
			l=l.strip()
			spl=l.split()
			if spl[0] in ids:
				print >> fout, " ".join(l.split()[1:])

def CoreFilt_b(s,w,o):
	seqs=[]
	with open(s,'r') as f:
		for l in f:
			l=l.strip()
			seqs+=l.split()
	seqs = set(seqs)

	with open(w,'r') as fin, open(o,'w') as fout:
		for l in fin:
			l=l.rstrip()
			A=set(l.split()[1:])
			if not A & seqs:
				print >> fout, l

def Tree_analysis(tree,tabla,out,analysis_type,out2,compute_core=1,cutoff=0,taxo_p=0.9):  

	###Al subsequents variables could be modified
	binomial_value = float(0.05) #Default value for the option 2 of the core evaluation method for the tree
	p_value = float(0.05) #p-value threeshold for the binomial method (2 method) 
	percentage = float(compute_core) #minimum percentage threeshold of subjects requiered to defined a core
	cutoff = float(cutoff) #minimum relative abundance needed to consider a OTU to be 'present' in one sample; OTU abundance filter.
	taxo_p = float(taxo_p) #minimum percentage of the same taxonomic group within all OTUs contained into the same Node
	
	output_file=open(out, 'w')
	output_file_2=open(out2, 'w')	

	tree = Tree(tree, quoted_node_names=True, format=1) #Here we load the 97_otus/99_otus tree
	table     = {} ## abund
	table2    = {} ## tax
	
	for line in open(tabla):
		if (line.startswith('#')):
			output_file_2.write(str(line))
		else:
	  		## Abund:
			fields = list(map(str.strip, line.split('\t'))) #We create a dictionary with all the keys and values of the OTU table against reference
			table[fields[0]] = list(map(float, fields[1:-1]))
			
			## Tax:
			table2[fields[0]] = list(map(str, fields[(len(fields)-1):len(fields)])) # Here we load a dictionary with the taxonomy information from the picked OTUs


	table_final_res = [0] * len(fields[1:-1]) # empty; it's as long as samples there are. Will contain the info from the last line of results.txt
	table_final_res = ([float(i) for i in table_final_res])
	sum_abun_rela = 0
	cores = 0
	
	for leaf in tree:
		if leaf.name not in table:
			leaf.vector = None
		else:
			leaf.vector = table[leaf.name] # [rows/OTUs from abs abund table] Create value vectors for each of the tree tips of the tree with the values of the OTU table previously generated

	node2content = tree.get_cached_content()

	flag=0
	for node in tree.traverse(): #This loop is used to add values into the vectors created before
		if not node.is_leaf(): # if it's a leaf, it already has a .vector attribute.
			leaf_vectors     = np.array([leaf.vector for leaf in node2content[node] if leaf.vector is not None])
			node.vector      = leaf_vectors.sum(axis=0) # sum all the abundances of all leaves from that node
			
			if(flag == 0): # when flag==0, we're looking at the root. that means adding ALL leaves of the tree.
			  	save_node1=node.vector
	  			total_saved_leaves = np.array([leaf.name for leaf in node2content[node]])
		  		flag=1			

	if(analysis_type==4): #This method only prints the information of the tree, only for information of the tree purpouse
		print(tree.get_ascii(show_internal=True))
		output_file.write(tree.get_ascii(show_internal=True) + '\n' + '\n')
		for node in tree.traverse("preorder"):
			print (node.name, node.vector)
			output_file.write(node.name + '\t' + str(node.vector) + '\n')

	if(analysis_type==1 or analysis_type==2 or analysis_type==3): #Here we evaluate the tree traversally using the chosen method: 100% core, binomial or percentage

		output_file.write("Core" + '\t' + "Prevalence" + '\t' + "Abundance" + '\t' + "Relative abundances" + '\t' + "Min" + '\t' + "Max" + '\t' + "Average" + '\t' + "SD" + '\t' + "Leaves" + '\t' + "Taxonomy" + '\t' + "Leaves number" + '\n') 

		for node in tree.traverse("postorder"):
			if type(node.vector) == list: # empty node.vectors are not filterable but we want them to be checked...
				node.vector=([float(i) for i in node.vector]) #Transform all the values contained in node.vector to float, to perform operations efficiently 
				abundance=node.vector/save_node1 #Relative abundance of each subject in the node over the terminal node (sum of all nodes)
				abundance =([float(i) for i in abundance])
				abundance=node.vector/save_node1 #Relative abundance of each subject in the node over the terminal node (sum of all nodes)
			else:
				continue

			filtered_n_vector=[] # Apply cutoff: a node won't be considered present at a given sample if the total abundance of its leaves is smaller than the cutoff
			for i, v in zip(node.vector, abundance):
				if v < cutoff:
					filtered_n_vector.append(0.0)
				else:
					filtered_n_vector.append(i)
		
			tot_cont=float(np.count_nonzero(filtered_n_vector)) #Count the number of subjects (samples) in this study with one ore more ocurrence in the vector for a certain node 
			tot_cont2=float(np.asarray(filtered_n_vector).size) #Count the total vector array size
			
			if analysis_type==2:
				a=stats.binom_test(tot_cont, n=tot_cont2, p=binomial_value, alternative='greater') #Binomial test that uses the binomial_value
			rela=(tot_cont/tot_cont2) # percentage of samples in which the node is present
			
			if((analysis_type==1 and np.all(node.vector)) or (analysis_type==2 and a <= p_value) or (analysis_type==3 and rela >= percentage)): #Depending on the method used to go through the tree, we will evaluate different parameters to check if the node should be or not taken into account
				node.vector=([float(i) for i in node.vector]) #Transform all the values contained in node.vector to float, to perform operations efficiently 
				abundance=node.vector/save_node1 #Relative abundance of each subject in the node over the terminal node (sum of all nodes)
				abundance =([float(i) for i in abundance])
				mean_abun=np.mean([float(i) for i in abundance]) #Mean abundance of the node
				std_abun=np.std([float(i) for i in abundance]) #Standard deviation of the node
				abundance_rela=sum(node.vector)/sum(save_node1) #Global relative abundance of the node over the terminal node
				table_final_res=list(map(sum, zip(table_final_res, abundance))) #Getting all the results for each node into a final result table
				sum_abun_rela=sum_abun_rela+abundance_rela #The sum of all global relative abundance
				cores=cores+1 #Total number of cores
				
				node2content = tree.get_cached_content()
				
				output_file_2.write(str(node.name) + '\t')
				for x in range(len(abundance)): 
					output_file_2.write(str(abundance[x]) + '\t'), output_file_2.write('\n')
				
				output_file.write(node.name + '\t' +  str(rela) + '\t' + str(node.vector) + '\t' + str(abundance) + '\t' + str(min(abundance)) + '\t' + str(max(abundance)) + '\t' + str(mean_abun) + '\t' + str(std_abun) + '\t')
								
				conteo_hojas=nodes_eval(node,tree,output_file,table2,taxo_p,total_saved_leaves) #With this line we can assign a taxonomy to each node based in the taxonomy of each OTU, depending on the minimum taxonomy percentage level stablished before 
				
				output_file.write(str(conteo_hojas) + '\n') #Print the total number of leaves of this node
				
				tree=erase_node(node, tree) #Once a node has been evaluated, this line deletes that node from the tree to simplify the calculations of the next nodes
			
				G = tree.search_nodes(name=node.name)[0]
				removed_node = G.detach()
		# This is the last line, where we don't have the node name but instead the "number of nodes" (cores)
		output_file.write(str(cores) + '\t' + '\t' + '\t' + str(table_final_res) + '\t' + str(min(table_final_res)) + '\t' + str(max(table_final_res)) + '\t' + str(np.mean([float(i) for i in table_final_res])) + '\t' + str(np.std([float(i) for i in table_final_res])) + '\n')
		
class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __eq__(self, other):
        return self.start <= other <= self.end
					
def erase_node(node,tree):	
	node = tree.search_nodes(name=node.name)[0]
	node_borrar = node
	while node:
		node2 = node.up
		if(node2):
			node2.vector=node2.vector-node_borrar.vector
		node = node.up
	return tree
	
def nodes_eval(node,tree,output_file,table2,taxo_p,total_saved_leaves):
	node = tree.search_nodes(name=node.name)[0]
	node2content = tree.get_cached_content()
	total_leaves = list(node2content[node])
	
	cont_leaves=0 #Counter of the total number of leaves of the node
	dicc_tax = {} #This dictionary will keep the information of the number of OTUs that have a certain taxonomy group assigned

	for node_im in total_leaves[0:(len(total_leaves))]:
			
		if(node_im.name in total_saved_leaves and node_im.name in table2):
			cont_leaves=cont_leaves+1 #Summing total leaves number
	
			output_file.write(node_im.name + ';')
			
			taxa_list=list(map(str.strip, str(table2[node_im.name]).split("'")))[1]
			taxa_list=list(map(str.strip, taxa_list.split(';'))) #Here me split the table that has the taxonomy information
			
			
			for taxo in taxa_list[0:(len(taxa_list))]:  #Here we count how many times we find a certain a taxonomic range
			
				dicc_tax[taxo]=dicc_tax.get(taxo, 0) + 1
		
	
	output_file.write('\t')
	
	highest_valu_k = 0 #We will evaluate here which of the ranges is the most abundant in each case
	highest_valu_p = 0
	highest_valu_c = 0
	highest_valu_o = 0
	highest_valu_f = 0
	highest_valu_g = 0
	highest_valu_s = 0
	taxa_final = {}
			
	for keys,values in dicc_tax.items():

		letter=list(map(str.strip, str(keys).split('_')))[0]
				
		if(letter == 'k'): #Each range is determined by a letter indicating the taxonomic range
			if((float(values)/float(cont_leaves)) >= taxo_p):
				if(float(values) > highest_valu_k):
					taxa_final[keys]=float(values)
					highest_valu_k = float(values)
		if(letter == 'p'):
			if((float(values)/float(cont_leaves)) >= taxo_p):
				if(float(values) > highest_valu_p):
					taxa_final[keys]=float(values)
					highest_valu_p = float(values)
		if(letter == 'c'):
			if((float(values)/float(cont_leaves)) >= taxo_p):
				if(float(values) > highest_valu_c):
					taxa_final[keys]=float(values)
					highest_valu_c = float(values)
		if(letter == 'o'):
			if((float(values)/float(cont_leaves)) >= taxo_p):
				if(float(values) > highest_valu_o):
					taxa_final[keys]=float(values)
					highest_valu_o = float(values)
		if(letter == 'f'):
			if((float(values)/float(cont_leaves)) >= taxo_p):
				if(float(values) > highest_valu_f):
					taxa_final[keys]=float(values)
					highest_valu_f = float(values)
		if(letter == 'g'):
			if((float(values)/float(cont_leaves)) >= taxo_p):
				if(float(values) > highest_valu_g):
					taxa_final[keys]=float(values)
					highest_valu_g = float(values)
		if(letter == 's'):
			if((float(values)/float(cont_leaves)) >= taxo_p):
				if(float(values) > highest_valu_s):
					taxa_final[keys]=float(values)
					highest_valu_s = float(values)
			
	king = 0			
	phyl = 0
	clas = 0
	ord = 0
	fam = 0			
	gene = 0
			
	array_letters = ['k','p','c','o','f','g','s']
	
	for item_letters in array_letters: #at this point, we evalue, in order from kingdom to species, if each taxonomic range is greater than the minimum percentage stablished. If it is greater, then then most abundante range will be printed. If not, the loop would finish indicating only the taxonomic ranges that passed the stablished threeshold
	
		for keys,values in taxa_final.items():
			
			letter=list(map(str.strip, str(keys).split('_')))[0]
			if(letter == item_letters and (float(values)/float(cont_leaves)) >= taxo_p):
				output_file.write(keys + ';')		
		
	output_file.write('\t')
	return(float(cont_leaves))

##############################################################################################
if __name__ == "__main__":

	main()



