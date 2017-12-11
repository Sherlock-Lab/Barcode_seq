# python bartender_BC1_BC2_lucas.py fastq_f fastq_r multiplex_file

# you need to have the directory bowtie_file/ created with the
# Reference file (fasta file) here Ref_bowtie.fasta, and bowtie will
# create more files.

# the output of this script gives the file "corrected_final.txt" with
# the cluster_id \t the BC (center of the cluster) \t and its count to
# pool all the files per condtion you need to apply pooling_file.py

import os
import sys
import re
import time
import os.path
import subprocess
from subprocess import call

inFile1 = sys.argv[1]  # Forward fastq file
inFile2 = sys.argv[2]  # Reverse fastq file
infile3 = sys.argv[3]  # multiplex file
out_file = inFile1[0:(len(inFile1)-6)] # define output for PEAR (fastq files name wihtout extension)

os.system('mkdir interm_file') # create a directory to put the intermediate file

os.system('pear -f %s -r %s -o interm_file/%s -j 8 -n 150 ' %(inFile1, inFile2, out_file)) # -j threads, Merge F and R Fastq files when overlapping and adjust Qscore

os.system("/Users/lucasherissant/bowtie2-2.2.0/bowtie2-build bowtie_file/Ref_bowtie.fasta bowtie_file/ref_BC -q")  # Builds the new reference db using the ref fasta file												

os.system("/Users/lucasherissant/bowtie2-2.2.0/bowtie2 -U interm_file/%s.assembled.fastq -p 8 -5 17 -3 17 --reorder --no-hd --n-ceil L,0,0.50 -L 12 --np 0 -x bowtie_file/ref_BC -S interm_file/seq.sam --norc" %(out_file))  

# -p 8 (uses 8 threads), 
# -5 17 (cuts first 17, 5' side), 
# -reorder (keeps -in FASTQ order), 
# --no-hd (no 3 line header), 
# --n-ceil (ambiguous characters), 
# -L 12 (length of seed), 
# -x (uses database), 
# -fr -1 -2 (uses 1 and 2 Fastq files in F and R),
# --np sets the penalty to 0 if N in seq 
# --norc does not try to align to reverse complement


handle_sam =  open('interm_file/seq.sam','r+') #Reads the extracted parts of the SAM File and breaks down cig_str to find start & end point of the barcode		
out_cigar = open('interm_file/cigar_break.txt','w')
counter_star = 0
count_diff_size = 0 
count_26 = 0 
while True:
	#Split the SAM file, output of Bowtie2, to break the cigar string and extract the BC
	#	The output here is the seqID BC2 qualBC2  BC1 qualBC1
	read = handle_sam.readline()   #Splits SAM file line by line
	if not read: 
		break

	try:

		read_list = read.rstrip().split('\t')
		
		if len(read) != 1 :

			seq_id,flag,start,cigar_string,sequence,quality = read_list[0],read_list[1],read_list[3],read_list[5],read_list[9],read_list[10]

		pos = 0  # Defines our position in the merged read
		start_aln = int(start) - 1 # Defines the starting point of the alignement 
		
		bar_start_f = 46  #Defines where the first BC start (this the BC2 or high complexity)
		bar_start_f = (bar_start_f - start_aln)  #Defines where the first BC start (this the BC2 or high complexity) if the alignement does not start at 1
		bar_end_f = bar_start_f + 26  # Defines where the first BC end (based on the size of the BC 26 nucleotides)

		bar_start_r = 106  #Defines where the second BC start (this the BC1 or low complexity), named because this BC is at the 3' end of the read so mostly on reverse reads
		bar_start_r = (bar_start_r - start_aln)  #Defines where the second BC start (this the BC1 or low complexity) if the alignement does not start at 1
		bar_end_r = bar_start_r + 26  # Defines where the second BC end (based on the size of the BC 26 nucleotides)

		if cigar_string == '*':  #If cig_str = * , it adds + 1 to counter_* and continues reading the rest of the lines
			counter_star += 1 
			continue

		if cigar_string != '164M' and cigar_string != '163M':
			cigar_op = re.compile("([M,I,D])").split(cigar_string)
			for i in range(0, len(cigar_op) -1, 2):

				if pos < bar_start_f:
					#if we are before the BC

					if cigar_op[i+1] == 'M':																				
						pos = pos + int(cigar_op[i])  #Counts the number of ' M  ' and moves current position positively accordingly
					
					elif cigar_op[i+1] == 'I': 
						bar_start_f = bar_start_f + int(cigar_op[i])  #Counts number of ' I ' and skips current position positively accordingly
						bar_end_f = bar_start_f + 26
						
						bar_start_r = bar_start_r + int(cigar_op[i])  #Counts number of ' I ' and skips current position  accordingly
						bar_end_r = bar_start_r + 26

						pos = pos + int(cigar_op[i])  #moves current position positively 

					elif cigar_op[i+1] == 'D':  #Counts number of ' D ' and moves current position accordingly
						bar_start_f = bar_start_f - int(cigar_op[i])
						bar_end_f = bar_start_f + 26

						bar_start_r = bar_start_r - int(cigar_op[i])  #Counts number of ' D ' and skips current position  accordingly
						bar_end_r = bar_start_r + 26

						pos = pos - int(cigar_op[i])  #moves current position positively

				elif pos >= bar_start_f and pos < bar_end_f:
					#if we are within the first BC do not change the start of the first BC but the end and the second BC (named bar_r)

					if cigar_op[i+1] == 'I':  #Counts number of ' I ' within the forward BC and moves end point and reverse BC positively accordingly
						bar_end_f = (bar_start_f + 26) + int(cigar_op[i])
						bar_start_r = bar_start_r + int(cigar_op[i])  
						bar_end_r = bar_start_r + 26

						pos = pos + int(cigar_op[i])  #moves current position positively

					elif  cigar_op[i+1] == 'D':  #Counts number of ' D ' between position forward BC and moves end point and reverse BC negatively accordingly
						bar_end_f = (bar_start_f + 26) - int(cigar_op[i])
						bar_start_r = bar_start_r - int(cigar_op[i])  #Counts number of ' D ' and skips current position  accordingly
						bar_end_r = bar_start_r + 26

						pos = pos - int(cigar_op[i])  #moves current position positively

					elif cigar_op[i+1] == 'M':																				
						pos_= pos + int(cigar_op[i])  #Counts the number of ' M  ' and moves current position positively accordingly

				elif pos > bar_end_f and pos <= bar_start_r:
					#if we are after the first BC change the point of the second BC

					if cigar_op[i+1] == 'I':  #Counts number of ' I ' within the forward BC and moves end point and reverse BC positively accordingly
						bar_start_r = bar_start_r + int(cigar_op[i])  
						bar_end_r = bar_start_r + 26

					elif  cigar_op[i+1] == 'D':  #Counts number of ' D ' between position forward BC and moves end point and reverse BC negatively accordingly
						bar_start_r = bar_start_r - int(cigar_op[i])  #Counts number of ' D ' and skips current position  accordingly
						bar_end_r = bar_start_r + 26		

					elif cigar_op[i+1] == 'M':																				
						pos = pos + int(cigar_op[i])  #Counts the number of ' M  ' and moves current position positively accordingly

				elif pos >= bar_start_r and pos < bar_end_r:
					#if we are within the second BC change the end point of the second BC

					if cigar_op[i+1] == 'I':  #Counts number of ' I ' within the forward BC and moves end point and reverse BC positively accordingly
						bar_end_r = (bar_start_r + 26) + int(cigar_op[i])

					elif  cigar_op[i+1] == 'D':  #Counts number of ' D ' between position forward BC and moves end point and reverse BC negatively accordingly
						bar_end_r = (bar_start_r + 26) - int(cigar_op[i])

					elif cigar_op[i+1] == 'M':																				
						pos = pos + int(cigar_op[i])  #Counts the number of ' M  ' and moves current position positively accordingly


		f_sequence = sequence[bar_start_f:bar_end_f]
		f_quality = quality[bar_start_f:bar_end_f]
	
		r_sequence = sequence[bar_start_r:bar_end_r]
		r_quality = quality[bar_start_r:bar_end_r]
		
		if len(f_sequence) != 26 or len(r_sequence) != 26:
			count_diff_size += 1
		else:
			count_26 += 1

		w_cigar = seq_id,f_sequence,f_quality,r_sequence,r_quality
		out_cigar.write('\t'.join(w_cigar) + '\n')  #Outputs the seq_id, forward seq, forward qual, reverse seq, and reverse qual 
	except StopIteration:
		break
out_cigar.close()
handle_sam.close()																													

print('number of unaligned read with * = ' + str(counter_star))	#Prints out the counter(number) of cigar strings with ' * ' 
print('number reads with one BC at least with wrong size =' + str(count_diff_size))
print('number reads with both BC at good size = ' + str(count_26))

multiplex_file = open(infile3, 'r') # !! if fastq file is merged file put reverse sequence of reverse primers !
multiplex_dic = {}
outfile_dic = {}
multiplex_pool = {}
nb_line_dic = {}

for line in multiplex_file:

	# here we create directory and file for each multiplex pairs
	# based on the multiplex files that contains: name of sample
	# \t forward multiplex \t Reverse multiplex with (reversed
	# complemented) sequence compare to the primer then pool number but not used anymore.

	# we create the path to a file as a value of a dictionnary that we will call
	# later be carefull because we have to open and close the file to write in it. We use append and
	# by doing so if the file exists it will add at the end.
	# To avoid this we can create, open the file with write (it will create an empty file) and close it immediately.

	(dir_name, multiplex_f, multiplex_r, pool)  = line.rstrip('\n').split('\t') # split the line
	multiplex = multiplex_f + '_' + multiplex_r
	multiplex_dic[multiplex] = dir_name
	nb_line_dic[multiplex] = 0
	os.system('mkdir %s' % dir_name)

	outFile = '%s/table_file.txt' % dir_name

	File = open(outFile, 'w') # allows to create an empty file in case you have to run the same script.
	File.close()

	outfile_dic[multiplex] = outFile


multiplex_file.close()

count_unmatch_id = 0
line_numb = -3
count_line = 0
out_new_multi = open('No_multiplex.txt','w')

# take the merged file to reassociate the UMI and Multiplex

# Reads the extracted barcodes and adds UMI to both ends while keeping
# seq ID, umi, barcode, and qualities in sync

with open('interm_file/%s.assembled.fastq' % out_file,'r+') as fastq, open('interm_file/cigar_break.txt', 'r+') as cigar_break:  

	while True:

		try:
			# In the first step of this loop we
			# reassociate the UMI, multiplex to the
			# barcodes based on the order of the merged
			# files (output of the merger) we make sure we
			# reassociate the UMI and multiplex to the
			# good reads with the SeqID. If seq ID is
			# different it does not write in file.  We
			# have to do so because the bowtie failed to
			# align some reads to the ref and we get rid
			# of those reads so basically if seqID are !=,
			# it breaks the loop and start over and try to
			# associate the next reads in the merged file
			# with the same line in the cigar break file

			table = cigar_break.readline()
			if not table: 
				break																							
			(seq_id, barcode_f, q_bar_f, barcode_r, q_bar_r) = table.strip('\n').split('\t')  #Splits the file containing sequence id, barcode(start/end point), and quality(start/end point)

			while True:  #Defines FASTQ files as sequence ID, Sequence, + Sign, and Sequence Quality
				seq_id_fastq = fastq.readline().strip('\n').split(' ')[0]
				seq_id_fastq = seq_id_fastq[1:]  #Removes the ' @ ' out of the sequence ID

				if seq_id == seq_id_fastq:

					seq_fastq = fastq.readline().strip('\n')
					sign_fastq = fastq.readline().strip('\n')
					quality_fastq = fastq.readline().strip('\n')


					end_read = len(seq_fastq)
					multi_r_end = end_read - 8
					multi_r_start = multi_r_end - 9
				
					line_numb += 4

					umi_f = seq_fastq[:8]  #Extracts out the forward UMI
					mult_f = seq_fastq[8:14]  #Enable this command in order to join the multiplex with UMI and barcode
					umi_quality_f = quality_fastq[:8]  #Extracts out the quality for the forward UMI 
					mult_quality_f = quality_fastq[8:14]  #Enable this command in order to join the multiplex quality with UMI and barcode quality
					umi_r = seq_fastq[multi_r_end:end_read]  #Extracts out the reverse UMI								
					mult_r = seq_fastq[multi_r_start:multi_r_end]  #Enable this command in order to join the multiplex with UMI and barcode
					umi_quality_r = quality_fastq[multi_r_end:end_read]  #Extracts out the quality for the reverse UMI
					mult_quality_r = quality_fastq[multi_r_start:multi_r_end]  #Enable this command in order to join the multiplex quality with UMI and barcode quality

					break

				else:
					count_unmatch_id += 1

			F_R_multi = mult_f + '_' + mult_r
			numb_line = str(line_numb)

			w = str.join('\t',(numb_line,seq_id,umi_f,mult_f,barcode_f,barcode_r,mult_r,umi_r,umi_quality_f,mult_quality_f,q_bar_f,q_bar_r,mult_quality_r,umi_quality_r + '\n'))

			if F_R_multi in outfile_dic: # and len(barcode_r) < 33 and len(barcode_f) < 33:

				# check orientation of multiplex depending on the file we extracted it (if merged file put reverse_comp)

				# in the second part of this loop we associate the reads with its multiplex 

				# we identify the multiplex either directly or by Hamming distance and write the table with all information we got for that read UMI, multiplex, BC ...
				
				sample_file = open(outfile_dic[F_R_multi], 'a')
				count_line += 1
				line_multi = int(nb_line_dic[F_R_multi])
				new_line_multi = line_multi + 1
				nb_line_dic[F_R_multi] = new_line_multi
				

				sample_file.write(w) 
				sample_file.close()

			elif F_R_multi not in outfile_dic:# and len(barcode_r) < 33 and len(barcode_f) < 33:
				new_multi = True

				for j in outfile_dic:
					diffs = 0
					index_mm = 0
					index_mm_lst = []

					for chr1, chr2 in zip(F_R_multi, j):  # for each character in the multiplex
						index_mm += 1
						if chr1 != chr2:
							diffs += 1  # if the 2 character are different add 1 to diffs

					if diffs <= 2:
						sample_file = open(outfile_dic[j], 'a')						
						count_line += 1
						line_multi = int(nb_line_dic[j])
						new_line_multi = line_multi + 1
						nb_line_dic[j] = new_line_multi
						

						sample_file.write(w) 
						sample_file.close()
			
						new_multi = False
						break

				if new_multi == True:
					out_new_multi.write(w)


		except StopIteration:
			break			

	 

fastq.close()
cigar_break.close()
outfile_dic.clear()


print('unmatched ID when mergin BC and UMI Multiplex = ' + str(count_unmatch_id))

print(nb_line_dic)


summary_stat_file = open('summary_stat.txt', 'a')
w_header = str.join('\t',('timepoint','total_umis','uniq_umi','umi_count','total_bc','bc_count','percent_dupl','percent_lost_reads','percent_lost_threshold','percent_lost_BC'+ '\n'))
summary_stat_file.write(w_header)
#file to keep track of different information about the PCR duplicate, lost of BC and reads due to low frequency BC.

not_clust_count = 0
count_no_clust = 0
count_no_BC2clust = 0

for key in multiplex_dic:

	# Loop of each sample/timepoints one by one based on the multiplex with the dictionnary created earlier
	# We created a specific directory for each pool (= each multiplex, or = each sample)
	# we have tab formatted files to keep all the information we want. Now we want to infer cluster of BC1 with bartender and have to make a files with the BC1 only the BC2 only and the UMI only each files contains a second column with the line number

	directory = multiplex_dic[key]

	table_file = open('%s/table_file.txt' % directory,'r+')
	table_BC1 = open('%s/table_BC1_bartender.txt' % directory,'w')
	table_BC2 = open('%s/table_BC2_bartender.txt' % directory,'w')
	table_UMI = open('%s/table_UMI.txt' % directory,'w')

	for line_table in table_file:

		(line_numb, seq_id, umi_f, multi_f, bar_f, bar_r, multi_r, umi_r, umi_qual_f, multi_qual_f, bar_qual_f, bar_qual_r, multi_qual_r, umi_qual_r) = line_table.split('\t')

		if len (bar_f) > 3 : 
			wBC1 = str(bar_r) + ',' + str(line_numb) + '\n'
			table_BC1.write(wBC1)
			wBC2 = str(bar_f) + ',' + str(line_numb) + '\n'
			table_BC2.write(wBC2)
			wUMI = str(umi_f) + str(umi_r) + ',' + str(line_numb) + '\n'
			table_UMI.write(wUMI)


	table_BC2.close()
	table_BC1.close()
	table_file.close()
	table_UMI.close()

	print('++++++++++++++++++++++++++++++' + '\n')
	print(directory)
	print('\n' + '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++' + '\n')

	print(time.clock())
	print(time.asctime(time.localtime(time.time())))

	# run bartender with the command line, no cutoff to not mess up with the line number -c 1, disable statistical test  -z -1, set up hamming distance at 2 -d 2, seed length 8 -l, 8 threads.
	
	os.system('bartender_single_com -f %s/table_BC1_bartender.txt -o %s/clust_BC1 -c 1 -t 8 -s 1 -l 8 -z -1 -d 2' % (directory, directory)) # cutoff_BC1

	os.system('bartender_single_com -f %s/table_BC2_bartender.txt -o %s/clust_BC2 -c 1 -t 8 -s 1 -l 8 -z -1 -d 2' % (directory, directory)) # cutoff_BC2
	
	print('========================' + '\n')
	print(directory)
	print('Clustering with Bartender Done !!')
	print('\n' + '========================' + '\n')
	
	
	group1 = {}

	for line in open("%s/clust_BC1_barcode.csv" %(directory),'r'): # output of bartender gives the uniques BC, b, their frequency, f, and the cluster they belong c_id.
		b,f,c_id = line.strip().split(',')
		if c_id.isdigit():
			group1[b] = c_id
		
	group2 = {}

	for line in open("%s/clust_BC2_barcode.csv" %(directory),'r'):# output of bartender gives the uniques BC, b, their frequency, f, and the cluster they belong c_id.
		b,f,c_id = line.strip().split(',')
		if c_id.isdigit():
			group2[b] = c_id
		
	merged_umis = {}

				
	bc1_handle = open("%s/table_BC1_bartender.txt" % directory,'r')
	bc2_handle = open("%s/table_BC2_bartender.txt" % directory,'r')
	umi_handle = open('%s/table_UMI.txt' % directory,'r')

	for bc1 in bc1_handle:

		# go through the 3 original files at the same time containing the same lines and associate the cluster id and umi to each reads/line
		# make a dict for each cluster based on the cluster id and put all the umi in a list as value
		# for each cluster add every UMI seen, will check them the next loop

		(bc1, bc1_line) = bc1.strip().split(',')
		(bc2, bc2_line) = bc2_handle.readline().strip().split(',')
		(umi, umi_line) = umi_handle.readline().strip().split(',')

		id1 = group1[bc1]
		id2 = group2[bc2]
		merged_id = id1 + '_' + id2

		if merged_id in merged_umis:
			merged_umis[merged_id].append(umi)
		else:
			merged_umis[merged_id] = [umi]

	before_pcr = 0
	after_pcr = 0
	merged_id_umi_dic = {}

	for key_id in merged_umis:
		#loop over the list of UMI for each cluster
		umi_dict = {}
		duplicate = len(merged_umis[key_id])
		before_pcr += duplicate
		for u in merged_umis[key_id]: #loop over all the umis for one cluster_id and keep the unique using dict
			if u in umi_dict:
				umi_dict[u] += 1
			else:
				umi_dict[u] = 1
		merged_id_umi_dic[key_id] = len(umi_dict) 
		after_pcr += len(umi_dict)					

	print(len(merged_umis))

	print(len(merged_id_umi_dic))

	# storing cluster centers for each cluster id, for both BC1 and BC2

	clusterBC1_file = open('%s/clust_BC1_cluster.csv' % directory, 'r')
	cluster_BC1 = {}

	for line_BC1 in clusterBC1_file: #loop over the real cluster for BC1 using the id as key and the cluster center or "real BC" as value
		if c_id.isdigit():
			(cluster_id,BC1_clust,Score,point) = line_BC1.strip().split(',')
			cluster_BC1[cluster_id] = BC1_clust
	clusterBC1_file.close()

	clusterBC2_file = open('%s/clust_BC2_cluster.csv' % directory, 'r')
	cluster_BC2 = {}
	for line_BC2 in clusterBC2_file:  #loop over the real cluster for BC2 using the id as key and the cluster center or "real BC" as value
		if c_id.isdigit():
			(cluster_id_BC2,BC2_clust,Score,point) = line_BC2.strip().split(',')
			cluster_BC2[cluster_id_BC2] = BC2_clust
	clusterBC2_file.close()
	

	count_BC_threshold = 0
	total_count = 0
	output = open('%s/cluster_final.txt' % directory, 'w')

	for key_id in merged_id_umi_dic:
		# reassociate the BC1 and teh BC2 together and keep all the 2BC that have more than 2 different UMI
		if int(merged_id_umi_dic[key_id]) > 0:
			BC1_id, BC2_id = key_id.split('_')

			w_cluster = str(BC2_id) + str('_') + str(BC1_id) + '\t' + str(cluster_BC2[BC2_id]+ '\t' + str(cluster_BC1[BC1_id]) + '\t' + str(merged_id_umi_dic[key_id]) + str('\n'))
			output.write(w_cluster)
			count_BC_threshold += 1
			total_count += int(merged_id_umi_dic[key_id])

	print('total of umis = ')
	print(before_pcr)

	print('total of Unique umis = ')
	print(after_pcr)	
	
	pourcentage_pcr = 1 - after_pcr / before_pcr
	print('pourcentage of pcr duplicate =')
	print(pourcentage_pcr)
	

	print('total count of reads after threshold = ')
	print(total_count)	


	pourcentage_lost_pcr = 1 - total_count / before_pcr
	print('pourcentage of lost reads total (including pcr effect and threshold) =')
	print(pourcentage_lost_pcr)

	pourcentage_lost = 1 - total_count / after_pcr
	print('pourcentage of lost reads after the threshold at >5 =')
	print(pourcentage_lost)

	print('total of BC after threshold = ')
	print(count_BC_threshold)	

	pourcentage_bc = 1- count_BC_threshold / len(merged_id_umi_dic)
	print('percentage of lost BC with threshold at >5 =')
	print(pourcentage_bc)

	w2_data = str(directory) + str('\t') + str(before_pcr) + str('\t') + str(after_pcr) + str('\t') + str(total_count) + str('\t') + str(len(merged_id_umi_dic)) + str('\t') + str(count_BC_threshold) + str('\t') + str(pourcentage_pcr) + str('\t') + str(pourcentage_lost_pcr) + str('\t') + str(pourcentage_lost) + str('\t') + str(pourcentage_bc) + str('\n')

	summary_stat_file.write(w2_data)

