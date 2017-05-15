#2017-05-09
#A. Pendleton
#This script is used to make an AGP file for Zoey. Can be built upon as new forms
#   of contigs become assembled and incorporated into the AGP.

#this uses iPython magic to make plots appear inline
import subprocess
import os
import sys
import numpy as np
import glob

###################################################################################################
def get_chrom_lengths():
	chrom_length_file = open(wkDir + 'input/cyto_band.bed','r')
	for line in chrom_length_file:
		line=line.rstrip().split('\t')
		chrom,length = line[0],int(line[2])
		if 'chrUn' in chrom or 'chrM' in chrom:
			continue
		chrom_lengths[chrom] = length
	chrom_length_file.close()
	return chrom_lengths

###################################################################################################
# STEP TWO
###################################################################################################
def read_line(Dict,curr_index,i):
	#Contig 1
	chrom1,start1,end1,contigID1,Dir1,length1 = Dict[curr_index][0],int(Dict[curr_index][1]),int(Dict[curr_index][2]), Dict[curr_index][3], Dict[curr_index][4], int(Dict[curr_index][5])
	#Contig 2
	chrom2,start2,end2,contigID2,Dir2,length2 = Dict[i+1][0],int(Dict[i+1][1]),int(Dict[i+1][2]), Dict[i+1][3], Dict[i+1][4], int(Dict[i+1][5])
	return chrom1,start1,end1,contigID1,Dir1,length1,chrom2,start2,end2,contigID2,Dir2,length2
###################################################################################################
def read_last_line(Dict):
	#Contig 1
	chrom1,start1,end1,contigID1,Dir1,length1 = Dict[i-1][0],int(Dict[i-1][1]),int(Dict[i-1][2]), Dict[i-1][3], Dict[i-1][4], int(Dict[i-1][5])
	#Contig 2
	chrom2,start2,end2,contigID2,Dir2,length2 = Dict[i][0],int(Dict[i][1]),int(Dict[i][2]), Dict[i][3], Dict[i][4], int(Dict[i][5])
	return chrom1,start1,end1,contigID1,Dir1,length1,chrom2,start2,end2,contigID2,Dir2,length2
###################################################################################################
def process_last_contig(contigDict,curr_index,i):
	chrom1,start1,end1,contigID1,Dir1,length1,chrom2,start2,end2,contigID2,Dir2,length2 = read_last_line(contigDict)
	#print ('\n#%s (%i) -- %s (%i)' % (contigID1,i-1,contigID2,i))
	logFile.write(('\n#%s (%i) -- %s (%i)\n' % (contigID1,i-1,contigID2,i)))
	if start1 == start2 or end1 == end2:
		return
	if end2 > end1:
		curr_index+=1 
		curated_contigDict[i] = [chrom2,start2,end2,contigID2,Dir2,length2]
		inList.append(contigID2)
		#print('PASS - Last contig on chromosome overlap or is spaced correctly')
		logFile.write('PASS - Last contig on chromosome overlap or is spaced correctly\n')
		return
###################################################################################################
def process_first_contig(contigDict,curr_index,i):
	curr_index+=1
	#Always save first contig in the dictionary/dataset
	logFile.write('PASS - First call in dataset\n')
	curated_contigDict[i] = [chrom1,start1,end1,contigID1,Dir1,length1]
	
	#If contig #2 is fully within contig #1 -- do not add contig #2 to the dictionary
	if start1 < start2 and end1 > end2: 
		logFile.write('FAIL -- Contig2 fully within contig1' + '\nKeeping contig: ' + contigID1 +  '\n' + str(contigDict[curr_index]) + '\n' + str(contigDict[i+1]) +'\n')
		curr_index = curr_index - 1
		return curr_index 
	else:
		logFile.write('PASS - first and second contig overlap or are spaced correctly\n')
		curated_contigDict[i+1] = [chrom2,start2,end2,contigID2,Dir2,length2]
		inList.append(contigID2)
	
	return curr_index 
###################################################################################################
def process_same_starts(contigDict,curr_index,i):
	if start1 == start2 and end1 > end2:
		#print('FAIL -- Share start, contig1 longer' + '\n' + str(contigDict[i]) + '\n' + str(contigDict[i+1]))
		logFile.write('FAIL -- Share start, contig1 longer' + '\n' + str(contigDict[i]) + '\n' + str(contigDict[i+1]) + '\n')
		return
	if start1 == start2 and end1 < end2: #overwrite previous position if contig 2 is longer than contig 1
		curated_contigDict[i] = [chrom2,start2,end2,contigID2,Dir2,length2]
		#print('FAIL -- Share start, contig2 longer' + '\n' + str(contigDict[i]) + '\n' + str(contigDict[i+1]))
		logFile.write('FAIL -- Share start, contig2 longer' + '\n' + str(contigDict[i]) + '\n' + str(contigDict[i+1]) + '\n')
		return
###################################################################################################
def process_same_ends(contigDict,curr_index,i):    
	if start1 < start2 and end1 == end2:
		#print('Ends same, kept contig 1 only' + '\n' + str(contigDict[i]) + '\n' + str(contigDict[i+1]))
		logFile.write('Ends same, kept contig 1 only' + '\n' + str(contigDict[i]) + '\n' + str(contigDict[i+1]) + '\n')
		return 
	if start2 < start1 and end1 == end2: 
		print('ERROR - is this contig list sorted??!!')
		logFile.write('ERROR - is this contig list sorted??!!\n\n\n')
		sys.exit(1)
###################################################################################################
def rename_keys(curated_contigDict):	
	count, tempDict = 0, {}
	keyList = curated_contigDict.keys()
	keyList.sort()
	
	for key in keyList:
	#for key in curated_contigDict.sorted(keys):
		count+=1
		tempDict[count] = curated_contigDict[key]
	return tempDict
###################################################################################################

###################################################################################################
# STEP THREE
###################################################################################################
def identify_golden_path(curated_contigDict):
	posDict, offset_BLAT = {}, []
	for i in range(1,len(curated_contigDict)+1): 
		curr_index = i
		offset = 0 #set equal to zero at beginning, will change if there is one from parsing the blat hit(s)
		
		#If last entry, automatically goes in
		if i == len(curated_contigDict):
			chrom1,start1,end1,contigID1,Dir1,length1 = curated_contigDict[curr_index][0],int(curated_contigDict[curr_index][1]),int(curated_contigDict[curr_index][2]), curated_contigDict[curr_index][3], curated_contigDict[curr_index][4], int(curated_contigDict[curr_index][5])
			posDict[i] = [chrom1,start1,end1,contigID1,Dir1,length1,int('0'),int('0'),'na']
			continue
		
		#1. Reading in coordinates from contigList
		chrom1,start1,end1,contigID1,Dir1,length1,chrom2,start2,end2,contigID2,Dir2,length2 = read_line(curated_contigDict, curr_index, i)

		#2. Checking contigs are on same chromosome
		if chrom1 != chrom2:# or i == 1: #This is first contig on chromosome, need to process it first
			continue
	
		#3. Determines the overlap/orientation of the two contigs
		determine_contig_overlap(end1,contigID1,Dir1,length1,end2,contigID2,Dir2,length2)
		
		#4. BLAT the properly oriented contig ends against one another
		#     and parse the results
		if overlap > 0:
			run_blat(wkDir)
			parse_blat(wkDir,overlap,contigID1,Dir1,extract_coord1,contigID2,Dir2,extract_coord2)
		else:
			offset = 0
		if offset > 0:
			offset_BLAT.append([contigID1,Dir1,extract_coord1,contigID2,Dir2,extract_coord2,overlap])
		
		#save contig1 information to dictionary
		posDict[curr_index] = [chrom1,start1,end1,contigID1,Dir1,length1,overlap,offset,contigID2]

		#print('###################\n')
		#if i > 2:
		#    break
	return posDict
###################################################################################################
def process_first_contig_on_chrom(chrom1,start1,end1,contigID1,Dir1,length1):
	chrom2 = chrom1
	start2 = start1 - 10000
	end2 = start1 + 10000
	contigID2 = 'canFam3'
	Dir2 = 'fwd'
	length1 = end2-start2

	if start2 < 0:
		start2 = 0

	#Determine overlap 
	determine_contig_overlap(end1, contigID1, Dir1, length1, end2, contigID2, Dir2, length2)

	#Find coordinates to extract
	find_overlapping_coordinates(overlap,end1, contigID1, Dir1, length1, end2, contigID2, Dir2, length2)

	#Extract FASTA
	fastaRoot = '/home/ampend/links/kidd-lab/jmkidd-projects/zoey/contig-assignment/kmer-matches/eval1/'
	extract_fasta(fastaRoot,contigID1,Dir1,extract_coord1,contigID2,Dir2,extract_coord2)
###################################################################################################
def determine_contig_overlap(end1, contigID1, Dir1, length1, end2, contigID2, Dir2, length2):
	#How much do the contigs overlap?
	global overlap
	overlap = end1 - start2
	logFile.write('Overlap = %i\n\n' % overlap)
	#OVERLAPPING CONTIGS
	if overlap > 0:
		find_overlapping_coordinates(overlap,end1, contigID1, Dir1, length1, end2, contigID2, Dir2, length2)

	#DIRECTLY ADJACENT CONTIGS
	if overlap == 0:
		#print('Contigs are directly adjacent to one another')
		logFile.write('Contigs are directly adjacent to one another\n')
		#WRITE FUNCTION FOR THESE
	#CONTIGS WITH GAP BETWEEN THEM
	if overlap < 0:
		#print('Gap between contigs')
		logFile.write('Gap between contigs\n')
###################################################################################################
def find_overlapping_coordinates(overlap,end1, contigID1, Dir1, length1, end2, contigID2, Dir2, length2):
	###Determine coordinates to extract for BLAT
	global extract_coord1
	global extract_coord2
	#Contig #1
	if 'fwd' in Dir1:
		extract_coord1 = [length1 - overlap - 10000, length1]
	else:
		extract_coord1 = [0, overlap + 10000]
	#Contig #2
	if 'fwd' in Dir2:
		extract_coord2 = [0, overlap + 10000] 
	else:
		extract_coord2 = [length2 - overlap - 10000, length2]

	#print ('Coordinates to extract for BLAT: ',extract_coord1,extract_coord2,'\n')
	logFile.write('Coordinates to extract for BLAT: ' + str(extract_coord1) + str(extract_coord2) + '\n')
	print('Coordinates to extract for BLAT: ' + str(extract_coord1) + str(extract_coord2) + '\n')

	#safety check
	for i in range(0,2): #CHECKS 
		#If the region to extract extends beyond the length of the contig, if so then the coordinate changes to the length of contig
		if extract_coord1[i] > length1:
			extract_coord1[i] = length1
			#print('Fixed coord1 ' + '\t' + str(i) + '\t' + str(length1))	
			logFile.write('Fixed coord1: i = %s    length1 = %s\n' % (str(i), str(length1)))	
		if extract_coord2[i] > length2:
			extract_coord2[i] = length2
			#print('Fixed coord2 ' + '\t' + str(i) + '\t' + str(length2))	
			logFile.write('Fixed coord2: i = %s    length2 = %s\n' % (str(i), str(length1)))	
		# If the region extended too far, and the value is negative -- then make starting extraction coordinate = 0
		if extract_coord1[i] < 0: 
			extract_coord1[i] = 0
			#print('Fixed coord1 ' '\t' + str(i) + '\t0')	
			logFile.write('Fixed coord1: i = %s \t start = 0\n' % (str(i)))	
		if extract_coord2[i] < 0:
			extract_coord2[i] = 0
			#print('Fixed coord2 ' + '\t' + str(i) + '\t0')	
			logFile.write('Fixed coord2: i = %s \t start = 0\n' % (str(i)))	
	
	logFile.write('FIXED coordinates to extract for BLAT: ' + str(extract_coord1) + '\t' + str(extract_coord2) + '\n')
	print('FIXED coordinates to extract for BLAT:\t' + str(extract_coord1) + '\t' + str(extract_coord2) + '\n')
	
	if overlap < 0: #### SEND TO DIFFERENT FUNCTION LATER -- FOLLOW OVERLAPPING CONTIGS FOR NOW
		#print('Contigs do not overlap')
		logFile.write('Contigs do not overlap\n')
	#Where to find fasta of contigs
	fastaRoot = '/home/ampend/links/kidd-lab/jmkidd-projects/zoey/contig-assignment/kmer-matches/eval1/'
	extract_fasta(fastaRoot,contigID1,Dir1,extract_coord1, contigID2,Dir2,extract_coord2) 
###################################################################################################
def extract_fasta(fastaRoot,contigID1,Dir1,extract_coord1,contigID2,Dir2,extract_coord2):
	#Contig 1
	fasta_path1 = fastaRoot + contigID1 + '/' + contigID1 + '.fa'
	cmd = 'samtools faidx %s %s:%i-%i  > %stemp/contig1.fa' % (fasta_path1,contigID1,extract_coord1[0],extract_coord1[1],wkDir)
	#print(cmd)
	logFile.write(cmd + '\n')
	subprocess.call(cmd,shell=True)

	#Contig 2
	fasta_path2 = fastaRoot + contigID2 + '/' + contigID2 + '.fa'
	cmd = 'samtools faidx %s %s:%i-%i > %stemp/contig2.fa' % (fasta_path2,contigID2,extract_coord2[0],extract_coord2[1],wkDir)
	#print(cmd)
	logFile.write(cmd + '\n')
	subprocess.call(cmd,shell=True)   
###################################################################################################
def run_blat(wkDir):
	blatcmd = 'blat %stemp/contig1.fa %stemp/contig2.fa %stemp/temp.blat' % (wkDir,wkDir,wkDir)
	print(blatcmd)
	logFile.write(blatcmd + '\n')
	subprocess.call(blatcmd,shell=True)
###################################################################################################
def parse_blat(wkDir,overlap,contigID1,Dir1,extract_coord1,contigID2,Dir2,extract_coord2):
	inFile = open('%stemp/temp.blat' % (wkDir),'r')
	lineCount, cleanAlignment, offset = 0, False, 0
	matchHits = []
	
	for line in inFile:
		line=line.rstrip().split('\t')
		lineCount+=1
		if lineCount > 5:#skips the BLAT results headers
			score,strand = int(line[0]),line[8]
			percent_of_overlap = float(score)/overlap
			if percent_of_overlap < float(0.75):
				continue
			#print('\n#Parsing BLAT results','Percent of overlap matched in BLAT: ',format(percent_of_overlap, '.3f'))
			logFile.write('\n#Parsing BLAT results\n' + 'Percent of overlap matched in BLAT: ' + format(percent_of_overlap, '.3f') + ' \n')
			#CHECKS TO MAKE SURE WE DID THIS RIGHT - correct end vs correct end
			#if '+' not in strand:
			#    print('ERROR: Top hit not in proper orientation.... SKIPPING -- PLEASE CHECK')
			#    return
			#parse contig #1 (left contig)
			blat_length1, blat_start1, blat_end1 = int(line[14]),int(line[15]),int(line[16])
			#print(blat_length1, blat_start1, blat_end1)
			logFile.write('length = %i, start = %i, end = %i' % (blat_length1, blat_start1, blat_end1))
			#parse contig #2 (right contig)
			blat_length2, blat_start2, blat_end2 = int(line[10]),int(line[11]),int(line[12])
			#print(blat_length2, blat_start2, blat_end2)
			logFile.write('length = %i, start = %i, end = %i' % (blat_length2, blat_start2, blat_end2))
			#Calculate offset  - to compensate for non-clean alignments (extension of left contig that does not overlap with adjacent contig)
			offset = calculate_offset(overlap,blat_length1, blat_start1, blat_end1,Dir1,blat_length2, blat_start2, blat_end2,Dir2)
			cleanAlignment = True
			return offset, cleanAlignment
	if cleanAlignment is False:
		return offset, cleanAlignment
###################################################################################################
def calculate_offset(overlap,blat_length1, blat_start1, blat_end1,Dir1,blat_length2, blat_start2, blat_end2,Dir2):
	global offset
	if 'rc' in Dir1:
		offset = blat_start1 - 0
	else:
		offset = blat_length1 - blat_end1
	if offset < 3:
		offset = 0
	#print('\nBLAT offset = ', offset)
	logFile.write('\nBLAT offset = ' + str(offset) + '\n')
	return offset
###################################################################################################

###################################################################################################
# STEP FOUR
###################################################################################################
def process_first_agp_contig(agp_chrom, start, end, contigID, direction, length, overlap, offset, pairedContig):
	agp_prev_chrom, agp_prev_overlap, prev_offset = agp_chrom, 0, 0

	agp_start, agp_end = 1, length
	if direction == 'fwd':
		direction = '+' #changes notation of the direction
	if direction == 'rc':
		direction = '-'  #changes notation of the direction

	#Specially process those with offsets > 0
	if offset > 0 or prev_offset > 0:
		info = process_blat_offsets(agp_prev_chrom, agp_prev_overlap, prev_offset, agp_start, agp_end, overlap, offset,direction)
	#Those with BLAT offsets = 0
	else:
		agp_start, agp_end = 1, length
		contig_start, contig_end = 1, length
		info = [agp_prev_chrom,agp_prev_overlap,agp_start, agp_end,contig_start,contig_end,direction,offset]

	return info
###################################################################################################
def process_blat_offsets(agp_prev_chrom, agp_prev_overlap, prev_offset, agp_start, agp_end, overlap, offset,direction):
	if agp_prev_overlap < 5 and agp_prev_overlap > 0:
		agp_prev_overlap = 0
	
	if prev_offset ==0  and offset > 0: #contig left of offset
		agp_start = agp_end + 1
		agp_end = agp_start + length + agp_prev_overlap

		#Determine contig coordinates that align
		if direction == 'fwd':
			direction = '+' #changes notation of the direction
			contig_start = 1 + agp_prev_overlap
			contig_end = length
		if direction == 'rc':
			direction = '-'  #changes notation of the direction
			contig_start = 1
			contig_end = length - agp_prev_overlap

	if prev_offset > 0:# and offset == 0: #contig right of offset
		agp_start = agp_end + 500 + 1 #compensate for added gap between the two which has length = 500
		agp_end = agp_start + length + agp_prev_overlap

		if direction == 'fwd':
			direction = '+' #changes notation of the direction
			contig_start = 1 + agp_prev_overlap
			contig_end = length
		if direction == 'rc':
			direction = '-'  #changes notation of the direction
			contig_start = 1
			contig_end = length - agp_prev_overlap  
		
	#if prev_offset > 0 and offset > 0:
	#	print('This type of alignment needs to be readdressed...\n')
	#	sys.exit()
	
	info = [agp_prev_chrom,agp_prev_overlap,agp_start, agp_end,contig_start,contig_end,direction,offset]

	return info
###################################################################################################
def process_next_agp_contig(agp_chrom, direction, length, overlap, offset, pairedContig, agp_start, agp_end):
	agp_prev_chrom, agp_prev_overlap, prev_offset = posDict[i-1][0], int(posDict[i-1][6]),int(posDict[i-1][7])

	#Specially address those with BLAT offsets (unaligned sequence at the junctions)
	if prev_offset > 0 or offset > 0:
		info = process_blat_offsets(agp_prev_chrom, agp_prev_overlap, prev_offset, agp_start, agp_end, overlap, offset,direction)
		return info

	#Junctions without unaligned sequence(s):
	else:
		if agp_prev_overlap < 5:
			agp_prev_overlap = 0

		#Determine contig coordinates that align
		if direction == 'fwd':
			direction = '+' #changes notation of the direction
			contig_start = 1 + agp_prev_overlap
			contig_end = length
		if direction == 'rc':
			direction = '-'  #changes notation of the direction
			contig_start = 1
			contig_end = length - agp_prev_overlap
	
		shift = contig_end - contig_start

		#now calculate the AGP coordinates of where the contig goes
		agp_start = agp_end + 1
		agp_end = agp_start + shift #+ 1
	
		info = [agp_prev_chrom,agp_prev_overlap,agp_start, agp_end,contig_start,contig_end,direction,offset]
		return info
###################################################################################################
def process_agp_GAP(agp_chrom, direction, length, overlap, offset, pairedContig, agp_start, agp_end):
    gap_length = 10000
    agp_start = agp_end + 1
    agp_end = agp_start + gap_length - 1
    if 'CTG' in contigID:
        gap_type = 'contig'
    else:#To change once we add more types of contigs in here
        gap_type = 'OTHER'
    agp_prev_overlap = 0
    
    info = [agp_chrom,agp_start,agp_end,gap_type,gap_length,agp_prev_overlap]
    return info  
###################################################################################################
def write_AGP_header(agpFile): #alter for each assembly
    agpFile.write('##agp-version 1.0\n# ORGANISM: Canis lupus familiaris\n# TAX_ID: 9615\n')
    agpFile.write('# ASSEMBLY NAME: Zoey_v1\n# ASSEMBLY DATE: 19-April-2017\n')
    agpFile.write('# GENOME CENTER: University of Michigan - J.M. Kidd Lab\n')
    agpFile.write('# DESCRIPTION: AGP specifying the assembly of chromosome 18 from primary PacBio contigs from FALCON assembly\n')      
###################################################################################################

###################################################################################################
# STEP FIVE
###################################################################################################
def extract_contig_fasta(fastaRoot,contigID,orient,extract_coord1,extract_coord2):
    #Contig 1
    fasta_path1 = fastaRoot + contigID + '/' + contigID + '.fa'
    cmd = 'samtools faidx %s %s:%i-%i  > %stemp/contig1.fa' % (fasta_path1,contigID,extract_coord1,extract_coord2,wkDir)
    #print(cmd)
    subprocess.call(cmd,shell=True)
###################################################################################################
def extract_filledGap_fasta(fastaRoot,contigID,orient,extract_coord1,extract_coord2):
    #Contig 1
    fasta_path1 = gap_fastaRoot + '/' + contigID + '.fa'
    cmd = 'samtools faidx %s %s:%i-%i  > %stemp/contig1.fa' % (fasta_path1,contigID,extract_coord1,extract_coord2,wkDir)
    #print(cmd)
    subprocess.call(cmd,shell=True)
###################################################################################################    
def reverse_comp_contig_fasta(wkDir):
    cmd = 'fastarevcomp %stemp/contig1.fa > %stemp/contig1.fa.rc' % (wkDir,wkDir)
    #print(cmd)
    subprocess.call(cmd,shell=True)
    
    cmd = 'mv %stemp/contig1.fa.rc %stemp/seq.fa'  % (wkDir,wkDir)
    #print(cmd)
    subprocess.call(cmd,shell=True) 
###################################################################################################    
def read_in_fasta(wkDir):    
    inFile = open(wkDir + 'temp/seq.fa', 'r')
    for line in inFile:
        line=line.rstrip()
        if '>' in line:
            continue
        else:
            seq.append(line)
    inFile.close()
    return seq
###################################################################################################    
def reformat_fasta_file(wkDir,CHROM):
    agp_Fasta = open(wkDir + 'results/AGP_v1.0/AGP_Zoey_Assembly_v1_' + CHROM + '.fa', 'r')
    tmpFile = open(wkDir + 'results/AGP_v1.0/tmp.fa', 'w')

    for line in agp_Fasta:
        if '>' in line:
            tmpFile.write(line)#Writes chrom ID
            continue
        seq=line.rstrip()
        
        start = 0
        fasta_length = float(len(seq))
        fasta_line_length = 80
        Max = int(fasta_length/fasta_line_length)

        for i in range(0,Max+1):
            Line = seq[start:start+fasta_line_length]
            start = start + fasta_line_length
            tmpFile.write(Line + '\n')
            """if i > 10:
                break"""    

    agp_Fasta.close()
    tmpFile.close()
    
    cmd = 'mv %sresults/AGP_v1.0/tmp.fa %sresults/AGP_v1.0/AGP_Zoey_Assembly_v1_%s.fa'  % (wkDir,wkDir,CHROM)
    subprocess.call(cmd,shell=True) 
###################################################################################################    

###################################################################################################
# STEP SEVEN - CONTINUITY
###################################################################################################
def calculate_canFam_continuity(CHROM):
	gapFile = open(wkDir + 'input/gaps.bed.sorted','r')
	#to keep track of the previous call
	prev_chrom, prev_start, prev_end = '', '', ''
	#save continuous stretches to array
	cont_pos_list, per_chrom, count =  [], {}, 0 


	for line in gapFile:
		line=line.rstrip().split('\t')
		chrom,start,end=line[0],int(line[1]),int(line[2])
		if chrom != CHROM:
			continue
		#if 'chrUn' in chrom:
		#	continue
		count+=1

		#If first line in file
		if count == 1:
			per_chrom[chrom] = []
			if start == 0: #gap is at the beginning of the chromosome
				#Save for next
				cont_start = end
				prev_chrom, prev_start, prev_end = chrom, start, end
				continue
			#first continuous stretch on chromosome
			cont_start = 0
			cont_end = start - 1
			cont_pos_list.append([count-1,cont_start,cont_end])
			per_chrom[chrom].append([count-1,cont_start,cont_end])
		
			#next continuous stretch, define start
			cont_start = end
			#Save for next
			prev_chrom, prev_start, prev_end = chrom, start, end
			continue
		#If the next call is on a different chromosome, must save the previous gap's info
		if chrom not in prev_chrom:
			per_chrom[chrom] = []
			#previous continuous stretch
			cont_end = chrom_lengths[prev_chrom]
			cont_pos_list.append([count-1,cont_start,cont_end])
			per_chrom[chrom].append([count-1,cont_start,cont_end])
			#first continuous stretch on new chromosome
			cont_start = 0
			cont_end = start - 1
			#next adjacent stretch, define start
			cont_start = end      
		else:
			cont_end = start - 1
			cont_pos_list.append([count-1,cont_start,cont_end])
			per_chrom[chrom].append([count-1,cont_start,cont_end])
			#next adjacent strech, define start
			cont_start = end
		
		prev_chrom, prev_start, prev_end = chrom, start, end
	
	#print('There are %i gaps in the canFam3 assembly (not including chrUn)' % count)
	gapFile.close()
	

	#Calculating the continuity statistics
	continuous_contig_array = []
	for i in range(0,len(cont_pos_list)):
		start,end = int(cont_pos_list[i][1]),int(cont_pos_list[i][2])
		length = end-start
		continuous_contig_array.append(length)    
	Average_continuity = format(np.mean(continuous_contig_array), '.2f')
	#print('Average continuity (bp): ', Average_continuity)  

	return count, Average_continuity
###################################################################################################

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################





######################################################################################
#Defines the working directory.
wkDir = '/home/ampend/links/kidd-lab/ampend-projects/Zoey_Genome_Project/AGP/'
print('Current working directory is:\n', wkDir)


######################################################################################
#Saves the names of chromosomes so that we are processing each chromosome in a loop
#	to generate Zoey AGP version 1.0

chrom_lengths = {}
chrom_lengths = get_chrom_lengths()

######################################################################################

for i in chrom_lengths:
	CHROM = i #Chromosome being processed
	#if 'chr18' not in CHROM:
	#	continue
	#Each step will save the results and processes to a logfile also defined below.
	logFile = open(wkDir + 'results/AGP_v1.0/logs/ZoeyAGP_LogFile_' + CHROM + '_v1.txt','w')
	logFile.write('\n\n#########Processing %s\n' % CHROM)
	#Loci to investigate further, where they have been collapsed, but we should take a closer
	#	look at the relationship of these. Especially needed to check with the canu assemblies
	lociFile = open(wkDir + 'results/AGP_v1.0/check/Zoey_AGP_CheckLoci_' + CHROM + '_v1.txt','w')
	lociFile.write('\n\n#########Processing %s\n' % CHROM)	
	print('\n\n#########Processing %s...' % CHROM)
	
	####################################################
	#STEP ONE - EXTRACTING CONTIGS ON THIS CHROMOSOME
	####################################################
	
	#Reads in primary contig alignments (no canu processing)
	primaryContigCoord = open(wkDir + 'input/primary.2017-04-10.txt', 'r')

	contigDict,index = {}, 0

	for line in primaryContigCoord:
		if line.startswith("#") is True: #skips header line
			continue
		line = line.rstrip().split('\t')
		contigID, length, Dir, chrom, start, end = line[0],line[1],line[2],line[3],int(line[4]),int(line[5])    
		if CHROM not in chrom:  
			continue
		index += 1
		contigDict[index] = [chrom,start,end,contigID,Dir,length]
		#if index > 15: #FOR TESTING--- REMOVE LATER
		#    break
	
	primaryContigCoord.close()
	print ('Identified coordinates for %i primary contigs on %s' % (len(contigDict), CHROM))
	logFile.write('\nIdentified coordinates for %i primary contigs on %s' % (len(contigDict), CHROM))
	logFile.write('\ncontigDict for %s\n' % CHROM)
	for i in contigDict:
		logFile.write(str(contigDict[i]) + '\n')
	
	print('STEP ONE: DONE!')
	
		
	####################################################
	#STEP TWO - REMOVE UNWANTED CONTIGS FOR AGP 
	####################################################
	curated_contigDict, curr_index, inList = {}, 1, []
	chrom1,chrom2 = '0', '0' #initializing values

	for key in range(1,len(contigDict)+1):
		i = int(key)
		#1. If last contig in dataset 
		if i == len(contigDict):
			continue
	
		#2. If contigs are on different chromosomes:
		if chrom1 != chrom2:
			print('ERROR - this should not happen\n')
			#process_last_contig(contigDict,curr_index,i)
			sys.exit()
		
		chrom1,start1,end1,contigID1,Dir1,length1,chrom2,start2,end2,contigID2,Dir2,length2 = read_line(contigDict,curr_index,i)
		logFile.write('\n#%s (%i) -- %s (%i)\n' % (contigID1,curr_index,contigID2,i+1))

		#3. Automatically saves FIRST CONTIG in dataset
		if curr_index == 1 and len(curated_contigDict) == 0:
			curr_index = process_first_contig(contigDict,curr_index,i)
			continue
	
		#4. If contigs have the same start site, choose the longest contig
		if start1 == start2:
			process_same_starts(contigDict,curr_index,i)
			continue

		#5. If contigs have same end coordinate, choose longest contig (i.e. contig #1 if sorted)
		if end1 == end2: 
			process_same_ends(contigDict,curr_index,i)
			continue
		
		#6. If contigs #1 and #2 have same start AND end, keep contig #1
		if start1 == start2 and end1 == end2:
			logFile.write('FAIL -- Same start and end coordinate' + '\n' + str(contigDict[i]) + '\n' + str(contigDict[i+1]) + '\n')
			continue #because contig1 is already in the dictionary
		
		#7. If contig #2 is fully within contig #1 -- continue
		if start1 < start2 and end1 > end2: 
			logFile.write('FAIL -- Contig2 fully within contig1' + '\nKeeping contig: ' + contigID1 +  '\n' + str(contigDict[i]) + '\n' + str(contigDict[i+1]) +'\n')
			continue
		if start1 > start2 and end1 < end2:
			logFile.write('ERROR - is this contig list sorted??!!\n\n\n')
			sys.exit(1)
	
		#8. If contig passes all these - then automatically saves
		curr_index=i+1
		curated_contigDict[i] = [chrom2,start2,end2,contigID2,Dir2,length2]
		logFile.write('PASS - contigs overlap or are spaced correctly\n')

	#print ('Identified CURATED coordinates for %i primary contigs' % len(curated_contigDict))
	logFile.write('\n##Identified CURATED coordinates for %i primary contigs\n\n########\n\n' % len(curated_contigDict))

	#Re-naming dictionary keys by creating a list of keys and then sorting them
	tempDict = rename_keys(curated_contigDict)
	curated_contigDict = tempDict
	
	logFile.write('\ncurated_contigDict for %s\n' % CHROM)
	for i in curated_contigDict:
		logFile.write(str(curated_contigDict[i]) + '\n')
	print('STEP TWO: DONE!')
	#break
	####################################################
	#STEP THREE - BLAT CURATED CONTIG SET TO IDENTIFY GOLDEN PATH 
	####################################################

	#posDict = identify_golden_path(curated_contigDict)
	posDict, offset_BLAT = {}, []
	for i in range(1,len(curated_contigDict)+1): 
		curr_index = i
		offset = 0 #set equal to zero at beginning, will change if there is one from parsing the blat hit(s)
		
		#If last entry, automatically goes in
		if i == len(curated_contigDict):
			chrom1,start1,end1,contigID1,Dir1,length1 = curated_contigDict[curr_index][0],int(curated_contigDict[curr_index][1]),int(curated_contigDict[curr_index][2]), curated_contigDict[curr_index][3], curated_contigDict[curr_index][4], int(curated_contigDict[curr_index][5])
			posDict[i] = [chrom1,start1,end1,contigID1,Dir1,length1,int('0'),int('0'),'na']
			continue
		
		#1. Reading in coordinates from contigList
		chrom1,start1,end1,contigID1,Dir1,length1,chrom2,start2,end2,contigID2,Dir2,length2 = read_line(curated_contigDict, curr_index, i)
		logFile.write('\n##   Contig %s    ---    Contig %s\n\n' % (contigID1, contigID2))
		
		#2. Checking contigs are on same chromosome
		if chrom1 != chrom2:# or i == 1: #This is first contig on chromosome, need to process it first
			continue
	
		#3. Determines the overlap/orientation of the two contigs
		determine_contig_overlap(end1,contigID1,Dir1,length1,end2,contigID2,Dir2,length2)
		
		#4. BLAT the properly oriented contig ends against one another and parse the results
		if overlap > 0:
			run_blat(wkDir)
			offset, cleanAlignment = parse_blat(wkDir,overlap,contigID1,Dir1,extract_coord1,contigID2,Dir2,extract_coord2)
		else:
			offset, cleanAlignment = 0, True
			
		#5. More closely examine the contigs that did not have a strong BLAT hit
		#process_poor_BLAT_hits()
		if cleanAlignment is False and overlap > 0:
			if curr_index + 1 <= len(curated_contigDict):
				offset = 1000000 # DEFAULT HUGE NUMBER -- NEED TO PUT GAPS UPSTREAM OF THIS CONTIG
				
		#Keep track of those that have BLAT based offsets
		if offset > 0:
			offset_BLAT.append([contigID1,Dir1,extract_coord1,contigID2,Dir2,extract_coord2,overlap,offset])

		#save contig1 information to dictionary
		posDict[curr_index] = [chrom1,start1,end1,contigID1,Dir1,length1,overlap,offset,contigID2]
	
	logFile.write('\nposDict for %s\n' % CHROM)
	for i in posDict:
		logFile.write(str(posDict[i]) + '\n')
	
	print('%i elements in position dictionary from Zoey PacBio contigs on %s....' % (len(posDict),CHROM))

	print('STEP THREE: DONE!')
	
	
	####################################################
	#STEP FOUR - WRITE AGP FILE
	####################################################

	#Defining AGP outfile
	agpFile = open(wkDir + 'results/AGP_v1.0/AGP_Zoey_Assembly_v1_' + CHROM + '.txt', 'w')

	#Now putting them in the coordinates of the Zoey genome!! :) 
	count, agp = 0, [] 

	for i in posDict:
		count += 1
		agp_chrom, start, end, contigID, direction, length, overlap, offset, pairedContig = posDict[i][0:9]
		#print('\n#',contigID)
	
		#Processing first contig in dataset -OR- on a new chromosome
		if count == 1 or agp_chrom != agp_prev_chrom: 
			count = 1
			info = process_first_agp_contig(agp_chrom, start, end, contigID, direction, length, overlap, offset, pairedContig)
			agp_prev_chrom,agp_prev_overlap,agp_start, agp_end, contig_start,contig_end,direction, offset = info[0:10]
			agp.append([agp_prev_chrom, agp_start, agp_end, count, 'D', contigID, contig_start, contig_end, direction ])
			#print(agp_prev_chrom,agp_start, agp_end, contig_start,contig_end,direction,agp_prev_overlap)
	
		#Process next contig (not first contig on chrom)
		else:
			if overlap >= 0: #if contigs overlap --> NO GAP!
				info = process_next_agp_contig(agp_chrom, direction, length, overlap, offset, pairedContig, agp_start, agp_end)
				agp_prev_chrom,agp_prev_overlap,agp_start, agp_end, contig_start,contig_end,direction,offset = info[0:10]
				if offset > 0:
					#Add contig to the left of the gap
					agp.append([agp_chrom, agp_start, agp_end, count, 'D', contigID, contig_start, contig_end, direction])
					#print(agp_prev_chrom,agp_start, agp_end, contig_start,contig_end,direction,agp_prev_overlap)
					#Add introduced gap from unaligned sequences from contig-contig junction BLAT
					count += 1
					agp.append([agp_chrom, agp_end+1, agp_end+499, count, 'U', '500','BLAT_gap','no','na'])
					#print(agp_chrom, agp_end+1, agp_end+501, count, 'U', '500','BLAT_gap','no','na')
				else:
					agp.append([agp_chrom, agp_start, agp_end, count, 'D', contigID, contig_start, contig_end, direction])
					#print(agp_prev_chrom,agp_prev_overlap,agp_start, agp_end, contig_start,contig_end,direction)
			else: #if contigs do not overlap --> GAP!
				#process contig to the left of the gap
				info = process_next_agp_contig(agp_chrom, direction, length, overlap, offset, pairedContig, agp_start, agp_end)
				agp_prev_chrom,agp_prev_overlap,agp_start, agp_end, contig_start,contig_end,direction,offset = info[0:10]
				agp.append([agp_chrom, agp_start, agp_end, count, 'D', contigID, contig_start, contig_end, direction])
				#print(agp_chrom, agp_start, agp_end, count, 'D', contigID, contig_start, contig_end, direction)
			
				#processing gap
				count += 1
				gap_info = process_agp_GAP(agp_chrom, direction, length, overlap, offset, pairedContig, agp_start, agp_end)
				agp_chrom,agp_start,agp_end,gap_type,gap_length,agp_prev_overlap = gap_info[0:7]
				agp.append([agp_chrom,agp_start,agp_end,count,'U',gap_length,gap_type,'no','na'])
				#print(agp_chrom,agp_start,agp_end,count,'U',gap_length,gap_type,'no','na')

	#Write AGP header lines
	write_AGP_header(agpFile)
	
	logFile.write('AGP file for %s\n' % CHROM)
	#Write out to agp File the information for contigs/gaps
	for i in range(0,len(agp)):
		agpFile.write("\t".join(map(str,agp[i])) + '\n')
		logFile.write("\t".join(map(str,agp[i])) + '\n')
	agpFile.close()
	logFile.close()

	print('STEP FOUR: DONE!')

	####################################################
	#STEP FIVE - MAKE AGP FASTA FILE FROM AGP TXT FILE FROM STEP FOUR
	####################################################

	agpFile = open(wkDir + 'results/AGP_v1.0/AGP_Zoey_Assembly_v1_' + CHROM + '.txt', 'r')
	fastaRoot = '/home/ampend/links/kidd-lab/jmkidd-projects/zoey/contig-assignment/kmer-matches/eval1/'
	agp_Fasta = open(wkDir + 'results/AGP_v1.0/AGP_Zoey_Assembly_v1_' + CHROM + '.fa', 'w')

	gap_fastaRoot = wkDir + 'input/filled_gaps/' + CHROM + '/'

	contigID, count, prev_chrom = '', 0, ''

	for line in agpFile:
		if '#' in line:
			continue
		count+= 1
	
		line = line.rstrip().split('\t')
		contigID = line[5]
	
		seq = []
	
		if 'CTG' in contigID or 'chr' in contigID:
			chrom, contigID, extract_coord1, extract_coord2, orient = line[0],line[5],int(line[6]),int(line[7]), line[8]
			if 'CTG' in contigID:
				#Extract fasta
				extract_contig_fasta(fastaRoot,contigID,orient,extract_coord1,extract_coord2)
			if 'chr' in contigID:
				extract_filledGap_fasta(gap_fastaRoot,contigID,orient,extract_coord1,extract_coord2)
			#Reverse complement extracted FASTA if needed
			if '-' == orient:
				reverse_comp_contig_fasta(wkDir)
			else:
				cmd = 'mv %stemp/contig1.fa %stemp/seq.fa'  % (wkDir,wkDir)
				#print(cmd)
				subprocess.call(cmd,shell=True) 
		
			seq = read_in_fasta(wkDir)
		else: #These are gaps   
			gap_length = int(line[5])
			seq = 'N' * gap_length
	
		#Writes '>' chromosome to the fasta file
		if count == 1 or prev_chrom != chrom:
			if count > 1:
				agp_Fasta.write('\n')
			agp_Fasta.write('>' + chrom + '\n')
	
		for i in range(0,len(seq)):
			Seq = seq[i].upper()#.upper()print(seq[i])
			agp_Fasta.write(Seq)

		prev_chrom = chrom #reset prev_chrom identity


	print('Finished writing out the sequence for %i contigs and gaps to Zoey fasta file' % count)
	agpFile.close()
	agp_Fasta.close()

	#Reformat the fasta file to only have certain # of nucleotides per line
	reformat_fasta_file(wkDir,CHROM)
	print('STEP FIVE: DONE!')
	
	####################################################
	#STEP SIX - CALCULATE CHROM CONTINUITY FOR ZOEY
	####################################################

	agpFile = open(wkDir + 'results/AGP_v1.0/AGP_Zoey_Assembly_v1_' + CHROM + '.txt', 'r')

	contig_count,gap_count,line_count,blat_gap_count = 0,0,0,0
	prev_Type,start_stop_list,contig_start = '',[],1

	for line in agpFile:
		if '#' in line: #skips header
			continue    
		line = line.rstrip().split('\t')
		line_count += 1
		chrom,start,end,num,Type = line[0], int(line[1]),int(line[2]),int(line[3]),line[4]
		if 'D' in Type:
			#this is a contig
			contig_count += 1
			if 'U' in prev_Type:
				contig_start = start
		if 'U' in Type:
			#this is a gap
			gap_count += 1
			if 'BLAT' in line[6]:
				blat_gap_count+=1
			#print(start,end)
			if 'D' in prev_Type:
				contig_end = start - 1
				start_stop_list.append([contig_start,contig_end])
		prev_Type = Type

	print('##%s\nCounts\n%i contigs are in the Golden Path for Zoey on %s' % (CHROM,contig_count,CHROM))
	print('%i gaps remain in the Golden Path for Zoey on %s' % (gap_count,CHROM))
	print('Of these %i gaps, %i gap(s) correspond(s) to BLAT gaps (unaligned sequences where gaps were introduced)' %(gap_count,blat_gap_count))

	continuous_contig_array = []
	for i in range(0,len(start_stop_list)):
		start,end = int(start_stop_list[i][0]),int(start_stop_list[i][1])
		length = end-start
		continuous_contig_array.append(length)

	Average_continuity = format(np.mean(continuous_contig_array), '.2f') 
	print('\n##Statistics\nAverage continuity (bp): %s ' % Average_continuity)  
	print('STEP SIX: DONE!')



	####################################################
	#STEP SEVEN - CALCULATE CHROM CONTINUITY FOR CANFAM3
	####################################################

	count, Average_continuity = calculate_canFam_continuity(CHROM)
	print('There are %i gaps in the canFam3 assembly (not including chrUn)' % count)
	print('Average continuity genome-wide from canFam3 (bp): ', Average_continuity)  
	
	#break
	logFile.close()
	print('\n####\n%s DONE!!!\n#####\n\n\n\n' % (CHROM))