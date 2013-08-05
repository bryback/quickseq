import os
import sys
import csv
import urllib2
from StringIO import StringIO
from string import lower
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import Entrez


class Genescript:
	genome=None
	organism = 'Escherichia coli strain mg1655'
	upstream=1000
	downstream=1000
	ust_homol=500
	dst_homol=500
	min_ust=50
	min_dst=50
	start_len=21
	stop_len=21
	input_list=[]
	type2s_name=None
	blunt_name=None
	spacer='ACAC'
	cut_length=4
	ust_maxlen=1000
	ust_optlen=500
	ust_minlen=50
	ust_blunt_site=None
	ust_type2s_site=None
	dst_maxlen=1000
	dst_optlen=500
	dst_minlen=50
	dst_blunt_site=None
	dst_type2s_site=None
	results_dict={}
	results = {}
	
	def __init__(self):
		self.set_type2s_seq('AarI','CACCTGC')
		self.set_blunt_seq('SmaI','GGGCCC')
		self.results_dict['ust_error']='asdffdsa'
	def set_gene_names(self, user_input):
		self.input_list=[]
		for n in user_input.split(','):
			item = n.lstrip(' ')
			item = item.split(' ')
			item = item[0]
			self.input_list.append(item)
				
	def set_type2s_seq(self,type2s_name,type2s_seq):
		self.ust_type2s_site=Seq(type2s_seq,IUPAC.unambiguous_dna)
		self.dst_type2s_site=Seq(type2s_seq,IUPAC.unambiguous_dna)
		self.type2s_name=type2s_name
	
	def set_blunt_seq(self,blunt_name,blunt_seq):
		self.ust_blunt_site=Seq(blunt_seq,IUPAC.unambiguous_dna)
		self.dst_blunt_site=Seq(blunt_seq,IUPAC.unambiguous_dna)
		self.blunt_name=blunt_name
		
	def get_all(self):
		gene_list = self.input_list
		for gene_name in gene_list:
			gene_input=self.operon_handler(gene_name)
			print 'Generating sequences for',gene_name,'\n'
			result = {
				'ust': self.ust_output(gene_input[0]),
				'dst': self.dst_output(gene_input[1])
			}
			self.results[gene_name] = result
		print '\nProgram complete.'
		return self.results
		
	def ust_output(self,gene_name):
		ust=self.get_ust(gene_name)
		start=self.get_start_seq(gene_name)
		ust_length=len(ust)
		uststart=ust+start
		gene_id= self.get_gene_id(gene_name)
		ust_features=self.build_ust_features(gene_name,ust_length)
		ust_rec=SeqRecord(seq=uststart, id=gene_id, name=gene_name+'_UST_'+self.type2s_name, features=ust_features)  
		#SeqIO.write(ust_rec, "UST_"+gene_name+".gb", "genbank")
		out_handle=StringIO()
		SeqIO.write(ust_rec, out_handle, "genbank")
		#dataout=out_handle.getvalue()
		#self.results_dict['ust_data']=dataout
		return out_handle.getvalue()
		
	def dst_output(self,gene_name):
		dst=self.get_dst(gene_name)
		stop=self.get_stop_seq(gene_name)
		dst_length=len(dst)
		dststop=stop+dst
		gene_id=self.get_gene_id(gene_name)
		dst_features=self.build_dst_features(gene_name, dst_length)
		dst_rec=SeqRecord(seq=dststop, id=gene_id, name=gene_name+'_DST_'+self.type2s_name, features=dst_features)  
		out_handle=StringIO()
		SeqIO.write(dst_rec, out_handle, "genbank")
		#dataout=out_handle.getvalue()
		#self.results_dict['dst_data']=dataout
		return out_handle.getvalue()
		
	def get_ust(self,gene_name):
	# fetches sequence "maxlen" (maximum length) bases upstream of start codon and trims sequence until all illegal restriction sites are excluded.
	# Then locates the 'GGG' motif closest to "optlen" (optimal length) and returns the sequence between that GGG and the start codon, if it is not smaller than minlen.
		maxlen=self.ust_maxlen
		optlen=self.ust_optlen
		minlen=self.ust_minlen
		blunt_site=self.ust_blunt_site
		type2s_site=self.ust_type2s_site
		print 'Processing upstream homology region for', gene_name,'...\n'
		raw_ust=self.concat_seq(gene_name,maxlen)['raw_ust']
		target_blunt_site=blunt_site[:len(blunt_site)/2]
		x=0
		rsite=0
		lsite=0
		while raw_ust[x:].count(blunt_site)>0 or raw_ust[x:].count(blunt_site.reverse_complement())>0 or raw_ust[x:].count(type2s_site)>0 or raw_ust[x:].count(type2s_site.reverse_complement())>0:
			x+=1
		else:
			nosite_ust=raw_ust[x:]
		#print 'pruned UST length:',len(nosite_ust)
		if str(nosite_ust[:-optlen]).count(str(target_blunt_site))>0:
			rsite=str(nosite_ust[:-optlen]).rindex(str(target_blunt_site))
			#print 'Rsite: ',rsite
		if str(nosite_ust[-optlen:]).count(str(target_blunt_site))>0 and len(nosite_ust)-optlen>0:
			lsite=str(nosite_ust[-optlen:]).index(str(target_blunt_site))+(len(nosite_ust)-optlen)
		elif str(nosite_ust[-optlen:]).count(str(target_blunt_site))>0 and len(nosite_ust)-optlen<=0:
			lsite=str(nosite_ust[-optlen:]).index(str(target_blunt_site))
			#print 'Lsite: ',lsite
			
		if len(nosite_ust)>=optlen:
			if abs(optlen-rsite)<=abs(optlen-lsite):
				ust=nosite_ust[rsite:]
			else:
				ust=nosite_ust[lsite:]
		else:
			ust=nosite_ust[lsite:]

		if len(ust)<minlen:
			error='Error! Upstream sequence length below minimum. Please check sequence and parameters.\n' \
			+'Minimum length:' + str(minlen) + 'UST sequence length:' + str(len(ust))+'.'
			self.results_dict['ust_error']=error
			print error
		else:
			correct_output= 'UST length:'+str(len(ust))+'\n'+ust
			#self.results_dict['ust_data']=correct_output
			print correct_output
			return ust

	def get_dst(self,gene_name):
		# Same as get_ust, with adjusted indexing.
		maxlen=self.dst_maxlen
		optlen=self.dst_optlen
		minlen=self.dst_minlen
		blunt_site=self.dst_blunt_site
		type2s_site=self.dst_type2s_site
		print 'Processing downstream homology region for', gene_name,'...\n'
		raw_dst=self.concat_seq(gene_name,maxlen)['raw_dst']
		target_blunt_site=blunt_site[len(blunt_site)/2:]
		x=-1
		rsite=0
		lsite=0
		while raw_dst[:x].count(blunt_site)>0 or raw_dst[:x].count(blunt_site.reverse_complement())>0 or raw_dst[:x].count(type2s_site)>0 or raw_dst[:x].count(type2s_site.reverse_complement())>0:
			x-=1
		else:
			if x+1==0:
				nosite_dst=raw_dst
			else:
				nosite_dst=raw_dst[:x]
		#print 'pruned DST length:', len(nosite_dst)
		if str(nosite_dst[optlen:]).count(str(target_blunt_site))>0:
			rsite=str(nosite_dst[optlen:]).index(str(target_blunt_site))+optlen+len(target_blunt_site)
			#print 'Rsite: ',rsite
		if str(nosite_dst[:optlen]).count(str(target_blunt_site))>0:
			lsite=str(nosite_dst[:optlen]).rindex(str(target_blunt_site))+len(target_blunt_site)
			#print 'Lsite: ',lsite
			
		if len(nosite_dst)>=optlen:
			if abs(optlen-rsite)<=abs(optlen-lsite):
				dst=nosite_dst[:rsite]
			else:
				dst=nosite_dst[:lsite]
		else:
			dst=nosite_dst[:lsite]

		if len(dst)<minlen:
			print 'Error! Downstream sequence length below minimum. Please check sequence and parameters.'
			print 'Minimum length:', minlen
			print 'DST sequence length:', len(dst)
		else:
			print 'DST length:',len(dst)
			print dst
			return dst
	
	def count_genes(self): # counts the genome features with type='gene' and adds them to the gene_list as a list with [gene_name,synonyms,db_xref]
		genome=self.genome
		gene_list=[]
		for entry in range(len(genome.features)):
			if genome.features[entry].type=='gene':
				gene_name=self.genome.features[entry].qualifiers['gene'][0]
				altnames=self.genome.features[entry].qualifiers['gene_synonym'][0]
				dbxref=self.genome.features[entry].qualifiers['db_xref'][0]
				gene_list.append([gene_name,altnames,dbxref])
		
		return gene_list
	
	def gather_operons(self,gene_list): # takes list of all genes and counts returns a dictionary with
	# families, the keys are the first 3 letters of the gene names in that family, each value is a list
	# containing all the gene names as 4-tuples with their location data.
		import collections
		genome=self.genome
		short_list=[]
		families={}
		
		for gene in gene_list:
			short_list.append(gene[0][:3])
		
		duplicates = [x for x, y in collections.Counter(short_list).items() if y > 1]
		
		for entry in duplicates:
			families[entry]=[]
			for gene in gene_list:
				if entry == gene[0][:3]:
					gene_loc=self.locate_gene(gene[0])
					gene_info=(gene[0],gene_loc[0],gene_loc[1],gene_loc[2])
					families[entry].append(gene_info)
		
		return families
		
	def sort_operon(self,homonym_list):	# takes list of genes from family dict, sorts and returns a list of 
# pseudo-operons (lists of 4-tuples of genes and location data).	
		homonym_list.sort(key=lambda tup: tup[1], reverse=False) # sorts genes according to their start location on the genome
		gene_list=homonym_list
		operon_name=[]
		operon_list=[]
		
		for index in range(len(gene_list))[1:]:
			if gene_list[index][3]==gene_list[index-1][3] and abs(gene_list[index][1]-gene_list[index-1][2])<1000:						
				if operon_name.count(gene_list[index-1])==0:
					operon_name.append(gene_list[index-1])
				if operon_name.count(gene_list[index])==0:
					operon_name.append(gene_list[index])
			else:
				if operon_name!=[]:
					operon_list.append(operon_name)
					operon_name=[]
		if operon_list==[] and operon_name!=[]:
			operon_list.append(operon_name)
		#for operon in operon_list:
			
		return operon_list	
	
	def get_operon_names(self,families): # takes families dict {'ara':[('araD', 65854, 66550, -1), ('araA', 66834, 68337, -1)]...} returns
	# list of operon_strings ['allR', 'alaW', 'ptsI', 'scpBC', 'mngB', 'gadB', 'mcrB', 'fabDGF']
		famgenes=[]
		opstrings=[]
		for key in families.keys():
			famgenes=families[key]
			#print famgenes
			oplist=self.sort_operon(famgenes)
			#print oplist
			for operon in oplist:
				if operon!=[]:
					opname=self.name_operon(operon)
					opstrings.append(opname)
		return opstrings

	def name_operon(self,op_list): #takes a list of genes for a specific operon
	# and returns the complete operon name (e.g. 'LacIZYA')

		opname=op_list[0][0][:3]
				
		for i in range(len(op_list)):
			
			if op_list[i][3]==1:
				opname+=op_list[i][0][3:]
			elif op_list[i][3]==-1:
				opname+=op_list[-i-1][0][3:]
		
		return opname	
		
	def get_op_combs(self,opstrings): # takeslist of full operons, returns list with subslices autocomplete.
		op_combs=[]
		for operonstring in opstrings:
			op_comb=self.operon_combinations(operonstring)
			for opstr in op_comb:
				op_combs.append(opstr)
		return op_combs
		
	
	def operon_combinations(self,opstring): #input='yjiMLKJH', output= ['yjiMLKJH', 'yjiLKJH', 'yjiKJH', 'yjiJH']
		op_combs=[]
		if len(opstring)>4:
			for index in range(len(opstring)-3)[:-1]:
				op_substr=opstring[:3]+opstring[3+index:]
				op_combs.append(op_substr)
		
		return op_combs


		
		
	def get_gene_id(self,gene_name): # Analogous to concat_seq - calls fetch_id fo genes and operons and returns the GeneID. 
		nameholder=self.operon_handler(gene_name)
		if nameholder[0]==nameholder[1]:
			return self.fetch_id(gene_name)      
		else:
			gene1=self.fetch_id(nameholder[0])
			gene2=self.fetch_id(nameholder[1])
			operon_concat=gene1+'-'+gene2
			return operon_concat
	
	def operon_handler(self,gene_name): # Checks syntax of a target input (string). If input is an operon (e.g. 'araBAD') returns a list (e.g. ['araB','araD'])
	# if input is a single gene name with an uppercase first letter (e.g. 'AraB') returns string with lowercase first letter ('araB')
		if gene_name[:-1].islower():
			return [gene_name, gene_name]
		else:
			operon=[]
			x=0
			for n in gene_name:
				x+=1
				if n.isupper() and x>1:
					operon.append(gene_name[:x-1].lower()+gene_name[x-1:x])
					operon.append(gene_name[:x-1].lower()+gene_name[-1:])
					return operon        
	
	def fetch_id(self,gene_name): # Uses gene locus to extract GeneID of a gene.
		try:
			geneloc=self.locate_gene(gene_name)
			gene=self.genome[geneloc[0]:geneloc[1]]
			gene_id=gene.features[0].qualifiers['db_xref'][-1]
		except TypeError:
			print 'No valid gene id for was found. Please verify that the gene name is correct.'    
		return gene_id

	def locate_gene(self,gene): # searches local genome annotations for input gene name, returns a list with [start,end,orientation
		genome=self.genome
		geneloc=0
		for entry in range(len(genome.features)):
			if genome.features[entry].type=='gene' and genome.features[entry].qualifiers['gene'][0]==gene:
				geneloc = [genome.features[entry].location.start.position,genome.features[entry].location.end.position,genome.features[entry].location.strand]               
		if geneloc!=0:
			return geneloc
		else:
				print gene + ' not found! Please verify that a valid gene name was entered.'

	def build_ust_features(self,gene_name, ust_length):
		nameholder=self.operon_handler(gene_name)
		print "\nGenerating upstream feature annotations for %s\n" % nameholder[0]
		qual_dict=self.build_ust_qualifiers(nameholder[0], ust_length)
		ust_feature_list=['UST_name','gene_name','cutsite','spacer','bs']
		ust_features=[]
		for n in range(len(ust_feature_list)):
			ust_features.append(SeqFeature(FeatureLocation(qual_dict[ust_feature_list[n]][0],qual_dict[ust_feature_list[n]][1]),type=qual_dict[ust_feature_list[n]][3], strand=qual_dict[ust_feature_list[n]][4]))
		for n in range(len(ust_feature_list)):
			ust_features[n].qualifiers['label']=qual_dict[ust_feature_list[n]][2].strip('"')                                                                           
		return ust_features

	def build_ust_qualifiers(self,gene_name, seq_len):
		type2s_site=self.ust_type2s_site
		type2s_name=self.type2s_name
		spacer=self.spacer
		cut_length=self.cut_length
		start_len=self.start_len
		qual_dict={'UST_name' : [0, seq_len,'UST_'+gene_name, "5'UTR", 1],\
		   'gene_name' : [seq_len, seq_len+start_len, gene_name, "gene", 1],\
		   'cutsite' : [seq_len+start_len-cut_length/2, seq_len+start_len+cut_length/2,type2s_name+'_cutsite', "misc_feature", -1],\
		   'spacer' : [seq_len+start_len+cut_length/2, seq_len+start_len+cut_length/2+len(spacer),type2s_name+'_spacer', "misc_feature", -1], \
		   'bs' : [seq_len+start_len+cut_length/2+len(spacer),seq_len+start_len+cut_length/2+len(spacer)+len(type2s_site), type2s_name+'_binding_site', "protein_bind", -1]}
		return qual_dict

	def build_dst_features(self,gene_name, dst_length):
		nameholder=self.operon_handler(gene_name)
		print "\nGenerating downstream feature annotations for %s\n" % nameholder[1]
		qual_dict=self.build_dst_qualifiers(nameholder[1], dst_length)
		dst_feature_list=['DST_name','gene_name','cutsite','spacer','bs']
		dst_features=[]
		for n in range(len(dst_feature_list)):
			dst_features.append(SeqFeature(FeatureLocation(qual_dict[dst_feature_list[n]][0],qual_dict[dst_feature_list[n]][1]),type=qual_dict[dst_feature_list[n]][3], strand=qual_dict[dst_feature_list[n]][4]))
		for n in range(len(dst_feature_list)):
			dst_features[n].qualifiers['label']=qual_dict[dst_feature_list[n]][2].strip('"')                                                                           
		return dst_features

	def build_dst_qualifiers(self,gene_name, seq_len):
		type2s_site=self.dst_type2s_site
		type2s_name=self.type2s_name
		spacer=self.spacer
		cut_length=self.cut_length
		stop_len=self.stop_len
		qual_dict={'DST_name' : [len(type2s_site)+len(spacer)+cut_length/2+stop_len, len(type2s_site)+len(spacer)+cut_length/2+stop_len+seq_len,'DST_'+gene_name, "3'UTR", 1],\
		   'gene_name' : [len(type2s_site)+len(spacer)+cut_length/2, len(type2s_site)+len(spacer)+cut_length/2+stop_len, gene_name, "gene", 1],\
		   'cutsite' : [len(type2s_site)+len(spacer), len(type2s_site)+len(spacer)+cut_length,type2s_name+'_cutsite', "misc_feature", 1],\
		   'spacer' : [len(type2s_site),len(type2s_site)+len(spacer),type2s_name+'_spacer', "misc_feature", 1], \
		   'bs' : [0,len(type2s_site), type2s_name+'_binding_site', "protein_bind", 1]}
		return qual_dict
	

	def concat_seq(self,gene_name,upstream=1000,downstream=1000,start_len=21,stop_len=21): #takes gene name and integers and returns sequence with raw upstream (default = 1000 bp) + first few translated bases (default = 21 bp)
		# as well as the raw downstream region and stop sequence and returns them in a dictionary. For later, move defaults to class global.
		nameholder=self.operon_handler(gene_name)
		if nameholder[0]==nameholder[1]:
			whole_seq=self.fetch_sequence(gene_name,upstream,downstream)
			gene_concat={'raw_ust' : whole_seq[:upstream].seq,
						 'start' : whole_seq[upstream:upstream+start_len].seq,
						 'stop' : whole_seq[-(stop_len+downstream):-downstream].seq,
						 'raw_dst' : whole_seq[-downstream:].seq}
		else:
			gene1=self.fetch_sequence(nameholder[0],upstream,0)
			gene2=self.fetch_sequence(nameholder[1],0,downstream)
			gene_concat={'raw_ust' : gene1[:upstream].seq,
						 'start' : gene1[upstream:upstream+start_len].seq,
						 'stop' : gene2[-(stop_len+downstream):-downstream].seq,
						 'raw_dst' : gene2[-downstream:].seq}
		return gene_concat

	def fetch_sequence(self,gene_name,upstream,downstream): # fetches sequence of query gene with arguments defining length of flanking sequences. Returns list with sequence (coding strand).
		geneloc=self.locate_gene(gene_name)
		if geneloc[2]==1:
			gene_seq=self.genome[geneloc[0]-upstream:geneloc[1]+downstream]
		elif geneloc[2]==-1:
			gene_seq=self.genome[geneloc[0]-downstream:geneloc[1]+upstream].reverse_complement()        
		return gene_seq

	def get_start_seq(self,gene_name):
		# extracts coding sequence from raw sequence and adds typ2s restriction site as well as 2 nucleotides to merge the start seq with the stop seq in frame.
		type2s_site=self.ust_type2s_site
		spacer=self.spacer
		cut_length=self.cut_length
		nameholder = self.operon_handler(gene_name)
		raw_start = self.concat_seq(nameholder[0])['start']
		raw_stop=self.concat_seq(nameholder[1])['stop']
		start_seq=raw_start+raw_stop[:cut_length/2]+spacer+type2s_site.reverse_complement()
		stop_seq=type2s_site+spacer+raw_start[-cut_length/2:]+raw_stop
		return start_seq

	def get_stop_seq(self,gene_name):
		# Analagous to get_start_seq
		type2s_site=self.dst_type2s_site
		spacer=self.spacer
		cut_length=self.cut_length
		nameholder = self.operon_handler(gene_name)
		raw_start=self.concat_seq(nameholder[0])['start']
		raw_stop=self.concat_seq(nameholder[1])['stop']
		stop_seq=type2s_site+spacer+raw_start[-cut_length/2:]+raw_stop
		return stop_seq

	def load_genome(self,localgenome): #opens a genbank file containing a (bacterial) genome from the working directory.
                #This allows much faster queries and sequence retrieval than querying Entrez, but requires the target genome to be locally available.
		print '\nLoading', localgenome
		genomefile = "genomes/"+localgenome+".gb"
		for seq_record in SeqIO.parse(genomefile, "genbank"):
			print "Genome ID:",seq_record.id
			print "Genome size:",len(seq_record),"bp"
		self.genome = seq_record	

