import cgi
import webapp2
from genescript4 import Genescript
import os
import urllib
import jinja2
import StringIO
import zipfile
from google.appengine.ext import blobstore
from google.appengine.ext.webapp import blobstore_handlers
from Bio import SeqIO

JINJA_ENVIRONMENT = jinja2.Environment( # Loads Jinja2 templating engine - takes python variables and plugs values into HTML template.
    loader=jinja2.FileSystemLoader(os.path.dirname(__file__)),
    extensions=['jinja2.ext.autoescape'])

genomes={"ecolimg1655":"Escherichia coli K-12 MG1655","ecolibl21":"Escherichia coli BL21"}
type2s={"AarI":'CACCTGC',"DerII":'TTTCCG'}
blunts={"SmaI":'CCCGGG',"DasIII":'AAATTT'}
	
class MainPage(webapp2.RequestHandler):
	def get(self):
		self.redirect('/form/')

class FormHandler(webapp2.RequestHandler):
    def get(self, resource = None):
		upload_url = blobstore.create_upload_url('/upload')
		
		custom_genome_key = genbank_data = None
		if resource != None:
			custom_genome_key = str(urllib.unquote(resource))
			
		template_values = {
			'genomes':genomes,
			'type2s':type2s,
			'blunts':blunts,
			'upload_url': upload_url,
			'custom_genome_key': custom_genome_key,
			'genbank_data': genbank_data
		}
		template = JINJA_ENVIRONMENT.get_template('index.html')
		self.response.write(template.render(template_values))



class UploadHandler(blobstore_handlers.BlobstoreUploadHandler):
	def post(self):
		upload_files = self.get_uploads('custom_genome_name')
		blob_info = upload_files[0]
		self.redirect('/form/%s' % blob_info.key())

		
class Output(webapp2.RequestHandler):
    def post(self):
	
		script = Genescript()
		
		# load genome
		genome = self.request.get("genome")
		if genome[0:7] == "custom:":
			# load from blobstore
			blob_key = genome[7:]
			self.response.out.write('Loading from Blobstore:'+blob_key+'\n')
						
			blob_reader = blobstore.BlobReader(blob_key)
			genbank_data = blob_reader.read()
			seq_record = SeqIO.read(StringIO.StringIO(genbank_data), "genbank")
			script.genome = seq_record

		else:
			script.load_genome(genome)
		
		# set gene name(s)
		gene_names = self.request.get("gene_name")
		script.set_gene_names(gene_names)
		
		# set type2s
		type2sname = self.request.get("type2s_site")
		if type2sname != "_custom_":
			type2sseq = type2s[type2sname]
		else:
			type2sname = self.request.get("type2s_custom_name")
			type2sseq = self.request.get("type2s_custom_seq")
			
		script.set_type2s_seq(type2sname,type2sseq)
		
		# set blunt
		bluntname = self.request.get("blunt_site")
		if bluntname != "_custom_":
			bluntseq = blunts[bluntname]
		else:
			bluntname = self.request.get("blunt_custom_name")
			bluntseq = self.request.get("blunt_custom_seq")
			
		script.set_blunt_seq(bluntname,bluntseq)
		
		# optional params
		script.ust_maxlen=int(self.request.get('max_ust'))
		script.ust_optlen=int(self.request.get('opt_ust'))
		script.ust_minlen=int(self.request.get('min_ust'))

		script.dst_maxlen=int(self.request.get('max_dst'))
		script.dst_optlen=int(self.request.get('opt_dst'))
		script.dst_minlen=int(self.request.get('min_dst'))
	
		# run script
		script.get_all()
		
		# output
		action = self.request.get('formaction')
		first_gene_name = script.input_list[0]	
		
		if action == "save-ust":
			filename = first_gene_name+"_UST.gb"
			self.response.headers['Content-Type'] ='application/octet-stream'
			self.response.headers['Content-Disposition'] ='attachment; filename="'+str(filename)+'"'
			self.response.out.write(script.results[first_gene_name]['ust'])

		elif action == "save-dst":
			filename = first_gene_name+"_DST.gb"
			self.response.headers['Content-Type'] ='application/octet-stream'
			self.response.headers['Content-Disposition'] ='attachment; filename="'+str(filename)+'"'
			self.response.out.write(script.results[first_gene_name]['ust'])
		
		elif action == "save-zip":
		
			self.response.headers["Content-Type"] = "multipart/x-zip"
			self.response.headers['Content-Disposition'] = 'attachment; filename="genes_ust_dst.zip"'
			
			# create new zip file
			output = StringIO.StringIO()
			z = zipfile.ZipFile(output,'w')
			
			# add ust & dst-files per gene
			for gene in script.results.keys():
				z.writestr(gene+"_UST.gb", str(script.results[gene]['ust']))
				z.writestr(gene+"_DST.gb", str(script.results[gene]['dst']))

			# close zip
			z.close()
			
			# send it
			self.response.out.write(output.getvalue())
			output.close()
			
		else:
			self.response.headers['Content-Type'] ='text/plain'
			for gene in script.results.keys():
				self.response.out.write('----------------------- '+gene+' UST ------------------------------\n\n\n\n')
				self.response.out.write(script.results[gene]['ust'])
				self.response.out.write('----------------------- '+gene+' DST ------------------------------\n\n\n\n')
				self.response.out.write(script.results[gene]['dst'])
				self.response.out.write('################################################################\n\n\n\n')
		
class Autocomplete(webapp2.RequestHandler):
	def get(self):
		
		script = Genescript()
		
		genome = self.request.get("genome")
		
		
		# load genome
		genome = self.request.get("genome")
		if genome[0:7] == "custom:":
			# load from blobstore
			blob_key = genome[7:] # strip "custom:"						
			blob_reader = blobstore.BlobReader(blob_key) # init blobstore by key
			genbank_data = blob_reader.read() # read from blobstore
			seq_record = SeqIO.read(StringIO.StringIO(genbank_data), "genbank")
			script.genome = seq_record

		else:
			script.load_genome(genome)
		
		gene_list = script.count_genes()
		#families=script.gather_operons(gene_list)
		#operon_list=script.get_operon_names(families)
		#combo_list=script.get_op_combs(operon_list)

		output = '['
		
		list = []
		#operonout=[]
		
		for gene in gene_list:
			list.append('"'+gene[0] + ' (' + gene[1] + ')"')
					
		# for entry in combo_list:  #includes combinations of operons to autocomplete list. 
			# list.append('"'+ entry +'"')
			# #print entry
		
		output += str.join(',',list)
		
		
		output += ']'
		#output = '[ "aaa", "bbb"]'
		#for gene in gene_list:
		self.response.out.write(output)	
			

class Ajax(webapp2.RequestHandler):
	def get(self):	
		self.response.headers['Content-Type'] ='text/plain'
		for item in self.request.GET.items():
			self.response.out.write("GET: " + item[0] + " = " + item[1] + "\n")
		for item in self.request.POST.items():
			self.response.out.write("POST: " + item[0] + " = " + item[1] + "\n")
		
		#self.response.out.write("HELLO\n\nHow are u\t\txxxxx")
		output = Output(webapp2.RequestHandler)
		output.post()
		
		
app = webapp2.WSGIApplication([('/', MainPage),
								('/form/([^/]+)?', FormHandler),
								('/form', FormHandler),
							  ('/submit', Output),
							  ('/autocomplete', Autocomplete),
							  ('/ajax', Ajax),
							  ('/upload',UploadHandler)
							  ],
                              debug=True)
