__author__ = 'np'

# from django.template import Template
# from django.template import Context
import requests
import json
import logging
import os
#import redis
import datetime
import pytz


class Calculator(object):
	"""
	Skeleton class for calculators
	"""
	def __init__(self, calc=None):
		self.name = ''
		self.propMap = {}
		self.baseUrl = None
		self.urlStruct = ''
		self.results = ''
		self.headers = {'Content-Type': 'application/json'}
		self.request_timeout = 30  # default, set unique ones in calc sub classes
		self.max_retries = 3

		self.image_scale = 50

		self.default_ph = 7.0

		#self.redis_hostname = os.environ.get('REDIS_HOSTNAME')
		#self.redis_port = os.environ.get('REDIS_PORT')
		#self.redis_conn = redis.StrictRedis(host=self.redis_hostname, port=self.redis_port, db=0)

		self.jchem_server_url = os.environ.get('CTS_JCHEM_SERVER', 'localhost:8080')
		self.efs_server_url = os.environ.get('CTS_EFS_SERVER', 'localhost:8080')

		# jchem ws urls:
		self.export_endpoint = '/webservices/rest-v0/util/calculate/molExport'
		self.detail_endpoint = '/webservices/rest-v0/util/detail'
		self.type_endpoint = '/webservices/rest-v0/util/analyze'
		self.checker_endpoint = '/webservices/rest-v0/util/calculate/structureChecker'

		# CTSWS (formerlly EFS) metabolizer endpoints:
		self.efs_metabolizer_endpoint = '/ctsws/rest/metabolizer'
		self.efs_standardizer_endpoint = '/ctsws/rest/standardizer'

		# cts p-chem properties
		self.pchem_props = [
			'boiling_point',
			'melting_point',
			'water_sol',
			'vapor_press',
			'mol_diss',
			'ion_con',
			'henrys_law_con',
			'kow_no_ph',
			'kow_wph',
			'kow_ph',
			'kow'
		]

		# cts chemical information dict
		self.chemical_information = {
			'chemical': None,  # user-entered chemical (as-entered or drawn)
			'orig_smiles': None,  # original conversion to SMILES
			'smiles': None,  # SMILES after filtering, used for calculations
			'formula': None,
			'iupac': None,
			'mass': None,
			'structureData': None,  # drawn chemical structure format for MarvinSketch
			'exactMass': None,
		}

		# cts chemical information request
		self.chemical_information_request = {
			'chemical': None,
			'get_structure_data': False,
		}

		# cts api data object for p-chem data request
		self.data_obj = {
			'calc': None,
			'prop': None,
			'data': None,
			'chemical': None,
		}

		# cts p-chem request object with default key:vals.
		# can handle list of props (ws) or single prop (cts api)
		self.pchem_request = {
			'service': None,
			'chemical': None,
			'prop': None,
			'sessionid': None,
			'method': None,
			'ph': 7.0,
			'node': None,
			'calc': None,
			'run_type': None,
			'workflow': None,
			'mass': None,
			'props': [],
		}

		# cts p-chem response object with defaults, returns one prop per reponse
		self.pchem_response = {
			'chemical': None,
			'calc': None,
			'prop': None,
			'method': None,
			'run_type': None,
			'workflow': None,
			'node': None,
			'request_post': None,
			'data': None,
			'error': False,
		}


	def getUrl(self, prop):
		if prop in self.propMap:
			calcProp = self.propMap[prop]['urlKey']
			return self.urlStruct.format(calcProp)
		else:
			return "Error: url key not found"

	def getPropKey(self, prop):
		if prop in self.propMap:
			return self.propMap[prop]['propKey']
		else:
			return "Error: prop not found"

	def getResultKey(self, prop):
		if prop in self.propMap:
			return self.propMap[prop]['resultKey']
		else:
			return "Error: result key not found"

	def gen_jid(self):
		ts = datetime.datetime.now(pytz.UTC)
		localDatetime = ts.astimezone(pytz.timezone('US/Eastern'))
		jid = localDatetime.strftime('%Y%m%d%H%M%S%f')
		return jid


	# def get_melting_point(self, structure, sessionid, calc=None):
	def get_melting_point(self, structure, sessionid, calc_obj):
		"""
		Gets mass of structure from Measured, tries
		TEST if not available in Measured, and finally EPI.
		Returns MP as float or None
		"""
		melting_point_request = {
			'calc': "",
			'prop': 'melting_point',
			'chemical': structure,
			'sessionid': sessionid
		}

		melting_point = None
		calc = calc_obj.name

		# Attempt at MP workflow as loop..
		mp_request_calcs = ['measured', 'test']  # ordered list of calcs for mp request
		if calc != 'epi':
			# Note: EPI also requests MP, but gets it from itself if it can't from Measured or TEST.
			mp_request_calcs.append('epi')

		if calc == 'test':
			melting_point_request['method'] = "hc"  # method used for MP value

		for calc in mp_request_calcs:

			melting_point_request['calc'] = calc

			logging.info("Requesting melting point from {}..".format(calc))

			# Calls calculator's data_request_handler which makes request to calc server:
			response_obj = calc_obj.data_request_handler(melting_point_request)

			if calc == 'test':
				melting_point = response_obj['data']
			elif not response_obj.get('valid'):
				# epi or measured mp request not valid, sets mp to None
				melting_point = None
			else:
				# Finds mp data from list of data objects for epi or measured:
				for data_obj in response_obj['data']:
					if data_obj['prop'] == "melting_point":
						melting_point = data_obj['data']

			try:
				# results_object = json.loads(mp_response.content)
				melting_point = float(melting_point)
			except Exception as e:
				logging.warning("Unable to get melting point from {}\n Exception: {}".format(calc, e))
				logging.warning("Data returned from Measured that triggered exception: {}".format(response_obj.get('data')))
			if isinstance(melting_point, float):
				logging.info("Melting point value found from {} calc, MP = {}".format(calc, melting_point))
				return melting_point

		# if no MP found from all 3 calcs, returns None for MP
		return None




	################ JCHEM REST STUFF (WHERE SHOULD IT GO??? CTS_REST???) ###################

	def getChemDetails(self, request_obj):
		"""
		getChemDetails

		Inputs:
		chem - chemical name (format: iupac, smiles, or formula)
		Returns:
		The iupac, formula, mass, and smiles string of the chemical
		along with the mrv of the chemical (to display in marvinjs)
		"""
		chemical = request_obj.get('chemical')

		chemDeatsDict = {
			"structures": [
				{"structure": chemical}
			],
			"display": {
				"include": [
					"structureData"
				],
				"additionalFields": {
					"formula": "chemicalTerms(formula)",
					"iupac": "chemicalTerms(name)",
					"mass": "chemicalTerms(mass)",
					"exactMass": "chemicalTerms(exactMass)",
					"smiles": "chemicalTerms(molString('smiles'))",
					# "cas": "chemicalTerms(molString('name:cas#'))",
					"preferredName": "chemicalTerms(molString('name:t'))"
				},
				"parameters": {
					"structureData": "mrv"
				}
			}
		}

		url = self.jchem_server_url + self.detail_endpoint
		return self.web_call(url, chemDeatsDict)


	def smilesToImage(self, request_obj):
		"""
		smilesToImage

		Returns image (.png) url for a 
		given SMILES
		"""
		smiles = request_obj.get('smiles')
		imgScale = request_obj.get('scale', 100)
		imgWidth = request_obj.get('width')
		imgHeight = request_obj.get('height')
		imgType = request_obj.get('type')

		# NOTE: Requesting image without width or height but scale
		# returns an image that just fits the molecule with nice resoltuion.
		# Providing width and height without scale might return a higher
		# resolution image for metabolites!
		
		request = {
			"structures": [
				{"structure": smiles}
			],
			"display": {
				"include": ["image"],
				"parameters": {
					"image": {
					}
				}
			}
		}

		if imgType:
			request['display']['parameters']['image'].update({'type': imgType})
		else:
			request['display']['parameters']['image'].update({'type': 'png'})

		if imgWidth and imgHeight:
			# these are metabolites in the space tree:
			request['display']['parameters']['image'].update({"width": imgWidth, "height": imgHeight})
		elif imgWidth and not imgHeight:
			request['display']['parameters']['image'].update({'width': imgWidth, 'scale': imgScale})
			# request['display']['parameters']['image'].update({'scale': imgScale})
		else:
			request['display']['parameters']['image'].update({'scale': imgScale})

		url = self.jchem_server_url + self.detail_endpoint
		imgData = self.web_call(url, request)  # get response from jchem ws
		return imgData  # return dict of image data


	def convertToSMILES(self, request_obj):
		"""
		convertToSMILES

		Inputs: chemical as mrv, smiles, etc. (chemaxon recognized)
		Returns: SMILES string of chemical
		"""
		chemStruct = request_obj.get('chemical')  # chemical in <cml> format (marvin sketch)
		data = {
			"structure": chemStruct,
			"parameters": "smiles"
		}
		url = self.jchem_server_url + self.export_endpoint
		return self.web_call(url, data)  # get responset))


	def getStructInfo(self, structure):
		"""
		Appends structure info to image url
		Input: structure in .mrv format
		Output: dict with structure's info (i.e., formula, iupac, mass, smiles),
		or dict with aforementioned keys but None values
		"""
		structDict = self.getChemDetails({"chemical": structure, "addH": True})
		infoDictKeys = ['formula', 'iupac', 'mass', 'smiles','exactMass']
		infoDict = {key: None for key in infoDictKeys}  # init dict with infoDictKeys and None vals
		struct_root = {}  # root of data in structInfo
		if 'data' in structDict:
			struct_root = structDict['data'][0]
			infoDict.update({
				"formula": struct_root['formula'],
				"iupac": struct_root['iupac'],
				"mass": struct_root['mass'],
				"smiles": struct_root['smiles'],
				'exactMass': struct_root['exactMass']
			})
		return infoDict


	def getMass(self, request_obj):
		"""
		get mass of structure from jchem ws
		"""
		chemical = request_obj.get('chemical')
		logging.info("jchem_rest getting mass for {}".format(chemical))
		post_data = {
			"structures": [
				{"structure": chemical}
			],
			"display": {
				"include": [
					"structureData"
				],
				"additionalFields": {
					"mass": "chemicalTerms(mass)"
				},
				"parameters": {
					"structureData": "smiles"
				}
			}
		}
		url = self.jchem_server_url + self.detail_endpoint
		return self.web_call(url, post_data)


	def get_chemical_type(self, chemical):
		"""
		Returns type of chemical (e.g., smiles, name, cas, etc.)
		"""
		url = self.jchem_server_url + self.type_endpoint
		request_header = {'Content-Type': "*/*"}
		response, results = None, None
		try:
			response = requests.post(url, data=chemical.encode('utf-8'), headers=request_header, timeout=self.request_timeout)
			results = json.loads(response.content)
		except Exception as e:
			logging.warning("Exception at get_chemical_type: {}".format(e))
			return {'type': None}
		_type_response = {
			'type': None
		}
		check_result = self.check_response_for_errors(results)
		if not check_result.get('valid'):
			# errors found in response..
			_type_response['error'] = check_result['error']
			return _type_response
		if 'properties' in results:
			_type_response['type'] = results['properties'].get('type')
		if not _type_response['type'] and results.get('type'):
			_type_response['type'] = results['type']
		return _type_response

	def get_smiles_from_name(self, chemical):
		"""
		Calls jchem ws /molExport endpoint and assumes the chemical
		is a name being converted into a SMILES. This is to help determine
		the intended chemical format in the event that the chemical is a valid
		SMILES and/or name (e.g., PFOS).
		"""
		_url = self.jchem_server_url + self.export_endpoint
		_post = {
			'structure': chemical,
			'inputFormat': "name",
			'parameters': "smiles"
		}
		_response = {
			'smiles': None,
			'format': "smiles"
		}

		_results = self.web_call(_url, _post)
		_check_results = self.check_response_for_errors(_results)

		if not _check_results.get('valid'):
			_response['error'] = "Not a valid name"
			return _response

		_response['smiles'] = _results.get('structure')
		return _response



	def check_response_for_errors(self, results):
		"""
		Checks for errors in HTTP responses from web_call.
		Typically errors to check are from jchem webservices.
		Inputs:
			results - response content, a JSON object
		Output:
			_checked_response - results object with wrapper for
			handling any possible errors.
		"""
		_errors_list = ['errorMessage', 'errorCode', 'error']
		_check_response = {
			'valid': False,
			'error': None
		}
		_result_keys = []

		if not isinstance(results, dict):
			_check_response['error'] = "Error processing calc results"
			logging.warning("Excepted response object not a dict object.")
			return _check_response

		_result_keys = list(results.keys())  # get keys of result object..

		for key in _result_keys:
			if key in _errors_list:
				_check_response['error'] = self.handle_error_messages(results)

		if not _check_response.get('error'):
			_check_response['valid'] = True  # didn't find any errors..

		return _check_response



	def handle_error_messages(self, results):
		"""
		Handles known error that occurred requesting data from jchem
		"""
		if results.get('errorCode') == 3:
			# jchem ws can't read molecule file..
			return "Chemical cannot be standardized"
		else:
			return "Chemical not recognized"


	def web_call(self, url, data, headers=None):
		"""
		Makes the request to a specified URL
		and POST data. Returns resonse data as dict
		"""
		# TODO: Deal with errors more granularly... 403, 500, etc.

		if not headers:
			headers = self.headers
		try:
			if data == None:
				response = requests.get(url, timeout=self.request_timeout)
			else:
				response = requests.post(url, data=json.dumps(data), headers=headers, timeout=self.request_timeout)

			results = json.loads(response.content)

			valid_object = self.check_response_for_errors(results)

			if valid_object.get('valid'):
				results['valid'] = True
				return results

			else:
				error_response = {
					'error': valid_object.get('error'),
					'data': results,
					'valid': False
				}
				return error_response

			# return json.loads(response.content)
		except requests.exceptions.RequestException as e:
			logging.warning("error at web call: {} /error".format(e))
			raise e





	################# TRANSFORMATION PRODUCTS STUFF (DATA_WALKS) #######################

	def nodeWrapper(self, smiles, height, width, scale, key=None, img_type=None, isProduct=None):
		"""
		Wraps image html tag around
		the molecule's image source
		Inputs: smiles, height, width, scale, key
		Returns: html of wrapped image
		"""

		# 1. Get image from smiles
		post = {
			"smiles": smiles,
			"scale": scale,
			"height": height,
			"width": width,
			# "type": img_type
		}

		if img_type:
			post.update({'type': img_type})

		results = self.smilesToImage(post)

		# 2. Get imageUrl out of results
		img, imgScale = '', ''
		if 'data' in results:
			root = results['data'][0]['image']
			if 'image' in root:
				img = root['image']

			if not height:
				height = root.get('height')  # sets height if not provided

		# 3. Wrap imageUrl with <img>
		# <img> wrapper for image byte string:
		if img_type and img_type == 'svg':

			if key:
				html = '<div class="cts-chem-wrap" style="background-color:white;"' + 'id="' + str(key) + '">' + img + '</div>'            
			else:
				html = '<div class="cts-chem-wrap" style="background-color:white;">' + img + '</div>'
			
		else:

			context = {'smiles': smiles, 'img': img, 'height': height, 'width': width, 'scale': scale, 'key': key}
			# html = self.imgTmpl(isProduct).render(Context(context))
			html = self.imgTmpl2(context, isProduct)

		return html


	def imgTmpl2(self, data, isProduct):
		"""
		Creates <img> without Django templates.
		"""
		if isProduct:
			imgTmpl = """
			<img class="metabolite" id="{}"
				alt="{}" src="data:image/png;base64,{}"
				width={} height={} /> 
			""".format(data['key'], data['smiles'], data['img'], data['width'], data['height'])
		else:
			imgTmpl = """
			<img class="metabolite hidden-chem" id="{}"
				alt="{}" src="data:image/png;base64,{}"
				width={} height={} hidden /> 
			""".format(data['key'], data['smiles'], data['img'], data['width'], data['height'])
		imgTmpl = imgTmpl.replace('\t', '').replace('\n', '')
		return imgTmpl


	# def imgTmpl(self, isProduct):
	# 	"""
	# 	Uses Django templates to create <img> for molecule.
	# 	"""
	# 	if isProduct:
	# 		imgTmpl = """
	# 		<img class="metabolite" id="{{key|default:""}}"
	# 			alt="{{smiles}}" src="data:image/png;base64,{{img}}"
	# 			width={{width}} height={{height}} /> 
	# 		"""        
	# 	else:
	# 		imgTmpl = """
	# 		<img class="metabolite hidden-chem" id="{{key|default:""}}"
	# 			alt="{{smiles}}" src="data:image/png;base64,{{img}}"
	# 			width={{width}} height={{height}} hidden /> 
	# 		"""
	# 	imgTmpl = imgTmpl.replace('\t', '').replace('\n', '')
	# 	return Template(imgTmpl)


	def popupBuilder(self, root, paramKeys, molKey=None, header=None, isProduct=False):
		"""
		Wraps molecule data (e.g., formula, iupac, mass, 
		smiles, image) for hover-over popups in chemspec
		and gentrans outputs.

		Inputs:
		root - dictionary of items to wrap in table
		paramKeys - keys to use for building table
		molKey - (optional) add id to wrap table
		header - (optional) add header above key/values 

		Returns: dictionary where html key is 
		the wrapped html and the other keys are
		same as the input keys
		"""
		dataProps = {key: None for key in paramKeys}  # metabolite properties
		html = '<div id="{}_div" class="nodeWrapDiv"><div class="metabolite_img" style="float:left;">'.format(molKey)

		# smiles, height, width, scale, key=None, img_type=None
		if isProduct:
			html += self.nodeWrapper(root['smiles'], None, 250, self.image_scale, molKey, 'svg')  # svg popups for chemspec and gentrans outputs
			html += self.nodeWrapper(root['smiles'], None, None, self.image_scale, molKey, None)  # hidden png for pdf
		else:
			html += self.nodeWrapper(root['smiles'], None, None, self.image_scale, molKey, 'png', True)  # NOTE: testing just png for popups to fix missing lines in svgs
			html += self.nodeWrapper(root['smiles'], None, None, self.image_scale, molKey, None)  # hidden png for pdf

		html += '</div>'

		if molKey:
			html += '<table class="ctsTableStylin" id="{}_table">'.format(molKey)
		else:
			html += '<table class="ctsTableStylin">'

		if header:
			html += '<tr class="header"><th colspan="2">' + header + '</th></tr>'

		for key, value in root.items():
			if key in paramKeys:
				# Convert other types (e.g., float, int) to string
				if isinstance(value, float):
					if key == 'exactMass':
						value = str(value)
					else:
						value = str(round(float(value), 3))
				dataProps[key] = value
				html += '<tr><td>' + key + '</td>'
				html += '<td>' + value + '</td></tr>'
		html += '</table></div>'

		dataProps["html"] = html

		return dataProps