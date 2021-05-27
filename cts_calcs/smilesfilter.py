import requests
import json
import logging
import os
from .calculator import Calculator
from .jchem_properties import Tautomerization, ElementalAnalysis



class SMILESFilter(object):
	"""
	This is the smilesfilter.py module as a class and
	clumped in with other classes related to chem info.
	"""

	def __init__(self):
		self.max_weight = 1500  # max weight [g/mol] for epi, test, and sparc
		self.excludestring = [".","[Ag]","[Al]","[As","[As+","[Au]","[B]","[B-]","[Br-]","[Ca]",
						"[Ca+","[Cl-]","[Co]","[Co+","[Fe]","[Fe+","[Hg]","[K]","[K+","[Li]",
						"[Li+","[Mg]","[Mg+","[Na]","[Na+","[Pb]","[Pb2+]","[Pb+","[Pt]",
						"[Sc]","[Si]","[Si+","[SiH]","[Sn]","[W]"]
		self.return_val = {
			"valid" : False,
			"smiles": "",
			"processedsmiles" : ""
		}
		self.baseUrl = os.environ.get('CTS_EFS_SERVER')
		self.is_valid_url = self.baseUrl + '/ctsws/rest/isvalidchemical'



	def is_valid_smiles(self, smiles):
		"""
		Makes request to ctsws /isvalidchemical endpoint to check
		if user smiles is valid. Returns boolean.
		"""
		is_valid_response = requests.post(self.is_valid_url, data=json.dumps({'smiles': smiles}), headers={'Content-Type': 'application/json'}, timeout=5)
		is_valid = json.loads(is_valid_response.content).get('result')  # result should be "true" or "false"
		if is_valid == "true":
			return True
		else:
			return False



	def check_for_carbon(self, smiles):
		"""
		Makes request to jchem_properties's ElementalAnalysis class,
		which returns the composition of a chemical from JchemWS
		elemental analysis endpoint.
		"""

		# Makes request to get chemical composition:
		analysis_class = ElementalAnalysis()
		analysis_class.make_data_request(smiles, analysis_class)  # sets 'results' attr to json object of response
		chemical_composition = analysis_class.get_elemental_analysis()  # returns list of chemical components

		# Looks through composition data until carbon is found:
		for composite_data in chemical_composition:
			# Gets element out of composite result (e.g., "C (44.34%)"):
			element = composite_data.split("(")[0].replace(" ", "")
			if element == "C":
				return True

		return False



	def check_smiles_against_exludestring(self, smiles):
		"""
		Checks if user's smiles contains any characters from
		the excludestring list. Returns True if smiles is valid,
		and False if the smiles has one of the exluded characters.
		"""
		# for exclude_char in self.excludestring:
		# NOTE: Now just checking for ionic bond (i.e., ".") since CTSWS checks metals:
		if self.excludestring[0] in smiles:
			return False
		return True



	def singleFilter(self, request_obj):
		"""
		Calls single EFS Standardizer filter
		for filtering SMILES
		"""
		smiles = request_obj.get('smiles')
		action = request_obj.get('action')
		post_data = {
			"structure": smiles,
			"actions": [
				action
			]
		}
		calc = Calculator()
		url = calc.efs_server_url + calc.efs_standardizer_endpoint
		return calc.web_call(url, post_data)



	def filterSMILES(self, smiles, is_node=False):
		"""
		cts ws call to jchem to perform various
		smiles processing before being sent to
		p-chem calculators
		"""
		calc_object = Calculator()

		# Performs carbon check (but not for transformation products):
		if not is_node and not self.check_for_carbon(smiles):
			return {'error': "CTS only accepts organic chemicals"}

		# Checks SMILES for invalid characters:
		if not self.check_smiles_against_exludestring(smiles):
			return {'error': "Chemical cannot be a salt or mixture"}

		# Calls CTSWS /isvalidchemical endpoint:
		if not self.is_valid_smiles(smiles):
			logging.warning("User chemical contains metals, sending error to client..")
			return {'error': "Chemical cannot contain metals"}

		# Updated approach (todo: more efficient to have CTSWS use major taut instead of canonical)
		# 1. CTSWS actions "removeExplicitH" and "transform".
		url = calc_object.efs_server_url + calc_object.efs_standardizer_endpoint
		post_data = {
			'structure': smiles,
			'actions': [
				"removeExplicitH",
				"transform"
			]
		}
		response = calc_object.web_call(url, post_data)

		filtered_smiles = response['results'][-1] # picks last item, format: [filter1 smiles, filter1 + filter2 smiles]
		
		# 2. Get major tautomer from jchem:
		taut_obj = Tautomerization()
		taut_obj.postData.update({'calculationType': 'MAJOR'})
		taut_obj.make_data_request(filtered_smiles, taut_obj)

		# todo: verify this is major taut result smiles, not original smiles for major taut request...
		major_taut_smiles = None
		try:
			major_taut_smiles = taut_obj.results['result']['structureData']['structure']
		except KeyError as e:
			# logging.info("Jchem error requesting major tautomer from {}.".format(filtered_smiles))
			# logging.info("Using smiles {} for next step..".format(filtered_smiles))
			pass

		if major_taut_smiles:
			filtered_smiles = major_taut_smiles

		# 3. Using major taut smiles for final "neutralize" filter:
		post_data = {
			'structure': filtered_smiles, 
			'actions': [
				"neutralize"
			]
		}
		response = calc_object.web_call(url, post_data)

		final_smiles = response['results'][-1]

		return final_smiles



	def checkMass(self, chemical):
		"""
		returns true if chemical mass is less
		than 1500 g/mol
		"""
		try:
			json_obj = Calculator().getMass({'chemical': chemical}) # get mass from jchem ws
		except Exception as e:
			logging.warning("!!! Error in checkMass() {} !!!".format(e))
			raise e
		struct_mass = json_obj['data'][0]['mass']

		if struct_mass < 1500  and struct_mass > 0:
			return True
		else:
			return False



	def clearStereos(self, smiles):
		"""
		clears stereoisomers from smiles
		"""
		try:
			response = self.singleFilter({'smiles':smiles, 'action': "clearStereo"})
			filtered_smiles = response['results'] # get stereoless structure
		except Exception as e:
			logging.warning("!!! Error in clearStereos() {} !!!".format(e))
			raise e
		return filtered_smiles



	def transformSMILES(self, smiles):
		"""
		N(=O)=O >> [N+](=O)[O-]
		"""
		try:
			response = self.singleFilter({'smiles':smiles, 'action': "transform"})
			filtered_smiles = response['results'] # get stereoless structure
		except Exception as e:
			logging.warning("!!! Error in transformSMILES() {} !!!".format(e))
			raise e
		return filtered_smiles



	def untransformSMILES(self, smiles):
		"""
		[N+](=O)[O-] >> N(=O)=O
		"""
		try:
			response = self.singleFilter({'smiles':smiles, 'action': "untransform"})
			filtered_smiles = response['results'] # get stereoless structure
		except Exception as e:
			logging.warning("!!! Error in untransformSMILES() {} !!!".format(e))
			raise e
		return filtered_smiles



	def parseSmilesByCalculator(self, structure, calculator):
		"""
		Calculator-dependent SMILES filtering!
		"""
		filtered_smiles = structure

		#1. check structure mass..
		if calculator != 'chemaxon':
			if not self.checkMass(structure):
				# raise "Structure too large, must be < 1500 g/mol.."
				raise Exception({'data': "structure too large"})

		#2-3. clear stereos from structure, untransform [N+](=O)[O-] >> N(=O)=O..
		if calculator == 'epi' or calculator == 'sparc' or calculator == 'measured':
			try:
				# clear stereoisomers:
				filtered_smiles = self.clearStereos(structure)

				# transform structure:
				filtered_smiles = str(filtered_smiles[-1])
				filtered_smiles = str(self.untransformSMILES(filtered_smiles)[-1])
			except Exception as e:
				logging.warning("!!! Error in parseSmilesByCalculator() {} !!!".format(e))
				raise {'data': "error filtering chemical"}

		# 4. Check for metals and stuff (square brackets):
		if calculator == 'epi' or calculator == 'measured':
			if '[' in filtered_smiles or ']' in filtered_smiles:
				# bubble up to calc for handling error
				# raise Exception("{} cannot process metals..".format(calculator))
				raise Exception({'data': "cannot process metals or charges"})

		return filtered_smiles