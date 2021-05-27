"""
Metadata about calculators used by CTS.
"""

class MetaData:
	def __init__(self):
		self.model = ""
		self.collection = ""
		self.modelVersion = ""
		self.description = ""
		self.status = ""
		self.timestamp = ""
		self.url = ""
		self.props = []
		self.availableProps = []

	def create_metadata_object(self, calc_info):
		"""
		Creates metadata object for calculator.
		"""
		




	'metaInfo': {
		'model': "chemaxon",
		'collection': "qed",
		'modelVersion': "Jchem Web Services 15.3.23.0",
		'description': "Cheminformatics software platforms, applications, and services to optimize the value of chemistry information in life science and other R&D.",
		'status': '',
		'timestamp': gen_jid(),
		'url': {
			'type': "application/json",
			'href': "http://qedinternal.epa.gov/cts/rest/chemaxon"
		},
		'props': ['water_sol', 'ion_con', 'kow_no_ph', 'kow_wph'],
		'availableProps': [
			{
				'prop': 'water_sol',
				'units': 'mg/L',
				'description': "water solubility"
			},
			{
				'prop': 'ion_con',
				'description': "pKa and pKa values"
			},
			{
				'prop': 'kow_no_ph',
				'units': "log",
				'description': "Octanol/water partition coefficient",
				'methods': ['KLOP', 'PHYS', 'VG']
			},
			{
				'prop': 'kow_wph',
				'units': "log",
				'description': "pH-dependent octanol/water partition coefficient",
				'methods': ['KLOP', 'PHYS', 'VG']
			}
		]
	}