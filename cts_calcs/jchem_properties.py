import requests
import json
import logging
import os
from .calculator import Calculator


class JchemProperty(Calculator):
    """
    Jchem Webservices Physicochemical Properties API handling.
    This module in particular handles Jchem's "Advanced Properties"
    requests (section 4.2.x in docs) using the /calculate endpoint.
    """

    def __init__(self):

        Calculator.__init__(self)  # inherit calculator base class
        
        self.request_timeout = 20
        self.headers = {'Content-Type': 'application/json'}
        self.max_retries = 3
        self.baseUrl = os.environ['CTS_JCHEM_SERVER']
        self.url_pattern = '/webservices/rest-v0/util/calculate/{}'
        self.results = ''  # json result
        self.name = ''  # name of property
        self.url = ''  # url to jchem ws endpoint
        self.postData = {}  # POST data in json
        self.ph = 7.0



    # def getJchemPropData(self, chemical, prop, ph=7.0, method=None, mass=None):
    def getJchemPropData(self, request_dict):
        """
        Calls jchem web services from chemaxon and
        wraps data in a CTS data object (keys: calc, prop, method, data)
        """
        prop_obj = self.getPropObject(request_dict.get('prop'))
        prop_obj.results = self.make_data_request(request_dict.get('chemical'), prop_obj, request_dict.get('method'))

        prop_obj.results = prop_obj.get_data(request_dict)

        _result_dict = {
            'calc': 'chemaxon',
            'prop': request_dict.get('prop'),
            'data': prop_obj.results
            # 'data': result
        }

        # ADD METHOD KEY:VALUE IF LOGD OR LOGP...
        if request_dict.get('method'):
            _result_dict['method'] = request_dict.get('method')

        return _result_dict



    @classmethod
    def getPropObject(self, prop):
        """
        For getting prop objects
        """
        if prop == 'pKa' or prop == 'ion_con':
            return Pka()
        elif prop == 'isoelectricPoint':
            return IsoelectricPoint()
        elif prop == 'majorMicrospecies':
            return MajorMicrospecies()
        elif prop == 'tautomerization':
            return Tautomerization()
        elif prop == 'stereoisomer':
            return Stereoisomer()
        elif prop == 'solubility' or prop == 'water_sol' or prop == 'water_sol_ph':
            return Solubility()
        elif prop == 'logP' or prop == 'kow_no_ph':
            return LogP()
        elif prop == 'logD' or prop == 'kow_wph':
            return LogD()
        elif prop == 'elementalAnalysis':
            return ElementalAnalysis()
        else:
            raise ValueError("Error initializing jchem property class..")



    def make_data_request(self, structure, prop_obj, method=None):
        url = self.baseUrl + prop_obj.url
        prop_obj.postData.update({
            "result-display": {
                "include": ["structureData", "image"],
                "parameters": {
                    "structureData": "smiles"
                }
            }
        })
        post_data = {
            "structure": structure,
            "parameters": prop_obj.postData
        }

        if method:
            post_data['parameters']['method'] = method

        _valid_result = False  # for retry logic
        _retries = 0
        while not _valid_result and _retries < self.max_retries:
            # retry data request to chemaxon server until max retries or a valid result is returned
            try:
                response = requests.post(url, data=json.dumps(post_data), headers=self.headers, timeout=self.request_timeout)
                _valid_result = self.validate_response(response)
                if _valid_result:
                    prop_obj.results = json.loads(response.content)
                    return json.loads(response.content)
                _retries += 1
            except Exception as e:
                logging.warning("Exception in jchem_calculator.py: {}".format(e))
                _retries += 1

            logging.info("Max retries: {}, Retries left: {}".format(self.max_retries, _retries))
        return None



    def validate_response(self, response):
        """
        Validates jchem response.
        Returns False if data is null, or any other
        values that indicate an error
        """
        if response.status_code != 200:
            logging.warning("cts_celery calculator_chemaxon -- jchem server response: {} \n {}".format(response, response.content))
            return False
        return True



    # def booleanize(self, value):
    #     """
    #     django checkbox comes back as 'on' or 'off',
    #     or True/False depending on version, so this
    #     makes sure they're True/False
    #     """
    #     if value == 'on' or value == 'true':
    #         return True
    #     if value == 'off' or value == 'false':
    #         return False
    #     if isinstance(value, bool):
    #         return value



    def getSpeciationResults(self, jchemResultObjects):
        """
        Loops jchemPropObjects (speciation results) from chemaxon,
        grabs the results and creates an object, jchemDictResults, that's
        used for chemspec_tables and data downloads.
        """
        jchem_results_obj = {}
        for key, value in jchemResultObjects.items():
            
            if not value:
                continue

            if key == 'pKa':
                jchem_results_obj.update({
                    'pka': jchemResultObjects['pKa'].getMostAcidicPka(),
                    'pkb': jchemResultObjects['pKa'].getMostBasicPka(),
                    'pka_parent': jchemResultObjects['pKa'].getParent(),
                    'pka_microspecies': jchemResultObjects['pKa'].getMicrospecies(),
                    'pka_chartdata': jchemResultObjects['pKa'].getChartData()
                })
            elif key == 'majorMicrospecies':
                jchem_results_obj.update({key: jchemResultObjects['majorMicrospecies'].getMajorMicrospecies()})
            elif key == 'isoelectricPoint':
                jchem_results_obj.update({
                    key: jchemResultObjects['isoelectricPoint'].getIsoelectricPoint(),
                    'isopt_chartdata': jchemResultObjects['isoelectricPoint'].getChartData()
                })
            elif key == 'tautomerization':
                jchem_results_obj.update({'tautomers': jchemResultObjects['tautomerization'].getTautomers()})
            elif key == 'stereoisomers':
                jchem_results_obj.update({key: jchemResultObjects['stereoisomers'].getStereoisomers()})

        return jchem_results_obj



class Pka(JchemProperty):
    def __init__(self):
        JchemProperty.__init__(self)
        self.name = 'pKa'
        self.url = self.url_pattern.format('pKa')
        self.postData = {
            "pHLower": 0.0,
            "pHUpper": 14.0,
            "pHStep": 0.1,
            "temperature": 298.0,
            "micro": False,
            "considerTautomerization": True,
            "pKaLowerLimit": 0.0,
            "pKaUpperLimit": 14.0,
            "prefix": "DYNAMIC"
        }

    def getMostAcidicPka(self):
        """
		Picks out pKa acidic value(s), returns list
		"""
        pkaValList = []
        if 'mostAcidic' in self.results:
            for pkaVal in self.results['mostAcidic']:
                pkaValList.append(pkaVal)
            return pkaValList
        else:
            logging.warning("key: 'mostAcidic' not in self.results")
            return pkaValList

    def getMostBasicPka(self):
        """
		Picks out pKa Basic value(s), returns list
		"""
        pkaValList = []
        if 'mostBasic' in self.results:
            for pkaVal in self.results['mostBasic']:
                pkaValList.append(pkaVal)
            return pkaValList
        else:
            logging.warning("no key 'mostBasic' in results")
            return pkaValList

    def getParent(self, test=False):
        """
		Gets parent image from result and adds structure
		info such as formula, iupac, mass, and smiles.
		Returns dict with keys: image, formula, iupac, mass, and smiles
		"""
        try:
            parentDict = {'image': self.results['result']['image']['image'], 'key': 'parent'}
            if not test:
                # Adds additional chem info from jchemws:
                parentDict.update(self.getStructInfo(self.results['result']['structureData']['structure']))
            return parentDict
        except KeyError as ke:
            logging.warning("key error: {}".format(ke))
            return None

    def getMicrospecies(self, test=False):
        """
		Gets microspecies from pKa result and appends 
		structure info (i.e., formula, iupac, mass, and smiles)
		Returns list of microspeceies as dicts
		with keys: image, formula, iupac, mass, and smiles
		"""
        if 'microspecies' in self.results:
            try:
                msList = []
                for ms in self.results['microspecies']:
                    msStructDict = {}  # list element in msList
                    msStructDict.update({'image': ms['image']['image'], 'key': ms['key']})
                    if not test:
                        structInfo = self.getStructInfo(ms['structureData']['structure'])
                        msStructDict.update(structInfo)
                    msList.append(msStructDict)
                return msList
            except KeyError as ke:
                logging.info("> key error: {}".format(ke))
                return None
        else:
            logging.info("no microspecies in results")
            return None

    def getChartData(self):
        if 'chartData' in self.results:
            microDistData = {}  # microspecies distribution data
            for ms in self.results['chartData']:
                valuesList = []  # format: [[ph1,con1], [ph2, con2], ...] per ms
                for vals in ms['values']:
                    xy = []  # [ph1, con1]
                    xy.append(vals['pH'])
                    xy.append(100.0 * vals['concentration'])  # convert to %
                    valuesList.append(xy)
                microDistData.update({ms['key']: valuesList})
            return microDistData
        else:
            return None

    def get_data(self, request_dict, data_type=None):
        """
        Gets pKa data based on type.
        Types:
          +  pchem - most acidic and most basic),
          +  speciation - acidic/basic, major microspecies, isoelectric point
        """

        # Now want it back to pKa and pKb like it used to be:
        pka_values = {
            'pKa': self.getMostAcidicPka(),
            'pKb': self.getMostBasicPka()
        }

        if not pka_values['pKa'] and not pka_values['pKb']:
            # if both lists are empty, return "none"
            pka_values = None

        return pka_values



class IsoelectricPoint(JchemProperty):
    def __init__(self):
        JchemProperty.__init__(self)
        self.name = 'isoelectricPoint'
        self.url = self.url_pattern.format('isoelectricPoint')
        self.postData = {
            "pHStep": 0.1,
            "doublePrecision": 2
        }

    def getIsoelectricPoint(self):
        """
		Returns isoelectricPoint value from results
		"""
        try:
            return self.results['isoelectricPoint']
        except KeyError:
            logging.warning("key 'isoelectricPoint' not in results")
            return None

    def getChartData(self):
        """
		Returns isoelectricPoint chart data
		"""
        valsList = []
        try:
            for pt in self.results['chartData']['values']:
                xyPair = []
                for key, value in pt.items():
                    xyPair.append(value)
                valsList.append(xyPair)
        except KeyError as ke:
            logging.warning("key error: {}".format(ke))
            return valsList
        else:
            return valsList



class MajorMicrospecies(JchemProperty):
    def __init__(self):
        JchemProperty.__init__(self)
        self.name = 'majorMicrospecies'
        self.url = self.url_pattern.format('majorMicrospecies')
        self.postData = {
            "pH": 7.0,
            "takeMajorTautomericForm": True
        }

    def getMajorMicrospecies(self, test=False):
        majorMsDict = {}
        try:
            majorMsDict.update({'image': self.results['result']['image']['image'], 'key': 'majorMS'})
            if not test:
                structInfo = self.getStructInfo(self.results['result']['structureData']['structure'])
                majorMsDict.update(structInfo)  # add smiles, iupac, mass, formula key:values
            return majorMsDict
        except KeyError as ke:
            logging.warning("key error: {}".format(ke))
            return None



class Tautomerization(JchemProperty):
    def __init__(self):
        JchemProperty.__init__(self)
        self.name = 'tautomerization'
        self.url = self.url_pattern.format('tautomerization')
        self.postData = {
            "calculationType": "DOMINANT",
            "maxStructureCount": 1000,
            "considerPH": True,
            "pH": 7.0,
            "enableMaxPathLength": True,
            "maxPathLength": 4,
            "rationalTautomerGenerationMode": False,
            "singleFragmentMode": True,
            "protectAromaticity": True,
            "protectCharge": True,
            "excludeAntiAromaticCompounds": True,
            "protectDoubleBondStereo": False,
            "protectAllTetrahedralStereoCenters": False,
            "protectLabeledTetrahedralStereoCenters": False,
            "protectEsterGroups": True,
            "ringChainTautomerizationAllowed": False
        }

    def getTautomers(self, test=False):
        """
        returns dict w/ key 'tautStructs' and
        value is a list of tautomer images
        """
        tautDict = {'tautStructs': [None]}
        tautImageList = []
        try:

            tauts = self.results['result']  # for DOMINANT tautomers

            for taut in tauts:
                tautStructDict = {'image': taut['image']['image'], 'key': 'taut'}
                if not test:
                    structInfo = self.getStructInfo(taut['structureData']['structure'])
                    tautStructDict.update(structInfo)
                tautStructDict.update({'dist': 100 * round(taut['dominantTautomerDistribution'], 4)})
                tautImageList.append(tautStructDict)

            tautDict.update({'tautStructs': tautImageList})
            return tautImageList
        except KeyError as ke:
            logging.warning("key error: {}".format(ke))
            return None



class Stereoisomer(JchemProperty):
    def __init__(self):
        JchemProperty.__init__(self)
        self.name = 'stereoisomer'
        self.url = self.url_pattern.format('stereoisomer')
        self.postData = {
            "stereoisomerismType": "TETRAHEDRAL",
            "maxStructureCount": 100,
            "protectDoubleBondStereo": False,
            "protectTetrahedralStereo": False,
            "filterInvalid3DStructures": False
        }

    def getStereoisomers(self, test=False):
        stereoList = []
        try:
            for stereo in self.results['result']:
                stereoDict = {'image': stereo['image']['image'], 'key': 'stereo'}
                if not test:
                    structInfo = self.getStructInfo(stereo['structureData']['structure'])
                    stereoDict.update(structInfo)
                stereoList.append(stereoDict)
            return stereoList
        except KeyError as ke:
            logging.warning("key error: {} @ jchem rest".format(ke))
            return None



class Solubility(JchemProperty):
    def __init__(self):
        JchemProperty.__init__(self)
        self.name = 'solubility'
        self.url = self.url_pattern.format('solubility')
        self.postData = {
            "pHLower": 0.0,
            "pHUpper": 14.0,
            "pHStep": 0.1,
            "unit": "MGPERML"
        }

    def getIntrinsicSolubility(self):
        """
		Gets water solubility for chemaxon
		"""
        try:
            return 1000.0 * self.results['intrinsicSolubility']
        except KeyError as ke:
            logging.warning("key error: {}".format(ke))
            return None

    def getPHDependentSolubility(self, ph=7.0):
        """
        Gets ph-dependent water solubility
        """
        try:
            ws_list = self.results['pHDependentSolubility']['values']
            ph = float(ph)
            for item in ws_list:
                item_ph = item['pH']
                item_ws = item['solubility']
                if item_ph == ph:
                    return item_ws
            return "N/A"
        except KeyError as ke:
            logging.warning("key error: {}".format(ke))
            return None

    def convertLogToMGPERL(self, log_val, mass):
        """
        Converts WS values of Log into mg/L.
        Units of mass in mol/g
        """
        try:
            return 1000 * float(mass) * 10**(log_val)
        except TypeError as e:
            logging.warning("exception converting mass to float in jchem_properties: {}".format(e))
            return None

    def get_data(self, request_dict):
        if request_dict.get('prop') == 'water_sol_ph':
            # pH dependent water solubility
            _result = self.getPHDependentSolubility(request_dict.get('ph'))
            # _result = self.convertLogToMGPERL(_result, request_dict.get('mass'))  # (jchem v15.3.16)
            _result = 1000.0 * _result  # converts g/L -> mg/L  # ( jchem v16.10.17)
            return _result
        elif request_dict.get('prop') == 'water_sol':
            _result = self.getIntrinsicSolubility()
            return _result
        else:
            return None



class LogP(JchemProperty):
    def __init__(self):
        JchemProperty.__init__(self)
        self.name = 'logP'
        self.url = self.url_pattern.format('logP')
        self.methods = ['KLOP', 'VG', 'PHYS']
        self.postData = {
            "wVG": 1.0,
            "wKLOP": 1.0,
            "wPHYS": 1.0,
            "Cl": 0.1,
            "NaK": 0.1,
            "considerTautomerization": False
        }

    def getLogP(self):
        """
		Gets pH-independent kow
		"""
        try:
            return self.results['logpnonionic']
        except KeyError as ke:
            logging.warning("ker error: {}".format(ke))
            return None

    def get_data(self, request_dict):
        return self.getLogP()



class LogD(JchemProperty):
    def __init__(self):
        JchemProperty.__init__(self)
        self.name = 'logD'
        self.url = self.url_pattern.format('logD')
        self.methods = ['KLOP', 'VG', 'PHYS']
        self.postData = {
            "pHLower": 0.0,
            "pHUpper": 14.0,
            "pHStep": 0.1,
            "wVG": 1.0,
            "wKLOP": 1.0,
            "wPHYS": 1.0,
            "Cl": 0.1,
            "NaK": 0.1,
            "considerTautomerization": False
        }

    def getLogD(self, ph):
        """
		Gets pH-dependent kow
		"""
        try:
            ph = float(ph)
            chartDataList = self.results['chartData']['values']
            for xyPair in chartDataList:
                if xyPair['pH'] == round(ph, 1):
                    value = xyPair['logD']
                    break
            return value
        except KeyError as ke:
            logging.warning("key error: {}".format(ke))
            return None

    def get_data(self, request_dict):
        return self.getLogD(request_dict.get('ph'))



class ElementalAnalysis(JchemProperty):
    def __init__(self):
        JchemProperty.__init__(self)
        self.name = 'elementalAnalysis'
        self.url = self.url_pattern.format('elementalAnalysis')
        self.methods = None
        self.result_key = 'composition'
        self.postData = {
            "singleFragmentMode": True,
            "symbolID": True
        }

    def get_elemental_analysis(self):
        """
        Returns value from expected result_key for
        elemental analysis.
        """
        try:
            return self.results[self.result_key]
        except KeyError as ke:
            logging.warning("ker error: {}".format(ke))
            return None

    def get_data(self, request_dict):
        return self.get_elemental_analysis()