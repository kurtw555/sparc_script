import json
import logging
import requests
import os
from .calculator import Calculator
#from .smilesfilter import SMILESFilter


class SparcCalc(Calculator):
    def __init__(self, smiles=None, melting_point=0.0, pressure=760.0, temperature=25.0):

        Calculator.__init__(self)  # inherit Calculator base class

        #self.base_url = os.environ['CTS_SPARC_SERVER']
        self.base_url = 'https://n2626ugath802.aa.ad.epa.gov'
        self.multiproperty_url = '/sparc-integration/rest/calc/multiProperty'
        self.name = "sparc"
        self.smiles = smiles
        self.solvents = dict()
        self.pressure = pressure
        self.melting_point = melting_point
        self.temperature = temperature
        self.props = ["water_sol", "vapor_press", "henrys_law_con", "mol_diss", "boiling_point"]
        self.sparc_props = {
            "SOLUBILITY": "water_sol",
            "VAPOR_PRESSURE": "vapor_press",
            "WATER_DIFFUSION": "mol_diss",
            "AIR_DIFFUSION": "mol_diss_air",
            "FULL_SPECIATION": "ion_con",
            "HENRYS_CONSTANT": "henrys_law_con",
            "DISTRIBUTION": "kow_no_ph",
            "LOGD": "kow_wph",
            "BOILING_POINT": "boiling_point" 
        }
        self.request_timeout = 10

    def get_sparc_query(self):
        query = {
            'pressure': self.pressure,
            'meltingPoint': self.melting_point,
            'temperature': self.temperature,
            'calculations': self.getCalculations(),
            'smiles': self.smiles,
            'userId': None,
            'apiKey': None,
            'type': 'MULTIPLE_PROPERTY',
            'doSolventInit': False
        }
        return query

    def get_calculation(self, sparc_prop=None, units=None):
        calc = {
            'solvents': [],
            'units': units,
            'pressure': self.pressure,
            'meltingPoint': self.melting_point,
            'temperature': self.temperature,
            'type': sparc_prop
        }
        return calc

    def get_solvent(self, smiles=None, name=None):
        solvent = {
            'solvents': None,
            'smiles': smiles,
            'mixedSolvent': False,
            'name': name
        }
        return solvent


    def getCalculations(self):
        calculations = []
        calculations.append(self.get_calculation("VAPOR_PRESSURE", "Torr"))
        calculations.append(self.get_calculation("BOILING_POINT", "degreesC"))
        calculations.append(self.get_calculation("DIFFUSION", "NO_UNITS"))
        calculations.append(self.get_calculation("VOLUME", "cmCubedPerMole"))
        calculations.append(self.get_calculation("DENSITY", "gPercmCubed"))
        calculations.append(self.get_calculation("POLARIZABLITY", "angCubedPerMolecule"))
        calculations.append(self.get_calculation("INDEX_OF_REFRACTION", "dummy"))

        calcHC = self.get_calculation("HENRYS_CONSTANT", "AtmPerMolPerM3")
        calcHC["solvents"].append(self.get_solvent("O", "water"))
        calculations.append(calcHC)

        calcSol = self.get_calculation("SOLUBILITY", "mgPerL")
        calcSol["solvents"].append(self.get_solvent("O", "water"))
        calculations.append(calcSol)

        calcAct = self.get_calculation("ACTIVITY", "dummy")
        calcAct["solvents"].append(self.get_solvent("O", "water"))
        calculations.append(calcAct)

        calculations.append(self.get_calculation("ELECTRON_AFFINITY", "dummy"))
        calcDist = self.get_calculation("DISTRIBUTION", "NO_UNITS")
        calcDist["solvents"].append(self.get_solvent("O", "water"))
        calcDist["solvents"].append(self.get_solvent("OCCCCCCCC", "octanol"))

        calculations.append(calcDist)

        return calculations


    def data_request_handler(self, request_dict):

        for key, val in self.pchem_request.items():
            if not key in request_dict.keys():
                logging.info("request key {} not in request, using default value: {}".format(key, val))
                request_dict.update({key: val})


        """ _filtered_smiles = ''
        try:
            _filtered_smiles = SMILESFilter().parseSmilesByCalculator(request_dict['chemical'], request_dict['calc']) # call smilesfilter
        except Exception as err:
            logging.warning("Error filtering SMILES: {}".format(err))
            request_dict.update({'data': 'Cannot filter SMILES'})
            return request_dict


        self.smiles = _filtered smiles  # set smiles attribute to filtered smiles
        """
        # Gets melting point for sparc calculations.
        """ if request_dict.get('prop') == 'water_sol' or request_dict.get('prop') == 'vapor_press':                
            self.melting_point = self.get_melting_point(self.smiles, 
                                    request_dict.get('sessionid'), self)
        else:
            self.melting_point = None """

        logging.info("Using melting point: {} for SPARC calculation".format(self.melting_point))

        _response_dict = {}
        for key in request_dict.keys():
            if not key == 'nodes':
                _response_dict[key] = request_dict.get(key)  # fill any overlapping keys from request1

        _response_dict.update({'request_post': request_dict, 'method': None})

        try:
            # Runs ion_con endpoint if it's user's requested property
            if request_dict.get('prop') == 'ion_con':
                response = self.makeCallForPka() # response as d ict returned..
                pka_data = self.getPkaResults(response)
                _response_dict.update({'data': pka_data, 'prop': 'ion_con'})
                return _response_dict

            # Runs kow_wph endpoint if it's user's requested property
            elif request_dict.get('prop') == 'kow_wph':
                response = self.makeCallForLogD() # response as dict returned..
                _response_dict.update({'data': self.getLogDForPH(response, request_dict['ph']), 'prop': 'kow_wph'})
                return _response_dict

            # Runs multiprop request if request prop is not kow_wph or ion_con
            else:
                _post = self.get_sparc_query()
                _url = self.base_url + self.multiproperty_url
                _multi_response = self.makeDataRequest()

                if 'calculationResults' in _multi_response:
                    _multi_response = self.parseMultiPropResponse(_multi_response['calculationResults'], request_dict)
                    return _multi_response

        except Exception as err:
            logging.warning("Exception occurred getting SPARC data: {}".format(err))
            _response_dict.update({
                'data': "request timed out",
                'prop': request_dict.get('prop')
            })
            return _response_dict


    def makeDataRequest(self):
        _post = self.get_sparc_query()
        _url = self.base_url + self.multiproperty_url
        return self.request_logic(_url, _post)


    def request_logic(self, url, post_data):
        """
        Handles retries and validation of responses
        """
        _valid_result = False  # for retry logic
        _retries = 0
        while not _valid_result and _retries < self.max_retries:
            # retry data request to chemaxon server until max retries or a valid result is returned
            try:
                #req = requests.Request(method='POST',url=url,data=json.dumps(post_data), headers=self.headers)
                #prepared = req.prepare()
                #self.pretty_print_POST(req)
                #print(prepared)
                response = requests.post(url, data=json.dumps(post_data), headers=self.headers, timeout=self.request_timeout,verify=False)
                print(response.request.url)
                print(response.request.body)
                print(response.request.headers)
                
                _valid_result = self.validate_response(response)
                if _valid_result:
                    self.results = json.loads(response.content)
                    return self.results
                _retries += 1
            except Exception as e:
                logging.warning("Exception in calculator_sparc.py: {}".format(e))
                _retries += 1
            logging.info("Max retries: {}, Retries left: {}".format(self.max_retries, _retries))
        self.results = "calc server not found"
        return self.results


    def validate_response(self, response):
        """
        Validates sparc response.
        Returns False if data is null, or any other
        values that indicate an error
        """
        if response.status_code != 200:
            logging.warning("sparc server response status: {}".format(response.status_code))
            return False
        
        try:
            response_obj = json.loads(response.content)
        except Exception as e:
            logging.warning("Could not convert response to json object, sparc validate_response: {}".format(e))
            return False

        response_type = response_obj.get('type')  # get SPARC property name

        # prop specific validation:
        if response_type == 'LOGD':
            # check 'plotCoordinates' key, should be list and not None
            if not isinstance(response_obj.get('plotCoordinates'), list):
                # retry if no 'plotCoordinates' key or 'plotCoordinates' is None
                logging.warning("SPARC LOGD 'plotCoordinates' not list as expected...Retrying request...")
                return False

        return True


    def parseMultiPropResponse(self, results, request_dict):
        """
        Loops through data grabbing the results
        and building {calc, prop, data} objects
        for front end.

        TODO: Add more info to returned data (object instead of value)
        """
        if not results or not isinstance(results, list):
            raise Exception("(sparc) no results given to parse or not a list..")

        sparc_response = []
        sparc_response_props = []  # list of multi response props to make sure all props made it
        for prop_data in results:
            sparc_prop = prop_data['type']
            logging.info("sparc prop: {}".format(sparc_prop))
            cts_prop_name = self.sparc_props.get(sparc_prop)
            data = prop_data['result']
            logging.info("cts prop name: {}".format(cts_prop_name))
            if cts_prop_name:
                data_obj = {'calc': 'sparc', 'prop': cts_prop_name, 'data': data}
                logging.info("data obj: {}".format(data_obj))
                sparc_response.append(data_obj)
                sparc_response_props.append(cts_prop_name)

        for prop in request_dict['props']:
            if not prop in sparc_response_props and prop != 'ion_con' and prop != 'kow_wph':
                # if sparc response doesn't have user request prop from multi-response, PANIC!
                logging.info("requested prop {} missing from sparc multi response...".format(prop))

                # add data obj for missing prop with error message for user:
                data_obj = {'calc': 'sparc', 'prop': prop, 'data': "prop not found"}
                sparc_response.append(data_obj)

        return sparc_response


    def makeCallForPka(self):
        """
        Separate call for SPARC pKa
        """
        _pka_url = "/sparc-integration/rest/calc/fullSpeciation"
        _url = self.base_url + _pka_url
        logging.info("URL: {}".format(_url))
        _sparc_post = {
            "type":"FULL_SPECIATION",
            "temperature":25.0,
            "minPh":0,
            "phIncrement":0.5,
            "smiles": str(self.smiles),
            "username":"browser1",
            "elimAcid":[],
            "elimBase":[],
            "considerMethylAsAcid": True
        }
        _post_string = json.dumps(_sparc_post)

        return self.request_logic(_url, _sparc_post)


    def getPkaResults(self, results):
        """
        Gets pKa values from SPARC pKa request
        """
        try:
            pka_results = results['macroPkaResults'] # list of pka..
            pka_data = [] # return list of pka values..
            pka_results_obj = {'pKa': [], 'pKb': []}
            for pka_item in pka_results:
                self.getPkaAndPkbFromResults(pka_item, pka_results_obj)

            if len(pka_results_obj['pKa']) + len(pka_results_obj['pKb']) == 0:
                # pka_results_obj = {'pKa': []}  # Return just "none" if no pKa or pKb
                pka_results_obj = None

        except Exception as e:
            logging.warning("Error getting pka from SPARC response: {}".format(e))
            raise Exception("error parsing sparc request")
        
        # return {"pKa": pka_data}
        return pka_results_obj


    def getPkaAndPkbFromResults(self, pka, pka_obj):
        _pka_type = pka.get('macroPkaType')
        _pka_data = pka.get('macroPka')

        if pka.get('macroPka') != -1000:

            if _pka_type == 'Acid' or _pka_type == 'Both':
                pka_obj['pKa'].append(_pka_data)
            elif _pka_type == 'Base':
                pka_obj['pKb'].append(_pka_data)

        return pka_obj


    def makeCallForLogD(self):
        """
        Seprate call for octanol/water partition
        coefficient with pH (logD?)
        """
        _logd_url = "/sparc-integration/rest/calc/logd"
        _url = self.base_url + _logd_url
        _post = {
           "type":"LOGD",
           "solvent": {
              "solvents": None,
              "smiles": "OCCCCCCCC",
              "mixedSolvent": False,
              "name": "octanol"
           },
           "temperature": 25.0,
           "pH_minimum": 0,
           "pH_increment": 0.1,
           "ionic_strength": 0.0,
           "smiles": self.smiles
        }

        logd_results = self.request_logic(_url, _post)
        return logd_results


    def getLogDForPH(self, results, ph=7.0):
        """
        Gets logD value at ph from
        logD response data
        TODO: add ph functionality
        """
        # logging.info("getting sparc logd at ph: {}".format(ph))
        try:
            plot_data = results['plotCoordinates'] # list of [x,y]..
            for xypair in plot_data:
                if xypair[0] == float(ph):
                    return xypair[1]
        except Exception as e:
            logging.warning("Error getting logD at PH from SPARC: {}".format(e))
            raise

    def pretty_print_POST(self, req):
        """
        At this point it is completely built and ready
        to be fired; it is "prepared".

        However pay attention at the formatting used in 
        this function because it is programmed to be pretty 
        printed and may differ from the actual request.
        """
        print('{}\n{}\r\n{}\r\n\r\n{}'.format(
        '-----------START-----------',
            req.method + ' ' + req.url,
            '\r\n'.join('{}: {}'.format(k, v) for k, v in req.headers.items()),
            req.body,
        ))