import os
import requests
from cts_calcs import calculator_sparc


def main():
    calc = 'sparc'
    smiles='CCC'
    melting_point=0.0
    pressure=760.0
    temperature=25.0
    ph = 7.0

    props = ["water_sol", "vapor_press", "henrys_law_con", "mol_diss", "boiling_point"]
    params = {"chemical": smiles, "calc": calc, "props": props, "ph": ph}

    sparc_url = 'https://n2626ugath802.aa.ad.epa.gov/sparc-integration/rest/calc/multiProperty'

    sparc = calculator_sparc.SparcCalc(smiles)
    data = sparc.data_request_handler(params)
    print(data)


if __name__ == '__main__':
    main()
    
