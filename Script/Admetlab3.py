import json

import requests
import pandas as pd

baseUrl = 'https://admetlab3.scbdd.com'

def transform(data):
    resultList = []
    for mol in data['data']:
        if not mol['data']:
            # Invalid SMILES
            tmp = {'smiles': mol['smiles']}
        else:
            tmp = dict({'smiles': mol['smiles']})
            for _, admet in mol['data'].items():
                for endpoint in admet:
                    # endpoint is a dict
                    tmp[endpoint['name']] = endpoint['value']
        resultList.append(tmp)
    return pd.DataFrame(resultList).fillna('Invalid SMILES')

if __name__ == '__main__':
    api = '/api/admet'
    url = baseUrl + api
    param = {
        'SMILES': []
    }
    data = pd.read_csv('input.csv')  # Read CSV file
    smiles_list = data['SMILES'].tolist()  # Getting data from the SMILES column
    param['SMILES'] = smiles_list  # Assign the SMILES list to the SMILES key of param

    response = requests.post(url, json=param)

    if response.status_code == 200:  # If access is successful
        data = response.json()['data']
        # transform to csv file
        result = transform(data)
        result.to_csv('result.csv', index=False)