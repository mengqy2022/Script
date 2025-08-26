import json
import argparse
import requests
import pandas as pd

baseUrl = 'https://admetlab3.scbdd.com'

def transform(data, symbols):
    resultList = []
    for idx, mol in enumerate(data['data']):
        if not mol['data']:
            # Invalid SMILES
            tmp = {'Symbol': symbols[idx], 'smiles': mol['smiles']}
        else:
            tmp = {'Symbol': symbols[idx], 'smiles': mol['smiles']}
            for _, admet in mol['data'].items():
                for endpoint in admet:
                    # endpoint is a dict
                    tmp[endpoint['name']] = endpoint['value']
        resultList.append(tmp)
    return pd.DataFrame(resultList).fillna('Invalid SMILES')

def main():
    parser = argparse.ArgumentParser(description='Predict ADMET properties from SMILES strings with Symbol IDs')
    
    # Input file argument
    parser.add_argument('-i', '--input', 
                        required=True,
                        help='Input CSV file containing Symbol IDs and SMILES strings')
    
    # Symbol column name argument
    parser.add_argument('-s', '--symbol_column',
                        default='Symbol',
                        help='Column name containing Symbol IDs (default: Symbol)')
    
    # SMILES column name argument
    parser.add_argument('-c', '--smiles_column',
                        default='SMILES',
                        help='Column name containing SMILES strings (default: SMILES)')
    
    # Output file argument
    parser.add_argument('-o', '--output',
                        default='result.csv',
                        help='Output CSV file name (default: result.csv)')
    
    args = parser.parse_args()
    
    api = '/api/admet'
    url = baseUrl + api
    param = {
        'SMILES': []
    }
    
    try:
        data = pd.read_csv(args.input)  # Read CSV file
        symbols = data[args.symbol_column].tolist()  # Get Symbol IDs
        smiles_list = data[args.smiles_column].tolist()  # Get SMILES strings
        param['SMILES'] = smiles_list  # Assign the SMILES list to the SMILES key of param

        response = requests.post(url, json=param)

        if response.status_code == 200:  # If access is successful
            data = response.json()['data']
            # transform to csv file with Symbols
            result = transform(data, symbols)
            result.to_csv(args.output, index=False)
            print(f"Success! Results saved to {args.output}")
        else:
            print(f"Error: API request failed with status code {response.status_code}")
    
    except FileNotFoundError:
        print(f"Error: Input file '{args.input}' not found")
    except KeyError as e:
        print(f"Error: Column '{str(e)}' not found in the input file")
    except Exception as e:
        print(f"An unexpected error occurred: {str(e)}")

if __name__ == '__main__':
    main()
