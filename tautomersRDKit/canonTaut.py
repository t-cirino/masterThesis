"""To test this script, make sure to be in an environment 
    with RDKit  (2023.09.01 or later) installed"""

import csv
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.MolStandardize.rdMolStandardize import TautomerEnumerator

# Canonical SMILES: This function returns the canonical SMILES string
def canonSmiles(smiles_str) -> str:
    smiles = smiles_str.replace('"','').split(',')[0]
    #From SMILES string to molecular graph and back to canonical SMILES string:
    mol = Chem.MolFromSmiles(smiles)
    canon_smiles = Chem.MolToSmiles(mol, canonical=True)

    return canon_smiles

"""
2. Canonical Tautomer: The functions below are to identify tautomers 
and convert molecules to their "canonical" tautomeric form
"""

taut_enumerator = rdMolStandardize.TautomerEnumerator()

# 2.1 This function returns True if a molecule has multiple tautomeric forms
def isTautomer(smiles_str: str) -> bool:

    molec = Chem.MolFromSmiles(smiles_str.replace('"',''))
    if molec is None:
        return False
    else:
        tauts = taut_enumerator.Enumerate(molec)
        tauts_smiles = [Chem.MolToSmiles(taut, canonical=True) for taut in tauts]

        if len(tauts_smiles) > 1:
            return True
        else:
            return False

#2.2 This function returns the canonical tautomer SMILES string
def canonicTaut(smiles_str: str) -> str:

    molec = Chem.MolFromSmiles(smiles_str.replace('"',''))
    tautomers = taut_enumerator.Enumerate(molec)
    original_smiles = canonSmiles(smiles_str) 

    # to get the best scored tautomeric forms:
    scored_tautomers = [(tautomer, taut_enumerator.ScoreTautomer(tautomer)) for tautomer in tautomers]
    max_score = max(scored_tautomers, key=lambda x: x[1])[1]
    best_tautomers = [tautomer for tautomer, score in scored_tautomers if score == max_score]

    # is True if the input molecule is already among the best scored tautomeric form
    is_original_smiles_best = any(Chem.MolToSmiles(tautomer) == original_smiles for tautomer in best_tautomers)

    # if the molecule is among the best scored, it returns itself, otherwise a canonical form
    if is_original_smiles_best:
        return original_smiles
    else:
        molec = taut_enumerator.Canonicalize(molec)
        canon_taut = Chem.MolToSmiles(molec, canonical=True)   
        return canon_taut

"""
2.3 This function takes a .csv file with SMILES strings and their ID and 
 returns another file with canonical SMILES (but non canonical tautomer) 
 in the 1st col and canonical tautomer in the 3rd col. 4th collumn tells
 if the input SMILES was already canonical.
"""
def canonicTautcsv(input_file, output_file):

    with open (input_file,'r') as infile, open(output_file, 'w', newline='', encoding='utf-8') as outfile:
        next(infile) #skip first line
        data_lst = []
        for line in infile:
            data_lst.append(line.rstrip('\n').split(','))

        writer = csv.writer(outfile)
        writer.writerow(["SMILES","ID","SMILES_RDKit"])
    
        for SMILES,RECORDID in data_lst:

            if isTautomer(SMILES):
                canon_smiles = canonSmiles(SMILES)
                canon_taut = canonicTaut(canon_smiles)

                if canon_smiles != canon_taut:
                    writer.writerow([canon_smiles,RECORDID,canon_taut,"non canonical"])
                else:
                    writer.writerow([canon_smiles,RECORDID,canon_taut,"canonical"])

canonicTautcsv("example_tautomers.csv", "example_output.csv")

