"""
This code ware designed to handle the result of the NetMHCII and NetMHCIIPan prediction.  
Results from other predictors may not work here and need different treatment process.
"""

import pandas as pd
from Bio import SeqIO
from prediction_treatment_functions import create_header_name,  write_fasta_file, viral_protein_treatment
from peptides_both_predictors import count_proteins_results, plot_venn
import re
from io import StringIO


### Protein names
new_names = {
    'Myelin proteolipid protein':'PLP',
    'Periaxin':'Periaxin',
    'Amyloid beta precursor like protein 1':'APLP1', 
    'Myelin protein P0':'P0 protein', 
    'Hepatocyte cell adhesion molecule':'GlialCam',
    "2',3'-cyclic-nucleotide 3'-phosphodiesterase":'CNP',
    'Myelin-associated glycoprotein':'MAG', 
    'Myelin-oligodendrocyte glycoprotein':'MOG',
    'Myelin basic protein':'MBP',
    'Myelin P2 protein':'P2 protein',
    'Myelin-associated oligodendrocyte basic protein':'MOBP'}
path = "../files/canonical_myelin_protein.fasta"
sequences = SeqIO.parse(path, "fasta")


# This code gets the name of the protein in the fasta sequence
protein_names_sn = {
    re.search(r'(?<=\|).*(?=\|)', sequence.description).group() + 
    '_':new_names[re.search(r'(?<= ).*(?= OS)', sequence.description).group()] 
    for sequence in sequences}
    

ebv_dict = {
    'Apoptosis regulator BHRF1':'BHRF1', 'Epstein-Barr nuclear antigen 3':'EBNA3',
    'Epstein-Barr nuclear antigen 6':'EBNA6', 'Epstein-Barr nuclear antigen 4':'EBNA4',
    'Envelope glycoprotein B':'BALF4', 'Envelope glycoprotein GP350':'BLLF1',
    'Latent membrane protein 1':'LMP1', 'Trans-activator protein BZLF1':'BZLF1',
    'DNA polymerase catalytic subunit':'BALF5',
    'Epstein-Barr nuclear antigen 1':'EBNA1', 'Epstein-Barr nuclear antigen 2':'EBNA2',
    'DNA polymerase processivity factor BMRF1':'BMRF1',
    'Major DNA-binding protein':'BALF2', 'Latent membrane protein 2':'LMP2',
    'Large tegument protein deneddylase':'BPLF1', 'Inner tegument protein':'BOLF1',
    'mRNA export factor ICP27 homolog':'BMLF1', 'Envelope glycoprotein L':'BKRF2',
    'Envelope glycoprotein H':'BXLF2',
    'Replication and transcription activator':'BRLF1', 'Glycoprotein 42':'BZLF2',
    'Nuclear egress protein 2':'NEC2',
    'Ribonucleoside-diphosphate reductase small subunit':'RIR2',
    'Epstein-Barr nuclear antigen leader protein':'EBNA-LP',
    'Secreted protein BARF1':'BARF1', 'Major tegument protein':'BNRF1',
    'Triplex capsid protein 1':'TRX1', 'Major capsid protein':'MCP'}
path = "../files/EBV_proteins.fasta"
sequences = SeqIO.parse(path, "fasta")

# This code gets the name of the protein in the fasta sequence
protein_names_ebv = {'Identity':[],'protein_names':[],'variant':[]}
for sequence in sequences:
    protein_name = (re.search(r'(?<= ).*(?= OS)',sequence.description).group())
    protein_name = ebv_dict[protein_name]
    variant = (re.search(r'(?<=\(strain ).*(?=\))',sequence.description).group())
    codigo = re.search(r'(?<=\|).*(?=\|)', sequence.description).group()
    
    protein_names_ebv['Identity'].append(f"{codigo}_")
    protein_names_ebv['protein_names'].append(protein_name)
    protein_names_ebv['variant'].append(variant)

protein_variants = pd.DataFrame(protein_names_ebv)


### netMHCII prediction treatment
def netMHCII(dir):
    file = open(dir)
    text = file.read()

    pattern = r"(\n<!-- saved from)(.*)(Strong binder threshold   2.00. Weak binder threshold  10.00.\n)"
    text = re.sub(pattern, "", text, flags=re.DOTALL)
    text = re.sub(r"-{3,}","", text)
    text = re.sub(r"(\nAllele:?)(.*?)(\n\n?)", "", text, flags=re.DOTALL)
    text = re.sub(r"(</pre>\n<hr>?)(.*)", "", text, flags=re.DOTALL)
    text = re.sub(r"(?<=\S) +", "\t", text)
    text = re.sub(r"sp_", "", text)

    return text

#### sn prediction
dir = "../prediction/prediction_sn_netMHCII_dq.txt"
dq = netMHCII(dir)

dir = "../prediction/prediction_sn_netMHCII_dr.txt"
dr = netMHCII(dir)

dr_df = StringIO(dr)
dr_df = pd.read_table(dr_df)

dq_df = StringIO(dq)
dq_df = pd.read_table(dq_df)

prediction_sn = pd.concat([dr_df, dq_df]).reset_index(drop=True)
prediction_sn.query("Bind == 'SB' | Bind == 'WB'", inplace=True)
prediction_sn.reset_index(inplace=True, drop=True)
prediction_sn['Identity'] = prediction_sn['Identity'].map(protein_names_sn)

#### ebv prediction
dir = "../prediction/prediction_ebv_netMHCII_dq.txt"
dq = netMHCII(dir)

dir = "../prediction/prediction_ebv_netMHCII_dr.txt"
dr = netMHCII(dir)

dr_df = StringIO(dr)
dr_df = pd.read_table(dr_df)

dq_df = StringIO(dq)
dq_df = pd.read_table(dq_df)

prediction_ebv = pd.concat([dr_df, dq_df]).reset_index(drop=True)
prediction_ebv.query("Bind == 'SB' | Bind == 'WB'", inplace=True)
prediction_ebv.reset_index(inplace=True, drop=True)
prediction_ebv["Bind"].value_counts()
prediction_ebv.reset_index(inplace=True, drop=True)
prediction_ebv = prediction_ebv.merge(protein_variants, on='Identity').copy()

### netMHCIIPan prediction treatment
def netMHCIIPan(dir):
    file = open(dir)
    text = file.read()

    pattern = r"(\n<!-- saved from)(.*)(5%\n)"
    text = re.sub(pattern, "", text,  flags=re.DOTALL)
    text = re.sub(r"-{3,}","", text)
    text = re.sub(r"(Number of strong binders: [0-9]* Number of weak binders: [0-9]*)?", "", text,  flags=re.DOTALL)
    text = re.sub(r"(# Allele: ?)(.*?)(D[RQ]B[1-5]_{0,1}\d{0,4}\n)", "", text,  flags=re.DOTALL)
    text = re.sub(r"(Link to output xls?)(.*)", "", text,  flags=re.DOTALL)
    text = re.sub("(?<=\S) +", "\t", text)
    text = re.sub("sp_", "", text)

    file.close()

    return text


#### sn prediction
dir = "../prediction/prediction_sn_netMHCIIPan_dq.txt"
dq = netMHCIIPan(dir)

dir = "../prediction/prediction_sn_netMHCIIPan_dr.txt"
dr = netMHCIIPan(dir)

dr_df = StringIO(dr)
dr_df = pd.read_table(dr_df)

dq_df = StringIO(dq)
dq_df = pd.read_table(dq_df)

prediction_pan_sn = pd.concat([dr_df, dq_df]).reset_index(drop=True)
prediction_pan_sn.query("BindLevel == '&lt;=SB' | BindLevel == '&lt;=WB'", inplace=True)
prediction_pan_sn.reset_index(inplace=True, drop=True)
prediction_pan_sn["Identity"] = prediction_pan_sn["Identity"].str[0:7]
prediction_pan_sn["Identity"] = prediction_pan_sn["Identity"].map(protein_names_sn)


#### ebv prediction
dir = "../prediction/prediction_ebv_netMHCIIPan_dq.txt"
dq = netMHCIIPan(dir)

dir = "../prediction/prediction_ebv_netMHCIIPan_dr.txt"
dr = netMHCIIPan(dir)

dr_df = StringIO(dr)
dr_df = pd.read_table(dr_df)

dq_df = StringIO(dq)
dq_df = pd.read_table(dq_df)

prediction_pan_ebv = pd.concat([dr_df, dq_df]).reset_index(drop=True)
prediction_pan_ebv.query("BindLevel == '&lt;=SB' | BindLevel == '&lt;=WB'", inplace=True)
prediction_pan_ebv.reset_index(inplace=True, drop=True)
prediction_pan_ebv["Identity"] = prediction_pan_ebv["Identity"].str[0:7]
prediction_pan_ebv = prediction_pan_ebv.merge(protein_variants, on='Identity').copy()


### Peptides predicted by both netMHCII and netMHCIIPan
prediction_sn.query("core == 'VHFFKNIVT'").drop_duplicates(subset="                Allele")
prediction_pan_sn.query("Core == 'VHFFKNIVT'").drop_duplicates("MHC")

columns_pan = ["Of", "Core_Rel", "Score_EL", "%Rank_EL", "Exp_Bind", "Score_BA", "Affinity(nM)", "%Rank_BA"]
prediction_pan_sn.drop(columns=columns_pan, inplace=True)
prediction_pan_ebv.drop(columns=columns_pan, inplace=True)

columns = ["Of", "1-log50k(aff)", "Relia", "%Rank", "Level", "affinity(nM)"]
prediction_ebv.drop(columns=columns, inplace=True)
prediction_sn.drop(columns=columns, inplace=True)

prediction_pan_sn.rename(mapper=str.strip, axis="columns", inplace=True)
prediction_pan_ebv.rename(mapper=str.strip, axis="columns", inplace=True)
prediction_sn.rename(mapper=str.strip, axis="columns", inplace=True)
prediction_ebv.rename(mapper=str.strip, axis="columns", inplace=True)

columns = {"Peptide":"peptide", "MHC":"Allele"}
prediction_pan_sn.rename(mapper=columns, axis="columns", inplace=True)
prediction_pan_ebv.rename(mapper=columns, axis="columns", inplace=True)

prediction_sn = prediction_sn.apply(lambda x: x.str.strip() if x.dtype == "object" else x).copy()
prediction_ebv = prediction_ebv.apply(lambda x: x.str.strip() if x.dtype == "object" else x).copy()

prediction_pan_sn = prediction_pan_sn.apply(lambda x: x.str.strip() if x.dtype == "object" else x).copy()
prediction_pan_ebv = prediction_pan_ebv.apply(lambda x: x.str.strip() if x.dtype == "object" else x).copy()

prediction_sn["id"] = (prediction_sn["peptide"] + " " + prediction_sn["Allele"])
prediction_pan_sn["id"] = (prediction_pan_sn["peptide"] + " " + prediction_pan_sn["Allele"])
peptides = prediction_sn["id"].unique()
both_sn = prediction_pan_sn.query("@peptides in id").reset_index(drop=True).copy()

prediction_ebv["id"] = (prediction_ebv["peptide"] + " " + prediction_ebv["Allele"])
prediction_pan_ebv["id"] = (prediction_pan_ebv["peptide"] + " " + prediction_pan_ebv["Allele"])
peptides = prediction_ebv["id"].unique()
both_ebv = prediction_pan_ebv.query("@peptides in id").reset_index(drop=True).copy()


def formata_hla(hla):
    hla = "HLA-" + (hla.\
        str.strip("HLA-").\
            str.replace("_", "*").\
                str.replace("-", "/").\
                    str.replace("DQA1", "DQA1*").\
                        str.replace("DQB1", "DQB1*").\
                            apply(lambda x: re.sub(r"(?<=\*\d{2}?)", ":", x)))
    return hla

both_ebv.Allele = formata_hla(both_ebv.Allele)
both_sn.Allele = formata_hla(both_sn.Allele)

#### Venn diagram
protein_names = protein_names_sn.values()

netMHCII_pan_sn = count_proteins_results(prediction_pan_sn, protein_names, col = "Identity")
netMHCII_sn = count_proteins_results(prediction_sn, protein_names, col = "Identity")
both_predictors_sn = count_proteins_results(both_sn, protein_names, col = "Identity")

protein_names = set(protein_names_ebv["protein_names"])

netMHCII_sn["EBV Proteins"] = prediction_ebv.shape[0]
netMHCII_pan_sn["EBV Proteins"] = prediction_pan_ebv.shape[0]
both_predictors_sn["EBV Proteins"] = both_ebv.shape[0]

plot_venn(protein_names_sn, netMHCII_sn, netMHCII_pan_sn, both_predictors_sn)

### Fasta file

both_sn['headers_name_allele'] = create_header_name(both_sn.Identity, both_sn.Allele)
path = "../results/15_mer_peptides_fastas/"
for protein in both_sn.Identity.unique():
    filtered_prediction = both_sn.groupby('Identity').get_group(protein).copy()
    write_fasta_file(f"{protein}_sn_peptide", filtered_prediction, 'peptide', path)

path = "../results/nonamers_fastas/"
for protein in both_sn.Identity.unique():
    filtered_prediction = both_sn.groupby('Identity').get_group(protein).copy()
    write_fasta_file(f"{protein}_sn_core", filtered_prediction, 'Core', path)

both_ebv['headers_name_allele'] = create_header_name(both_ebv.protein_names, both_ebv.Allele)
both_ebv['headers_name_allele'] = viral_protein_treatment(both_ebv)
write_fasta_file("ebv_core", both_ebv, 'Core', path)

both_ebv.to_csv("../prediction/ebv_both_predictors.csv", index=False)
both_sn.to_csv("../prediction/sn_both_predictors.csv", index=False)