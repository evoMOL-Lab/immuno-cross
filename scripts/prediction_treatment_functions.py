import pandas as pd
import os


def create_header_name(names, alleles):
    """
    It receives a column with the names of the proteins and another with the alleles.
    Concatenates the data from both columns, row by row, to generate a header for Fasta sequences.   

    Args:
        names (pandas.core.series.Series): pd.Series object with the names of the proteins
        alleles (pandas.core.series.Series): pd.Series object with the alleles
 
    Returns:
        headers_name_allele (pandas.core.series.Series): pd.Series object with header for Fasta sequences
        Ex.: Major capsid protein (HLA-DRB1*09:01)
    """
    
    headers_name_allele = names + ' ' + '(' + alleles + ')'
    return headers_name_allele

def write_fasta_file(name, prediction, col, path):
    # Check whether the specified path exists or not
    # path = os.path.join("fasta_files")

    fasta_file = prediction[['headers_name_allele', f'{col}']].copy()
    # print(prediction.shape[0])
    isExist = os.path.exists(path)
    if not isExist:
    # Create a new directory because it does not exist 
        os.makedirs(path)
        print("The new directory is created!")

    fasta_file.drop_duplicates(inplace=True)
    # print(fasta_file.shape[0], "\n")
    with open(f"{path}{name}.fasta", 'w') as f:
        f.write("\n".join(">" + fasta_file['headers_name_allele'] + "\n" + fasta_file[f'{col}']))
        

def viral_protein_treatment(prediction):
    """
    Viral proteins can be present in more than one variant. 
    Thus, this function identifies all variants that contain the sequence of certain nonamers.

    Args:
        prediction (pandas.core.frame.DataFrame): MHC-II Binding Predictions of viral proteins

    Returns:
        headers(pandas.core.series.Series): Modified headers containing all variants that have the nonamer sequence
        Ex.: Nonamer FLAGLTLSL of protein BHRF1 is present in the variants B95-8 and GD1. 
             So the header will look like this: BHRF1 (HLA-DRB1*01:01) (B95-8 - GD1)
    """
    unique_nonamers = prediction['Core'].unique()
    groups_nonamers = prediction.groupby('Core')
    headers = pd.Series(dtype='O')

    for nonamer in unique_nonamers:
        group = groups_nonamers.get_group(nonamer).copy()
        variants = " - ".join([variant.replace('strain ', '') for variant in group['variant'].unique()])
        variants = "(" + variants + ")"
        headers = pd.concat([headers, group['headers_name_allele'] + ' ' + variants])
    
    return headers