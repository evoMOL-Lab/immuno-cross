from xlsxwriter.workbook import Workbook
from Bio import pairwise2 as pw
from Bio.Align import substitution_matrices as sm
import glob, os, itertools
import regex as re


def create_excel_file(protein):
    global worksheet
    """ Create workbook and worksheet objects for a protein of interest

    Args:
        protein ('_io.TextIOWrapper'): File containing the sequences of the protein of interest
    Returns:
        Workbook and worksheet objects
    """
    path = os.path.join("..", "IdentityAnalysis")
    isExist = os.path.exists(path)
    if not isExist:
    # Create a new directory because it does not exist 
        os.makedirs(path)
        print("The new directory is created!")

    protein_name = re.search(r'.*(?=\.fasta)', protein.name).group()
    dir_name = os.path.join("..", "IdentityAnalysis", f"{protein_name}.xlsx")
    workbook = Workbook(dir_name)
    worksheet = workbook.add_worksheet()

    return workbook, worksheet


def style_text(workbook):
    global style1, style2
    """
    Defines a formatting style to highlight the identity between two amino acid residues.
    Two equal residues will be formatted with 'vivid blue' (#058ED9) color, bold and font size 14.
    Two different residues will be formatted with black color and font size 12.

    Args:
        workbook: xlsxwriter.workbook.Workbook type object

    Returns:
        Objects xlsxwriter.format.Format: style1 and style2
    """
    style1 = workbook.add_format({'color': '#058ED9', 'bold': True, 'font_size': 14, 'font':'arial'})
    style2 = workbook.add_format({'color': 'black', 'font_size': 12, 'font':'arial'})
    return style1, style2


def get_nonamers_and_headers(file):
    """Get amino acid sequences (nonamers) and headers from a fasta file

    Args:
        file (_io.TextIOWrapper): Fasta file containing sequences of the protein(s) of interest

    Returns:
        Lists containing the nonamers and headers of a fasta file
    """
    index = 1

    headers = []
    nonamers = []
    for line in file:
        line = line.rstrip("\n")
        """
        When index is odd, the corresponding line of text must be added to the header list. 
        Otherwise, it is added to the nonamers list.
        """
        if index % 2 == 1:
            # The first character of the header is the '>', so the content is stored from position 1 (second character)
            headers.append(line[1:])
        else:
            nonamers.append(line)
        index += 1

    return nonamers, headers


def select_rule(identity):
    """
    Identifies rule for which the pair of nonamers was selected.

    General rule: An identical residue in P5
    Rule 1: All residues at positions P2, P3 and P8 are identical.
    Rule 2: Two identical residues in P2, P3 or P8 and at least one in P4, P6 or P7.
    Rule 3: All residues in positions P4, P6 and P7 are identical and at least one in P2, P3 ou P8.

    Args:
        identity (int): identity between two nonamers
    Returns:
        identity (str): rule for which pair of nonamers was selected
    """
    rules = {30:'Rule 1', 33:'Rule 1', 36:'Rule 1', 39:'Rule 1', 
             28:'Rule 2', 31:'Rule 2', 34:'Rule 2', 29:'Rule 3'}
    rule = rules[identity]
    return rule


def identity_TCR_contacts(nonamer_1, nonamer_2):
    """
    Checks if there is identity between two residues located in the same position. 
    The residue is formatted with style 1 if identical, or style 2 if not identical.

    Args:
        nonamer_1 (str): nonamer from the protein of interesting 1
        nonamer_2 (str): nonamer from the protein of interesting 2
        style1 (xlsxwriter.format.Format): Formatting style for identical residues ('vivid blue' color, bold and font size 14)
        style2 (xlsxwriter.format.Format): Formatting style for different residues (black color and font size 12)

    Returns:
        Lists containing the stylized residues of nonamers from the proteins of interest 1 (formatting_style_1) and 2 (formatting_style_2).
        It also returns the value of the identity between the two nonamers.
    """
    formatting_style_1 = []
    formatting_style_2 = []
    identity = 0

    # The identity between residues is scored according to the relevance of their position for recognition by the TCR

    positions = {'P1': 0.0, 'P2': 5, 'P3': 5, 'P4': 3, 'P5': 15, 'P6': 3, 'P7': 3, 'P8': 5, 'P9': 0.0}

    values = list(positions.values())
    position = 0

    # Nonamer_1 and nonamer_2 are strings, so the zip function creates tuples containing characters located in the same position.
    # Ex.: zip('ABC', 'DEF') -> (('A','D'), ('B','E'), ('C','F'))

    for aa_1, aa_2 in zip(nonamer_1, nonamer_2):
        if aa_1 == aa_2:
            formatting_style_2.extend((style1, aa_2))
            formatting_style_1.extend((style1, aa_1))
            identity += values[position]

        else:
            formatting_style_2.extend((style2, aa_2))
            formatting_style_1.extend((style2, aa_1))

        position += 1

    return formatting_style_1, formatting_style_2, identity


def hla_analysis(haplotype_list, hla_1, hla_2):
    """ Checks if the two alleles are the same or if they are present in the same haplotype

    Args:
        haplotype_list (list): list containing the HLA-DQ and HLA-DR pairs that are present in the same haplotype haplotype
        hla_1 (str): HLA for which the nonamer of the protein of interest 1 was predicted
        hla_2 (str): HLA for which the nonamer of the protein of interest 2 was predicted

    Returns:
        (bool): Returns True if the two alleles are the same or if they are present in the same haplotype. Otherwise, it returns False.
    """
    same_alleles = (hla_1 == hla_2)

    if same_alleles:
        return True

    for haplotype in haplotype_list:
        allele_1_present = hla_1 in haplotype
        allele_2_present = hla_2 in haplotype

        if allele_1_present and allele_2_present:
            return True

    return False


def selected_haplotypes():
    """
    Load file containing the HLA-DQ and HLA-DR pairs that are present in haplotypes

    Returns:
        haplotype_list (list): list containing the HLA-DQ and HLA-DR pairs that form haplotypes
    """
    dir_name = os.path.join("..", "..", "files", "selected_haplotypes.txt")
    file_haplotypes = open(dir_name)
    haplotype_list = file_haplotypes.read().split('\n')
    file_haplotypes.close()
    return haplotype_list


def extract_headers_data(headers):
    """
    Find the HLA present in the nonamer header

    Args:
        headers (str): Header in the following format: 'protein name (HLA allele)'
        Ex.: CNP (HLA-DRB1*04:01)

    Returns:
        hlas(list): List of HLA alleles
        new_headers(list): List of new headers
    """
 
    # The code below finds the HLA in the header. If it is HLA-DQ, the '/' is replaced by '-'
    # Ex.: CNP (HLA-DRB1*04:01) -> DRB1*04:01
    # Ex.: MAG (HLA-DQA1*05:01/DQB1*02:01) -> DQA1*05:01-DQB1*02:01

    pattern_hla_find = re.compile(r'(?<=\(HLA-).*\d{2}:\d{2}(?=\))')
    hlas = [re.search(pattern_hla_find, headers).group().replace("/", "-") 
            if re.search(pattern_hla_find, headers)
            else ''
            for headers in headers]

    # The code below finds the HLA in the header and replaces it with nothing
    # Ex.: 'CNP (HLA-DRB1*04:01)' > CNP
    # Ex.: 'BHRF1 (HLA-DRB1*12:01) (B95-8 - GD1)' > 'BHRF1  (B95-8 - GD1)'

    pattern_hla_sub = re.compile(r'\(HLA-.*\d{2}:\d{2}\)')
    new_headers = [re.sub(pattern_hla_sub, '', headers) for headers in headers]
    
    return hlas, new_headers


def combines_alleles_and_nonamers(headers_1, headers_2, nonamers_1, nonamers_2):
    """ Receives headers and nonamers from the proteins being analyzed.
        Generates all possible combinations of pairs between the nonamers of the two proteins.
        Generates identifiers for each nonamers and the respective HLA for which it was predicted.

    Args:
        headers_1 (list): List containing the headers of the protein of interest 1.
        headers_2 (list): List containing the headers of the protein of interest 2.
        nonamers_1 (list): List containing the nonamers of the protein of interest 1.
        nonamers_2 (list): List containing the nonamers of the protein of interest 2.

    Returns:
        hlas_nonamers_combination (list): List containing all possible combinations between the nonamers of the proteins being analyzed.\n
        See the example below to understand what this variable represents:

            Nonamer from protein APLP1 presented by HLA-DRB1*01:01
                >APLP1 (HLA-DRB1*01:01)
                LLLLRAQPA

            Five nonamers from Epstein-Barr Virus
                >BHRF1 (HLA-DRB1*12:01) (B95-8 - GD1)
                FLAGLTLSL\n
                >BHRF1 (HLA-DRB1*11:01) (B95-8 - GD1)
                LLALCIRDS\n
                >BHRF1 (HLA-DRB1*11:01) (B95-8 - GD1)
                YLFISRGRH\n
                >BHRF1 (HLA-DRB4*01:01) (B95-8 - GD1)
                YHVLLEEII\n
                >BHRF1 (HLA-DRB4*01:01) (B95-8 - GD1)
                TLSLLVICS\n

            The HLA DRB1*01:01 of APLP1 is combined with all five HLAs from EBV\n
            HLA_pairs = [('DRB1*01:01', 'DRB1*12:01'),\n
                         ('DRB1*01:01', 'DRB1*11:01'),\n
                         ('DRB1*01:01', 'DRB1*11:01'),\n
                         ('DRB1*01:01', 'DRB4*01:01'),\n
                         ('DRB1*01:01', 'DRB4*01:01')]
            
            The nonamer LLLLRAQPA of APLP1 is combined with all five nonamers from EBV\n
            nonamers_pairs = [('LLLLRAQPA', 'FLAGLTLSL'),\n
                             ('LLLLRAQPA', 'LLALCIRDS'),\n
                             ('LLLLRAQPA', 'YLFISRGRH'),\n
                             ('LLLLRAQPA', 'YHVLLEEII'),\n
                             ('LLLLRAQPA', 'TLSLLVICS')]
            
            
            Each tuple in the list contains the nonamers and HLAs for which they were predicted. 
            These nonamers will be evaluated for the identity of their residues.\n
            hlas_nonamers_combination = [('DRB1*01:01', 'DRB1*12:01', 'LLLLRAQPA', 'FLAGLTLSL'),\n
                                         ('DRB1*01:01', 'DRB1*11:01', 'LLLLRAQPA', 'LLALCIRDS'),\n
                                         ('DRB1*01:01', 'DRB1*11:01', 'LLLLRAQPA', 'YLFISRGRH'),\n
                                         ('DRB1*01:01', 'DRB4*01:01', 'LLLLRAQPA', 'YHVLLEEII'),\n
                                         ('DRB1*01:01', 'DRB4*01:01', 'LLLLRAQPA', 'TLSLLVICS')]
        
        identifier_1 (dict): Dictionary containing the header keys and the values that identify them.
                                Ex.: Header -> CNP (HLA-DRB1*04:01)
                                     Nonamer -> FRKMSSSGA
                                     HLA -> DRB1*04:01
                                     Identifier -> DRB1*04:01-FRKMSSSGA
        identifier_2 (dict): Dictionary containing the header keys and the values that identify them.
    """

    HLAs_1, new_headers_1 = extract_headers_data(headers_1)
    HLAs_2, new_headers_2 = extract_headers_data(headers_2)

    HLA_pairs = itertools.product(HLAs_1, HLAs_2)
    nonamers_pairs = itertools.product(nonamers_1, nonamers_2)
    hlas_nonamers_combination = [
        alleles + nonamers for alleles, nonamers in (zip(HLA_pairs, nonamers_pairs))]

    identifier_1 = ["-".join([h, n])
                        for h, n in [*zip(HLAs_1, nonamers_1)]]
    identifier_1 = dict(zip(identifier_1, new_headers_1))

    identifier_2 = ["-".join([h, n])
                       for h, n in [*zip(HLAs_2, nonamers_2)]]
    identifier_2 = dict(zip(identifier_2, new_headers_2))

    return hlas_nonamers_combination, identifier_1, identifier_2


def calculate_similarity_score(nonamer_1, nonamer_2):
    """
    Calculates the blosum62 similarity score between the central region of the nonamers

    Args:
        nonamer_1 (str): Nonamer from protein of interest 1
        nonamer_2 (str): Nonamer from protein of interest 2

    Returns:
        similarity_score (int): Similarity score of the central region of the nonamers (P2 to P8)
    """
    blosum62 = sm.load("BLOSUM62")
    aln = pw.align.globalds(
        nonamer_1[1:8], nonamer_2[1:8], blosum62, -5, -1)[0]
    result = pw.format_alignment(*aln)
    similarity_score = int(re.search('\d+', result).group())
    return similarity_score


def write_line_in_excel(line, formatting_style_1, formatting_style_2, identifier_1, identifier_2, HLA_1, 
                        HLA_2, nonamer_1, nonamer_2, similarity_score, selection_rule):
    """
    Write a line in the worksheet object

    Args:
        worksheet (worksheet): Worksheet object in which the lines with the result of the identity analysis between the nonamers will be written
        line (int): Number that indicates the position of the line that should receive the results
        formatting_style_1 (list): List containing the formatted characters of nonamers
        formatting_style_2 (list): List containing the formatted characters of nonamers
        identifier_1 (dict): Dictionary containing the headers in the keys and the values that identify them.
        identifier_2 (dict): Dictionary containing the headers in the keys and the values that identify them.
        HLA_1 (str): HLA allele
        HLA_2 (str): HLA allele
        nonamer_1 (str): Nonamers from protein of interest 1
        nonamer_2 (str): Nonamers from protein of interest 2
        similarity_score (int): Similarity score between nonamers
    """
    identifier_1 = identifier_1[f'{HLA_1}-{nonamer_1}']
    identifier_2 = identifier_2[f'{HLA_2}-{nonamer_2}']

    # Defines the texts (the column names) of the first eight cells in the first row.
    columns = ['Protein_1', 'HLA_1', 'Nonamer_1', 'Nonamer_2', 'HLA_2',
               'Protein_2', 'Similarity Score BLOSUM62', 'Selection Rule']
               
    # Write column names from cell A1
    worksheet.write_row('A1', columns, style2)

    # Write the results of each analysis.
    worksheet.write(line, 0, identifier_1, style2)
    worksheet.write(line, 1, HLA_1, style2)
    worksheet.write_rich_string(line, 2, *formatting_style_1)
    worksheet.write_rich_string(line, 3, *formatting_style_2)
    worksheet.write(line, 4, HLA_2, style2)
    worksheet.write(line, 5, identifier_2, style2)
    worksheet.write(line, 6, similarity_score, style2)
    worksheet.write(line, 7, selection_rule, style2)


def nonamer_identity_analysis(file_1, file_2):
    """
    Function responsible for calling the other functions necessary for the analysis of identity between nonamers.

    Args:
        file_1 (_io.TextIOWrapper): File containing the nonamers of the protein of interest 1
        file_2 (_io.TextIOWrapper): File containing the nonamers of the protein of interest 2
        worksheet (worksheet): Object on which the analysis results will be written
        style1 (format): Formatting style for identical residues
        style2 (format): Formatting style for different  residues
    """
    line_excel = 1
    nonamers_1, headers_1 = get_nonamers_and_headers(file_1)
    nonamers_2, headers_2 = get_nonamers_and_headers(file_2)
    hlas_nonamers_combination, identifier_1, identifier_2 = combines_alleles_and_nonamers(
        headers_1, headers_2, nonamers_1, nonamers_2)

    haplotype_list = selected_haplotypes()

    for hla_1, hla_2, nonamer_1, nonamer_2 in hlas_nonamers_combination:
        equal_hlas_or_haplotype = hla_analysis(haplotype_list, hla_1, hla_2)

        if equal_hlas_or_haplotype:
            formatting_style_1, formatting_style_2, identity = identity_TCR_contacts(nonamer_1, nonamer_2)
        else:
            continue

        if identity >= 28:
            selection_rule = select_rule(identity)
            similarity_score = calculate_similarity_score(nonamer_1, nonamer_2)
            write_line_in_excel(line_excel, formatting_style_1, formatting_style_2, identifier_1, identifier_2, hla_1, 
                                hla_2, nonamer_1, nonamer_2, similarity_score, selection_rule)
            line_excel += 1


def main(peptides_1, peptides_2):
    """
    Function responsible for reading the files with the nonamers 
    and running the create_excel_file, style_text and analyze_identity_nonamers functions.

    Args:
        peptides_1 (str): Name of the fasta file with the nonamers of the protein of interest 1.
        peptides_2 (str): Name of the fasta file with the nonamers of the protein of interest 2.
    """
    file_1 = open(peptides_1)
    file_2 = open(peptides_2)
    workbook, worksheet = create_excel_file(file_1)
    style1, style2 = style_text(workbook)

    nonamer_identity_analysis(file_1, file_2)

    workbook.close()
    file_1.close()
    file_2.close()

try:
    path = os.path.join("..", "results", "nonamers_fastas")
    os.chdir(path)
except:
    pass

for proteina in glob.glob(r"*_sn_core.fasta"):
    print(proteina)
    main(proteina, 'ebv_core.fasta')