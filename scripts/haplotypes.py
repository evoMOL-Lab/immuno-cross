import pandas as pd
import regex as re
import itertools
from urllib.request import urlopen
from urllib import parse
from bs4 import BeautifulSoup
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import seaborn as sns
import os


# Load data files
#Table 1 Haplotype and phenotype frequencies of HLA class II alleles](https://link.springer.com/article/10.1007/s00251-011-0513-0/tables/1)
#The code below creates the dataframe directly from the table extracted from l's HTML document.

data = pd.read_html("https://link.springer.com/article/10.1007/s00251-011-0513-0/tables/1")[0]
data.query("Allele != 'Total'", inplace=True)
data.set_index("Locus", inplace=True)

# data.to_csv("/content/drive/MyDrive/mestrado/files/hla_frequencies.csv")
# After saving the table, you can read it directly from the files.
# data = pd.read_csv('/content/drive/MyDrive/mestrado/files/hla_frequencies.csv', index_col=0)

# Format HLA naming
# Modify the nomenclature of HLAs to the conventional format.

#  http://hla.alleles.org/nomenclature/naming.html

def nomenclature(alleles):
  alleles = [re.sub(r'(\*{1}\d{2})', r'\1:', allele) for allele in alleles]
  alleles = [allele.replace('/', '-') if allele[0:2]=='DQ' else allele for allele in alleles]
  return alleles
data['allele_nomenclature'] = nomenclature(data.loc[:,'Allele'])
data

# Separate different HLA class II groups

# Only the DR and DQ HLAs will be analyzed, as they are the main ones related to autoimmune diseases.
# DRB3/4/5 are not listed in the database
# https://onlinelibrary.wiley.com/doi/full/10.1002/iid3.416

path = os.path.join("..", "files", "hla_frequencies.csv")
data.filter(regex=r'D[RQ]', axis=0).to_csv(path)
DR = data.filter(regex='DRB1', axis=0).copy()
DQ = data.filter(regex='DQ', axis=0).copy()

# Generate HLA pairs
# http://www.allelefrequencies.net/hla6003a_scr.asp?hla_selection=A*01-B*08&hla_order=order_3

alleles_DR = DR['allele_nomenclature'].copy()
alleles_DQ = DQ['allele_nomenclature'].copy()

def hla_combiner(alleles_DR, alleles_DQ):
  pairs = list(itertools.product(alleles_DR, alleles_DQ))
  return pd.DataFrame({"pairs":["-".join(pair) for pair in pairs]})

DR_DQ = hla_combiner(alleles_DR, alleles_DQ)

# Generate urls to search haplotypes in the *Allele Frequency Net Database*
# Exemple
# http://www.allelefrequencies.net/hla6003a.asp?hla_selection=DRB1*03:01-DQA1*05:01-DQB1*02:01&hla_order=order_3
# The `hla_selection=DRB1*03:01-DQA1*05:01-DQB1*02:01` snippet defines the alleles that must be present in the haplotype.
# the `&hla_order=order_3` snippet defines that the selection will be shown in order from the most frequent haplotype to the least.
# Font: http://www.allelefrequencies.net/extaccess.asp

def url_create(pair):
  url_pairt_1 = "http://www.allelefrequencies.net/hla6003a.asp?hla_selection="
  url_pairt_3 = "&hla_order=order_3"
  return "{}{}{}".format(url_pairt_1, pair, url_pairt_3)

DR_DQ['Haplotype_Search'] = DR_DQ['pairs'].apply(lambda x: url_create(x))

# Extract data from *The Allele Frequency Net Database*

def scraping(url):
  response = urlopen(url)
  html = response.read()
  soup = BeautifulSoup(html, 'html.pairser')
  return soup

def haplotypes_number(url):
  soup = scraping(url)
  table = soup.find('div', {"id":"divGenNavig"})
  try:
    n_haplotypes = table.findAll('td')[0].text
  except:
    return 0
  n_haplotypes = int(re.search(r'(?<=from )\d*', n_haplotypes).group())
  return n_haplotypes

def haplotype_not_identified(search_result):
  if re.search(r"we did not find any results matching your criteria", search_result.text):
    return True

def format_url_search(page, url_search):
    pattern = r'(asp\?)'
    substitute = r'\1page={}&'.format(page+1)
    return re.sub(pattern, substitute, url_search)

def get_data_table(table):
  haplotypes = []
  frequency = []
  sample_size = []

  for line in table:
    data = line.findAll('td')

    haplotypes.append(data[1].find('a').text)
    frequency.append(data[4].text)
    sample_size.append(data[6].text)

  return pd.DataFrame({
      'haplotypes':haplotypes, 'frequency (%)':frequency, 'sample_size':sample_size
      })

def haplotype_search(lista_urls_pesquisa, pairs):
  haplotypes = pd.DataFrame(
      columns=['pair', 'haplotypes', 'frequency (%)', 'sample_size'])
  for url_search, pair in zip(lista_urls_pesquisa, pairs):

    soup = scraping(url_search)
    search_result = soup.find('div', {"id":"divGenNavig"})

    if haplotype_not_identified(search_result):
      continue

    page_number = search_result.findAll('td')[5].text
    page_number = int(re.search(r'(?<=of )\d*', page_number).group())

    for page in range(1, page_number+1):
      table = soup.find('div', {"id":"divGenDetail"}).find('table', {"class":"tblNormal"}).findAll('tr')[1:]
      result = get_data_table(table)
      result.insert(0, "pair", pair)
      haplotypes = pd.concat([haplotypes, result])
      url_formatted = format_url_search(page, url_search)
      soup = scraping(url_formatted)

  return haplotypes

# It takes a long time to run the following commands. Just run it once and then save the dataframe to avoid losing data. The DR_DQ.csv and haplotypes_DR_DQ.csv files are saved in the "files" folder.

# DR_DQ['haplotypes_number'] = [haplotypes_number(url) for url in DR_DQ['Haplotype_Search']]
# DR_DQ.to_csv("../files/DR_DQ.csv", index=False)
# DR_DQ.to_csv("DR_DQ.csv")
# DR_DQ = pd.read_csv("../files/DR_DQ.csv")
# DR_DQ
# haplotypes_DR_DQ = haplotype_search(DR_DQ.loc[:, 'Haplotype_Search'], DR_DQ.loc[:,'pairs'])
# haplotypes_DR_DQ.to_csv("../files/haplotypes_DR_DQ.csv", index=False)
# haplotypes_DR_DQ = pd.read_csv('/content/drive/MyDrive/mestrado/files/haplotypes_DR_DQ.csv')
# Haplotype study

## Load data
DR_DQ = pd.read_csv("../files/DR_DQ.csv")
DR_DQ.sort_values(by=['pairs'], ascending=False, inplace=True)

path = os.path.join("..", "files", "haplotypes_DR_DQ.csv")
haplotypes_DR_DQ = pd.read_csv(path)

## Haplotypes frequency

# Haplotypes with frequency less than 1% are rare and will be removed.
# https://www.nature.com/articles/s41598-019-42385-6#Sec1
haplotypes_DR_DQ_filter = haplotypes_DR_DQ.query('`frequency (%)` >= 1').copy()

## Number of studies that demonstrated the Haplotype

haplotypes_DR_DQ_by_pais = haplotypes_DR_DQ["pair"].value_counts().to_frame().reset_index()
haplotypes_DR_DQ_by_pais.rename(columns = {'index':'pair', 'pair':'Contagem'}, inplace = True)
haplotypes_DR_DQ_by_pais.sort_values(by = "Contagem", ascending=False, inplace=True)

haplotypes_DR_DQ_by_pais['Selected Pair'] = ["Yes" 
    if Pair in haplotypes_DR_DQ_filter['pair'].unique() 
    else "No" 
    for Pair in haplotypes_DR_DQ_by_pais['pair']]

df1 = haplotypes_DR_DQ_by_pais.query("Contagem > 3").copy()
df2 = haplotypes_DR_DQ_by_pais.query("Contagem <= 3").copy()

def show_values(axs, orient="v", space=.01):
    def _single(ax):
        if orient == "v":
            for p in ax.patches:
                _x = p.get_x() + p.get_width() / 2
                _y = p.get_y() + p.get_height() + (p.get_height()*0.01)
                value = '{:.0f}'.format(p.get_height())
                ax.text(_x, _y, value, ha="center") 
        elif orient == "h":
            for p in ax.patches:
                _x = p.get_x() + p.get_width() + float(space)
                _y = p.get_y() + p.get_height() - (p.get_height()*0.5)
                value = '{:.0f}'.format(p.get_width())
                ax.text(_x, _y, value, ha="left", size=12)

    if isinstance(axs, np.ndarray):
        for idx, ax in np.ndenumerate(axs):
            _single(ax)
    else:
        _single(axs)

mpl.rc('xtick', labelsize=12)
mpl.rc('ytick', labelsize=12) 

fig, ax1 = plt.subplots(figsize=(12,20))

palette = {"Yes":"C0", "No":"yellow"}

ax = sns.barplot(x="Contagem", y="pair", orient="h", palette=palette,  
            hue="Selected Pair", data=haplotypes_DR_DQ_by_pais, ax=ax1, dodge = False)

show_values(axs=ax, orient="h")

ax1.set(xlabel=None)
ax1.set(ylabel=None)

ax1.set_xscale("log")

# https://www.statology.org/seaborn-legend-position/
ax1.legend(bbox_to_anchor=(1.005, 0.9), fontsize=12, title='Selected Pair', title_fontsize=12)

plt.subplots_adjust(left=0.1,
                    bottom=0.05, 
                    right=0.9, 
                    top=0.9, 
                    wspace=1.5, # Horizontal 
                    hspace=0.5) # Vertical

fig.text(0.4, 0.02, 'Total haplotype records per HLA pair (logaritmic scale)', ha='center', size=14)
fig.text(-0.21, 0.5, 'HLA Pair', va='center', rotation='vertical', size=14)
fig.text(0.4, 0.91, 'Total haplotype records per HLA pair.\nPairs that have any record with a frequency greater than or equal to 1% were selected.', ha='center', size=16)

plt.savefig('../figures-tables/figure 2.jpg', 
            dpi=500, bbox_inches='tight')

plt.show()

## Select HLAs to HLA-binding prediction
haplotypes_DR_DQ_filter.sort_values('frequency (%)')
selected_haplotypes = haplotypes_DR_DQ_filter['pair'].drop_duplicates().reset_index(drop=True)

haplotypes_DR_DQ_filter.to_csv("../files/haplotypes_DR_DQ_filter.csv")
selected_haplotypes = '\n'.join(selected_haplotypes)

alleles_DR = DR['allele_nomenclature']
alleles_DR = pd.Series([allele for allele in alleles_DR.values if allele in selected_haplotypes])

alleles_DQ = DQ['allele_nomenclature']
alleles_DQ = pd.Series([allele for allele in alleles_DQ.values if allele in selected_haplotypes])

alleles_DR_IEDB = "HLA-" + alleles_DR.values
alleles_DQ_IEDB = "HLA-" + alleles_DQ.str.replace('-', '/').values
alleles = np.append(alleles_DR_IEDB, alleles_DQ_IEDB) + '\n'
alleles = "".join(alleles)

with open("../files/alleles.txt", 'w') as f:
  f.write(alleles)

with open("../files/selected_haplotypes.txt", 'w') as f:
  f.write(selected_haplotypes)