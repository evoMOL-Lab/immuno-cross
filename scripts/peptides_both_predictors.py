import itertools
import string
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

def count_proteins_results(prediction, protein_names, col):
    count_proteins = {}
    proteins_group = prediction.groupby(col)
    for protein in protein_names:
        protein_group = proteins_group.get_group(name=protein)
        count_proteins[protein] = protein_group.shape[0]

    if sum(count_proteins.values()) == prediction.shape[0]:
        return count_proteins
    else:
        raise ValueError


def venn(netMHCII, netMHCIIPan, both_predictors, protein, ax):
    subset = (netMHCII[protein], netMHCIIPan[protein], both_predictors[protein])
    v2 = venn2(subset, set_labels=("netMHCII", "netMHCIIPan"), 
               alpha=1, set_colors=('blue', 'yellow'), ax=ax)
    for text in v2.set_labels:     # the text outside the circle
        text.set_fontsize(14);
    for text in v2.subset_labels:  # the text inside the circle
        text.set_fontsize(12)


def plot_venn(protein_names_sn, netMHCII_sn, netMHCII_pan_sn, both_predictors_sn):
    position = 0
    lines = 4
    cols = 3
    positions = [*itertools.product(range(0, lines), range(0, cols))]

    # A4 canvas
    fig_width_cm = 21                                # A4 page
    fig_height_cm = 29.7
    inches_per_cm = 1 / 2.54                         # Convert cm to inches
    fig_width = fig_width_cm * inches_per_cm         # width in inches
    fig_height = fig_height_cm * inches_per_cm       # height in inches
    fig_size = [fig_width * 1.5, fig_height * 1.5]

    plt.rc('text', usetex=False) # so that LaTeX is not needed when creating a PDF with PdfPages later on
    fig, axs = plt.subplots(lines, cols, figsize=fig_size)

    protein_names = sorted(protein_names_sn.values())
    protein_names.append("EBV Proteins")

    for name in protein_names:
        position_graph = positions[position]
        ax = axs[position_graph]
        venn(netMHCII_sn, netMHCII_pan_sn, both_predictors_sn, name, ax)
        ax.set_box_aspect(4/len(ax.patches))

        ax.text(-0.1, .8, string.ascii_uppercase[position], transform=ax.transAxes, 
                size=16, weight='bold')

        ax.set_title(f"{name}", fontsize=14, y=0.8)
        position += 1

    for delete in range(position, lines * cols):
        fig.delaxes(ax=axs[position[delete]])

    plt.subplots_adjust(left=0.1,
                        bottom=0.15, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.4, # Horizontal 
                        hspace=0.01) # Vertical

    plt.savefig('../figures-tables/figure 3.jpg', dpi=500, bbox_inches='tight')
    plt.show()