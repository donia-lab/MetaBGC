# Script to produce a non-uniform grid heatmap
import matplotlib.pyplot as plt
import matplotlib
import os
import math
import numpy as np
from Bio import SeqIO
import pandas as pd
from matplotlib.colors import ListedColormap
import seaborn as sns
from scipy.cluster import hierarchy
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

def plot_pca(df_breath, sample_breath_matrix, version):
    pca = PCA()
    pipe = Pipeline([('scaler', StandardScaler()), ('pca', pca)])
    X = np.array(sample_breath_matrix)
    Xt = pipe.fit_transform(X)
    plot = plt.scatter(Xt[:,0], Xt[:,1])
    for i, txt in enumerate(df_breath.columns[1:]):
        short_label = txt.split(sep='-')
        plt.annotate('-'.join(short_label[1:]), (Xt[i,0], Xt[i,1]), fontsize=5)
        #plt.annotate(str(i), (Xt[i,0], Xt[i,1]))
    plt.xlabel("PCA_1")
    plt.ylabel("PCA_2")
    pca_plot_file_png = os.path.join(working_dir, 'pca_v' + version + '.png')
    pca_plot_file_eps = os.path.join(working_dir, 'pca_v' + version + '.eps')
    plt.savefig(pca_plot_file_eps)
    plt.savefig(pca_plot_file_png)
    plt.close()



def plot_stacked_bar(df_breath, sample_breath_matrix, version):
    x1 = []
    x2 = []
    y = []
    for idx, sample_col in enumerate(df_breath.columns[1:]):
        gene_pct = sum(sample_breath_matrix[idx]) / len(sample_breath_matrix[idx]) * 100.0
        x1.append(gene_pct)
        x2.append(100-gene_pct)
        y.append(idx)
    x1.sort()
    x2.sort(reverse=True)
    plt.barh(y, x1, color='skyblue')
    plt.barh(y, x2, left=x1, color='gray')
    plt.xlabel("% of Genes")
    plt.ylabel("Genomes")
    plt.xlim([0, 100])
    bar_file_eps = os.path.join(working_dir, 'bar_1_v' + version + '.eps')
    plt.savefig(bar_file_eps)
    plt.close()

    # x = []
    # y = []
    # # Loop over each gene
    # for idx_1 in range(len(sample_breath_matrix[0])):
    #     sum_genomes = 0
    #     for idx_2, sample_col in enumerate(df_breath.columns[1:]):
    #         sum_genomes = sum_genomes + sample_breath_matrix[idx_2][idx_1]
    #     y.append(idx_1)
    #     x.append(sum_genomes / len(df_breath.columns[1:]) * 100.0)
    # x.sort()
    # plt.scatter(x, y, c="blue", s=1)
    # plt.xlabel("genome %")
    # plt.ylabel("genes")
    # scatter_file_eps = os.path.join(working_dir, 'scatter_2_v' + version + '.eps')
    # plt.savefig(scatter_file_eps)
    # plt.close()
    #
    x = []
    y = []
    # Loop over each gene
    for idx_1 in range(len(sample_breath_matrix[0])):
        sum_genomes = 0
        for idx_2, sample_col in enumerate(df_breath.columns[1:]):
            sum_genomes = sum_genomes + sample_breath_matrix[idx_2][idx_1]
        p = sum_genomes / len(df_breath.columns[1:])
        if p == 1 or p == 0:
            e = 0
        else:
            e = -p * math.log2(p) - (1-p) * math.log2(1-p)
        y.append(idx_1)
        x.append(e)
    x.sort()
    plt.hist(x, bins=10)
    plt.xlabel("Variability Entropy")
    plt.ylabel("# of Genes")
    bar_file_eps = os.path.join(working_dir, 'bar_3_v' + version + '.eps')
    plt.savefig(bar_file_eps)
    plt.close()

def plot_percentile(df_breath, sample_breath_matrix, version):
    # x and y are arrays of percentage values
    x = np.arange(0, len(df_breath.columns[1:]), 1) / (len(df_breath.columns[1:])-1)
    y = []
    # Loop over each gene
    for idx_1 in range(len(sample_breath_matrix[0])):
        sum_genomes = 0
        for idx_2, sample_col in enumerate(df_breath.columns[1:]):
            sum_genomes = sum_genomes + sample_breath_matrix[idx_2][idx_1]
        y.append(sum_genomes)

    # create a 1D histogram with 80 bins for each genome
    # Hy has the number of genes that are found in a certain number of genomes
    Hy, _ = np.histogram(y, bins=len(df_breath.columns[1:]))
    Hy_norm = np.flip(np.sum(Hy) - np.cumsum(Hy)) / np.sum(Hy)

    # plot the percentile plot
    # plt.hist2d(x, Hy_norm, bins=80, cmap='gist_heat_r')
    # plt.colorbar()
    # plt.xlabel('% of Genes')
    # plt.ylabel('% of Genomes')
    # plt.xlim([x[0], x[-1]])
    # #plt.xticks(np.arange(x[0], x[-1]+0.1, 0.1))
    # plt.ylim([Hy_norm[0], Hy_norm[-1]])
    # plt.yticks(np.arange(Hy_norm[0], Hy_norm[-1]+0.1, 0.1))
    # plt.show()

    # Plot scatter
    plt.scatter(x, Hy_norm, c="blue", s=5)
    plt.xlabel('% of Genomes')
    plt.ylabel('% of Genes')
    #plt.show()
    percentile_file_eps = os.path.join(working_dir, 'percentile_1_v' + version + '.eps')
    plt.savefig(percentile_file_eps)
    plt.close()


def plot_scatter(df_breath, sample_breath_matrix, version):
    x = []
    y = []
    for idx, sample_col in enumerate(df_breath.columns[1:]):
        x.append(sum(sample_breath_matrix[idx]) / len(sample_breath_matrix[idx]) * 100.0,)
        y.append(idx)
    x.sort()
    plt.scatter(x, y, c ="blue", s=1)
    plt.xlabel("% of Genes")
    plt.ylabel("Genomes")
    plt.xlim([0, 100])
    scatter_file_eps = os.path.join(working_dir, 'scatter_1_v' + version + '.eps')
    plt.savefig(scatter_file_eps)
    plt.close()

    x = []
    y = []
    # Loop over each gene
    for idx_1 in range(len(sample_breath_matrix[0])):
        sum_genomes = 0
        for idx_2, sample_col in enumerate(df_breath.columns[1:]):
            sum_genomes = sum_genomes + sample_breath_matrix[idx_2][idx_1]
        y.append(idx_1)
        x.append(sum_genomes / len(df_breath.columns[1:]) * 100.0)
    x.sort()
    plt.scatter(x, y, c="blue", s=1)
    plt.xlabel("% of Genomes")
    plt.ylabel("Genes")
    scatter_file_eps = os.path.join(working_dir, 'scatter_2_v' + version + '.eps')
    plt.savefig(scatter_file_eps)
    plt.close()

    x = []
    y = []
    # Loop over each gene
    for idx_1 in range(len(sample_breath_matrix[0])):
        sum_genomes = 0
        for idx_2, sample_col in enumerate(df_breath.columns[1:]):
            sum_genomes = sum_genomes + sample_breath_matrix[idx_2][idx_1]
        p = sum_genomes / len(df_breath.columns[1:])
        if p == 1 or p == 0:
            e = 0
        else:
            e = -p * math.log2(p) - (1-p) * math.log2(1-p)
        y.append(idx_1)
        x.append(e)
    x.sort()
    plt.scatter(x, y, c="blue", s=1)
    plt.xlabel("Variability Entropy")
    plt.ylabel("Genes")
    scatter_file_eps = os.path.join(working_dir, 'scatter_3_v' + version + '.eps')
    plt.savefig(scatter_file_eps)
    plt.close()


if __name__ == '__main__':

    working_dir = r'C:\Users\ab50\Documents\data\binning\cEK_Search\pcolormesh'
    seq_file = os.path.join(working_dir, '2716884990_genes.fasta')
    version = "13"
    df_breath_file = os.path.join(working_dir, 'Ga0173608_11_quantified_breath_filtered_LI.csv')
    heatmap_plot_file_png = os.path.join(working_dir, 'Ga0173608_11_quantified_breath_filtered_LI_v' + version + '.png')
    heatmap_plot_file_svg = os.path.join(working_dir, 'Ga0173608_11_quantified_breath_filtered_LI_v' + version + '.svg')
    heatmap_plot_file_pdf = os.path.join(working_dir, 'Ga0173608_11_quantified_breath_filtered_LI_v' + version + '.pdf')
    heatmap_plot_file_eps = os.path.join(working_dir, 'Ga0173608_11_quantified_breath_filtered_LI_v' + version + '.eps')
    dendrogram_plot_file_png = os.path.join(working_dir, 'Ga0173608_11_quantified_breath_filtered_LI_dendrogram_v' + version + '.png')
    dendrogram_plot_file_eps = os.path.join(working_dir, 'Ga0173608_11_quantified_breath_filtered_LI_dendrogram_v' + version + '.eps')

    gene_len_dict = {}
    for record in SeqIO.parse(seq_file, "fasta"):
        desc_str = record.description
        gene_name = desc_str.split(' ')[1]
        gene_len_dict[gene_name] = len(record)

    df_breath = pd.read_csv(df_breath_file)

    # Setup the non-uniform columns based on the length of the genes
    bounds_col = [0]
    prev_sum = 0
    gene_cum_len_dict = {}
    for gene_name in df_breath['locus_tag']:
        bounds_col.append(prev_sum + gene_len_dict[gene_name])
        prev_sum = prev_sum + gene_len_dict[gene_name]
        gene_cum_len_dict[gene_name] = prev_sum

    # Setup the rows based on number of samples
    #bounds_row = [x for x in range((len(df_breath.columns)-1))]
    bounds_row = []
    for x in range((len(df_breath.columns))):
        bounds_row.append(x)
        # Uncomment lines to add horizontal boundary
        #if x != len(df_breath.columns) - 1:
        #    bounds_row.append(x+0.99)
    # Print sample ID list
    with open(os.path.join(working_dir,'sample_ids.csv'), 'w') as f:
        for idx, sample_name in enumerate(df_breath.columns[1:]):
            f.write(str(idx)+','+sample_name+'\n')

    # Setup the list of lists for the breath data
    sample_breath_matrix = []
    for sample_col in df_breath.columns[1:]:
        #sample_breath_matrix.append(list(df_breath[sample_col]))
        # Thresholding
        sample_breath_matrix.append([1 if x >= 0.5 else 0 for x in df_breath[sample_col]])
        # Uncomment to add horizontal boundary
        # sample_breath_matrix.append([0 for x in range(len(df_breath[sample_col]))])

    # Perform PCA Analysis
    #plot_pca(df_breath, sample_breath_matrix, version)

    # Generate scatter plots
    plot_scatter(df_breath, sample_breath_matrix, version)
    plot_stacked_bar(df_breath, sample_breath_matrix, version)

    # Generate percentile plot
    plot_percentile(df_breath, sample_breath_matrix, version)


    sns_grid = sns.clustermap(sample_breath_matrix, col_cluster=False, cbar=False, method='weighted')
    # Reordered based on clustered rows
    clustered_breath_matrix = []
    for row_idx in sns_grid.dendrogram_row.reordered_ind:
        clustered_breath_matrix.append(sample_breath_matrix[row_idx])

    # Setup the alternating pattern for the gene map
    gene_map_matrix = []
    gene_range_list = []
    for x in range(len(df_breath)):
        gene_range_list.append(x % 2)
    gene_map_matrix.append(gene_range_list)

    # define colormap for the sample heatmap and the gene map
    N = 2  # number of desired color bins
    cmap = plt.cm.get_cmap('gray', N)
    cmap2 = ListedColormap(["red", "green"])

    # define the bins and normalize
    bounds = np.linspace(0, 1, N + 1)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    my_dpi = 300
    my_linewidth=my_dpi/(1024*32)
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex='row', figsize=(60, 40), dpi=my_dpi, gridspec_kw={'height_ratios': [60, 1]})
    fig.tight_layout()
    # Sample heatmap
    colormesh1 = ax1.pcolormesh(bounds_col, bounds_row, clustered_breath_matrix, cmap=cmap, norm=norm, linewidths=0.1)
    ax1.set_yticks(np.arange(len(sns_grid.dendrogram_row.reordered_ind)) + 0.5)
    ax1.set_yticklabels(sns_grid.dendrogram_row.reordered_ind)
    #colormesh = ax.pcolormesh(bounds2, bounds1, matrix, cmap=cmap, norm=norm, linewidths=my_linewidth, edgecolor='k', antialiased=False)
    #colormesh.set_edgecolor('black')
    # ax.invert_yaxis()
    #ax1.tick_params(axis='x', which='major', rotation=50)
    ax1.tick_params(labelbottom = False, bottom = False)
    # Gene lengths heatmap
    colormesh2 = ax2.pcolormesh(bounds_col, [0, 1], gene_map_matrix, cmap=cmap2, norm=norm, linewidths=0.1)
    ax2.tick_params(left = False, right = False , labelleft = False ,
                    labelbottom = True, bottom = True)
    #ax.set_xticks(bounds2)
    #ax.set_yticks(bounds1)
    #cbar = fig.colorbar(colormesh, ax=ax)
    #cbar.set_ticks(bounds)
    #ax.plot(x2, x1, color='black', marker='o')
    plt.savefig(heatmap_plot_file_png)
    #plt.savefig(heatmap_plot_file_svg)
    #plt.savefig(heatmap_plot_file_pdf)
    plt.savefig(heatmap_plot_file_eps)
    plt.close()

    # Plot the dendrogram separately
    Z = hierarchy.linkage(sample_breath_matrix, method='weighted', metric='euclidean')
    plt.figure()
    dn = hierarchy.dendrogram(Z)
    hierarchy.set_link_color_palette(['m', 'c', 'y', 'k'])
    fig, axes = plt.subplots(1, 1, figsize=(30, 10), dpi=my_dpi)
    #dn1 = hierarchy.dendrogram(Z, ax=axes[0], above_threshold_color='y',
    #                           orientation='top')
    dn2 = hierarchy.dendrogram(Z, ax=axes,
                               above_threshold_color='#bcbddc',
                               orientation='right')
    hierarchy.set_link_color_palette(None)  # reset to default after use
    plt.savefig(dendrogram_plot_file_png)
    plt.savefig(dendrogram_plot_file_eps)
    plt.close()



