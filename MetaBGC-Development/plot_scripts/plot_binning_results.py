import matplotlib.pyplot as plt
import pandas as pd



if __name__ == '__main__':
    read_breath = pd.read_csv('C:\\Users\\ab50\\Documents\\data\\binning\\plotting\\sample_read_breath.csv')
    read_depth = pd.read_csv('C:\\Users\\ab50\\Documents\\data\\binning\\plotting\\sample_read_depth.csv')

    metabat2_breath = pd.read_csv('C:\\Users\\ab50\\Documents\\data\\binning\\plotting\\synthetic_1_metabat2_reference_breath.csv')
    vamb_breath = pd.read_csv('C:\\Users\\ab50\\Documents\\data\\binning\\plotting\\synthetic_1_vamb_reference_breath.csv')

    read_breath_melt = pd.melt(read_breath, id_vars=["reference_genome"], var_name="sample_name", value_name="read_breath")
    read_depth_melt = pd.melt(read_depth, id_vars=["reference_genome"], var_name="sample_name", value_name="read_depth")

    # sample_read_breath = read_breath[['reference_genome','synthetic_1_algae_S0']]
    # sample_read_breath.rename(columns={"synthetic_1_algae_S0": "S0_read_breath"},inplace=True)
    # sample_read_depth = read_depth[['reference_genome','synthetic_1_algae_S0']]
    # sample_read_depth.rename(columns={"synthetic_1_algae_S0": "S0_read_depth"},inplace=True)

    metabat2_breath_melt = pd.melt(metabat2_breath, id_vars=["reference_genome"], var_name="sample_name", value_name="contig_breath")
    metabat2_breath_melt['binning'] = 'MetaBAT2'

    vamb_breath_melt = pd.melt(vamb_breath, id_vars=["reference_genome"], var_name="sample_name", value_name="contig_breath")
    vamb_breath_melt['binning'] = 'VAMB'

    binning_breath = pd.concat([metabat2_breath_melt, vamb_breath_melt])
    binning_breath['contig_breath'] = binning_breath['contig_breath'] / 100.0
    binning_breath['contig_breath'] = binning_breath['contig_breath'].where(binning_breath['contig_breath'] <= 1.0, 1.0)
    #binning_breath['contig_breath'][binning_breath['contig_breath'] > 1.0] = 1.0

    df_plot_breath = pd.merge(binning_breath, read_breath_melt, on=['reference_genome','sample_name'])
    df_plot_depth = pd.merge(binning_breath, read_depth_melt, on=['reference_genome','sample_name'])

    df_plot = pd.merge(df_plot_breath, df_plot_depth, on=['reference_genome', 'sample_name', 'contig_breath', 'binning'])

    #plt.scatter(df_plot_breath.S0_read_breath, df_plot_breath.synthetic_1_algae_S0, s=120, c=df_plot_breath.binning)

    df_plot.hist('read_depth', bins=50)
    plt.yscale('log')
    plt.savefig('C:\\Users\\ab50\\Documents\\data\\binning\\plotting\\read_depth_histo.png',dpi=600,format='png')
    plt.close()

    df_tmp = df_plot[df_plot['read_depth'] <= 50]
    df_tmp.hist('read_depth', bins=50)
    plt.yscale('log')
    plt.savefig('C:\\Users\\ab50\\Documents\\data\\binning\\plotting\\read_depth_histo_50.png',dpi=600,format='png')
    plt.close()

    df_tmp = df_plot.copy()
    df_tmp['read_depth'].values[df_tmp['read_depth'].values > 10] = 10
    df_tmp = df_plot[df_plot['read_depth'] <= 10]
    df_tmp =  df_tmp[df_tmp['binning'] == 'VAMB']
    df_tmp.plot.scatter("read_breath", "contig_breath", c="read_depth", colormap='jet')
    plt.xlabel('Breath Covered by Reads')
    plt.ylabel('Breath Covered by VAMB Binned Contigs')
    plt.savefig('C:\\Users\\ab50\\Documents\\data\\binning\\plotting\\vamb\\VAMB_depth.png',dpi=600,format='png')
    plt.close()

    groups = df_tmp.groupby("sample_name")
    for name, group in groups:
        df_tmp_plot = df_tmp[df_tmp['sample_name'] == name]
        ax = df_tmp_plot.plot.scatter("read_breath", "contig_breath", c="read_depth", colormap='jet')
        plt.xlabel('Breath Covered by Reads')
        plt.ylabel('Breath Covered by Contigs')
        df_tmp_plot[['read_breath','contig_breath','reference_genome']].apply(lambda row: ax.text(*row, size=3),axis=1)
        plt.savefig('C:\\Users\\ab50\\Documents\\data\\binning\\plotting\\vamb\\'+ name +'_depth.png',dpi=600,format='png')
        plt.close()

    df_tmp = df_plot.copy()
    df_tmp['read_depth'].values[df_tmp['read_depth'].values > 10] = 10
    df_tmp =  df_tmp[df_tmp['binning'] == 'MetaBAT2']
    df_tmp.plot.scatter("read_breath", "contig_breath", c="read_depth", colormap='jet')
    plt.xlabel('Breath Covered by Reads')
    plt.ylabel('Breath Covered by MetaBAT2 Binned Contigs')
    plt.savefig('C:\\Users\\ab50\\Documents\\data\\binning\\plotting\\metabat2\\MetaBAT2_depth.png',dpi=600,format='png')
    plt.close()

    groups = df_tmp.groupby("sample_name")
    for name, group in groups:
        df_tmp_plot = df_tmp[df_tmp['sample_name'] == name]
        ax = df_tmp_plot.plot.scatter("read_breath", "contig_breath", c="read_depth", colormap='jet')
        plt.xlabel('Breath Covered by Reads')
        plt.ylabel('Breath Covered by Contigs')
        df_tmp_plot[['read_breath','contig_breath','reference_genome']].apply(lambda row: ax.text(*row, size=3),axis=1)
        plt.savefig('C:\\Users\\ab50\\Documents\\data\\binning\\plotting\\metabat2\\'+ name +'_depth.png',dpi=600,format='png')
        plt.close()





