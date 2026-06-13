import matplotlib.pyplot as plt
import pandas as pd


if __name__ == '__main__':
    read_breath = pd.read_csv('C:\\Users\\ab50\\Documents\\data\\binning\\plotting\\sample_read_breath.csv')
    metabat2_breath = pd.read_csv('C:\\Users\\ab50\\Documents\\data\\binning\\plotting\\synthetic_1_metabat2_reference_breath.csv')
    vamb_breath = pd.read_csv('C:\\Users\\ab50\\Documents\\data\\binning\\plotting\\synthetic_1_vamb_reference_breath.csv')

    art_read_counts = pd.read_csv('C:\\Users\\ab50\\Documents\\data\\binning\\plotting\\synthetic_read_count.csv', index_col=0, header=None)
    art_read_counts = art_read_counts.transpose()
    art_read_counts.rename(columns={"sample_name": "reference_genome"}, inplace=True)

    art_read_cov = pd.read_csv('C:\\Users\\ab50\\Documents\\data\\binning\\plotting\\synthetic_read_coverage.csv', index_col=0, header=None)
    art_read_cov = art_read_cov.transpose()
    art_read_cov.rename(columns={"sample_name": "reference_genome"}, inplace=True)

    # Melt ops to make everything row based
    read_breath_melt = pd.melt(read_breath, id_vars=["reference_genome"], var_name="sample_name", value_name="read_breath")
    metabat2_breath_melt = pd.melt(metabat2_breath, id_vars=["reference_genome"], var_name="sample_name", value_name="contig_breath_pct")
    vamb_breath_melt = pd.melt(vamb_breath, id_vars=["reference_genome"], var_name="sample_name", value_name="contig_breath_pct")
    art_read_count_melt = pd.melt(art_read_counts, id_vars=["reference_genome"], var_name="sample_name", value_name="art_simulated_read_count")
    art_read_count_melt["art_simulated_read_count"] = art_read_count_melt["art_simulated_read_count"].astype(float)
    art_read_cov_melt = pd.melt(art_read_cov, id_vars=["reference_genome"], var_name="sample_name", value_name="art_simulated_read_coverage")
    art_read_cov_melt["art_simulated_read_coverage"] = art_read_cov_melt["art_simulated_read_coverage"].astype(float)

    # Metabat2 Chimeric assemblies
    df_chimeric_metabat2 = pd.merge(read_breath_melt, metabat2_breath_melt)
    df_chimeric_metabat2 = pd.merge(df_chimeric_metabat2, art_read_count_melt)
    df_chimeric_metabat2 = pd.merge(df_chimeric_metabat2, art_read_cov_melt)
    df_chimeric = df_chimeric_metabat2[(df_chimeric_metabat2['read_breath'] > 0.2) &
                                                 (df_chimeric_metabat2['contig_breath_pct'] > 20.0) &
                                                 (df_chimeric_metabat2['art_simulated_read_count'] == 0)]
    df_chimeric_metabat2.to_csv('C:\\Users\\ab50\\Documents\\data\\binning\\plotting\\metabat2_breaths.csv', index=False)

    # VAMB Chimeric assemblies
    df_chimeric_vamb = pd.merge(read_breath_melt, vamb_breath_melt)
    df_chimeric_vamb = pd.merge(df_chimeric_vamb, art_read_count_melt)
    df_chimeric_vamb = pd.merge(df_chimeric_vamb, art_read_cov_melt)
    df_chimeric = df_chimeric_vamb[(df_chimeric_vamb['read_breath'] > 0.2) &
                                                (df_chimeric_vamb['contig_breath_pct'] > 20) &
                                                (df_chimeric_vamb['art_simulated_read_count'] == 0)]
    df_chimeric_vamb.to_csv('C:\\Users\\ab50\\Documents\\data\\binning\\plotting\\vamb_breaths.csv', index=False)








