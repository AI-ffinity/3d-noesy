import pandas as pd
import numpy as np


def tidy_list(df):
    """Takes a SPARKY list with intensities and
    reformats it into a nice dataframe
    Arguments:
        - df (DataFrame): the saved sparky peak and intensity list read by pandas
    Returns:
        The nice data frame
    """

    df.drop(columns='Height', inplace=True)
    df.rename({
        'Assignment': 'label',
        'Data': 'height',
        'w1': 'N',
        'w2': 'Hn',
        'w3': 'H'
    }, axis=1, inplace=True)

    df.insert(0, 'noe', df.label.apply(lambda s: s.split('-')[-1]))
    df.insert(0, 'res', df.label.apply(lambda s: s.split('-')[0]))

    df['noe_res'] = df.noe.apply(lambda s: s.split('H')[0])
    df.loc[df.noe_res == '', 'noe_res'] = df.loc[df.noe_res == ''].res
    df['noe_res'] = df.noe_res.str.removesuffix('N')
    df['res'] = df.res.str.removesuffix('N')
    df.drop(columns='label', inplace=True)
    df['inter'] = df.noe_res != df.res

    df['resnum'] = df['res'].str.extract('(\d+)', expand=False).fillna(0).astype(int)
    df['noe_resnum'] = df['noe_res'].str.extract('(\d+)', expand=False).fillna(0).astype(int)

    df['res_diff'] = df['resnum'] - df['noe_resnum']
    df['atom_type'] = df['noe'].str.extract(r'(H[A-Za-z]?)')
    df = df.assign(atom_type_pos=lambda s: s.atom_type + '_i-' + s.res_diff.astype('str'))
    df['atom_type_pos'] = df.atom_type_pos.str.replace('--', '+')
    df['atom_type_pos'] = df.atom_type_pos.str.replace('-0', '')

    return df


def convert_heights_to_relative(df):
    """
    Calculate relative NOE peak intensities for each spin system. Negates the effects of the
    1) Overall dynamic of the residue; 2) Global scaling of signal intensity,
    thereby allowing to include peaks from different spectra in one data set.

    Arguments:
        - df: the result of tidy_list() function above
            (possibly filtered for backbone, diagonal, etc)
    Returns:
        - DataFrame of the same shape as df but with height ranging from 0 to 1 for each residue
    """
    ### df.loc[df.res_diff.abs() > 1, "atom_type_pos"] = df.loc[df.res_diff.abs() > 1, "atom_type"] + "_far"
    max_height = df[['res', 'height']].groupby('res').transform('max')
    # pandas makes strange complaints on the data types here, thus so many explicit type casting
    df.loc[:, 'height'] = df['height'].astype('float')
    df.loc[:, 'height'] = df['height'].astype('float64') / max_height['height'].astype('float64')

    return df



def compare_strongest_noes(df_intra, df_inter):
    """Get two lists of NOEs with inter- and intra-residual peaks and
    returns two height columns
    Arguments:
        - df_strong (DataFrame):  list of intra-residual NOE peaks, which are supposed to be the most intense (filtered by the NOE atom type)
        - df_weak (DataFrame): list of inter-residual list of NOE peaks, which are assumed to be less intense (filtered by the NOE atom type, must be the same as in df_i)
    Returns:
        DataFrame with 4 columns: residue number, height of the intra- and inter-residual NOE and the number of residue which gave that NOE
    """

    idx_inter = df_inter[['height', 'resnum']].groupby(['resnum']).idxmax()

    noes_intra = df_intra.loc[:, ['height', 'resnum']]
    noes_inter = df_inter.loc[idx_inter.height.to_list(), ['height', 'resnum', 'noe_resnum']]

    noes_intra.set_index('resnum', inplace=True)
    noes_inter.set_index('resnum', inplace=True)

    noe_compare = (noes_intra.join(noes_inter, how='outer',
                                   lsuffix='_intra', rsuffix='_inter') \
                   .fillna(0).astype('int'))
    return noe_compare.drop_duplicates() # Duplicates occur in Glycines with Ha1 and Ha2


def get_anomalies(df_strong, df_weak):
    """Get the details about the intensity anomaly"""

    noes_strong = df_strong[['height', 'res']].groupby('res', as_index=True).max('height')
    noes_weak = df_weak[['height', 'res']].groupby('res', as_index=True).max('height')

    noe_compare = noes_strong.join(noes_weak, how='left',
                                   lsuffix='_noe_intra', rsuffix='_noe_inter').fillna(0)

    print("There are some deviations, check those most prominent deviations closer:")
    return (noe_compare.loc[
        noe_compare['height_noe_intra'].apply(np.abs) < noe_compare['height_noe_inter'].apply(np.abs)])

    pass


def get_atoms_w_strongest_noes(df):
    """
    Returns the counts of the atoms giving rise to the strongest NOEs
    with their relative positions included into the atom name.
    Arguments:
        - df: the result of tidy_list() function above
            (possibly filtered for backbone, diagonal, etc)
    :return: DataFrame with columns "Atom type" and "count" (usually no longer than 20 rows)
    """

    idx_strongest_in_spinsys = df[['res', 'height']].groupby('res').idxmax() \
        .height.to_list()

    df_sss = df.loc[idx_strongest_in_spinsys, ['atom_type_pos']]

    result = df_sss.atom_type_pos.value_counts()
    result.index.name = ("Atom type")
    return result.to_frame()


def get_noe_ranks(df, exclude_sc=False):
    """Get the counts for how many times each of the
    H_i-1, HA_i and HA_i-1 appear on
    the 1st-3rd rank of NOE intensity
    Arguments:
        - df: the result of tidy_list() function above
        - exclude_sc: whether only include NOEs to backbone protons (HA and Hn)
                      and exclude sidechains for the rank calculation
        """

    if exclude_sc:
        dfh = df.loc[df['atom_type'].isin(['HA', 'H']), ['res', 'atom_type_pos', 'height']]
    else:
        dfh = df.loc[:, ['res', 'atom_type_pos', 'height']]
    dfh = dfh.drop_duplicates()  # This drops only the GlyHA(2,3) peaks which can not be distinguished

    dfh['rank'] = dfh[['res', 'atom_type_pos', 'height']].groupby(['res'], as_index=False)["height"].rank(
        method='dense', ascending=False)
    dfh_bb = dfh.loc[dfh.atom_type_pos.isin(["HA_i", "HA_i-1", "H_i-1"])]
    rank_counts = dfh_bb.groupby(['atom_type_pos', 'rank']).size().unstack(fill_value=0)
    rank_counts.insert(3, "4+", rank_counts.iloc[:, 3:].sum(axis=1))
    rank_counts = rank_counts.iloc[:, :4]
    rank_counts.index.name = "Atom name"
    rank_counts.columns = ["1st highest", "2nd highest", "3rd highest", "4th or lower"]
    return rank_counts






