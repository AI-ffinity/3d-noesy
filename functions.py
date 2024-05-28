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
    return df


def get_n_deviations(df_strong, df_weak):
    """Calcaulates how many interresidual NOEs are
    more intense than the correlation Hn-N-Hnoe

    Arguments:
        - df_strong (DataFrame):  list of intra-residual NOE peaks, which are supposed to be the most intense (filtered by the NOE atom type)
        - df_weak (DataFrame): list of inter-residual list of NOE peaks, which are assumed to be less intense (filtered by the NOE atom type, must be the same as in df_i)
    Returns:
        None
        """
    noes_strong = df_strong[['height', 'res']].groupby('res', as_index=True).max('height')
    noes_weak = df_weak[['height', 'res']].groupby('res', as_index=True).max('height')

    noe_compare = noes_strong.join(noes_weak, how='left',
                                   lsuffix='_noe_intra', rsuffix='_noe_inter').fillna(0)
    n_devs = (noe_compare['height_noe_intra'].apply(np.abs)
              < noe_compare['height_noe_inter'].apply(np.abs)).sum()
    print(n_devs)
    if n_devs == 0:
        print("Brilliant result, we have 0 deviations from the theory!")
    else:
        print("There are some deviations, check those most prominent deviations closer:")
        print(noe_compare.loc[
                  noe_compare['height_noe_intra'].apply(np.abs) < noe_compare['height_noe_inter'].apply(np.abs)])

def get_deviations(df_strong, df_weak):
    """Get the details about the mismatch"""
    pass