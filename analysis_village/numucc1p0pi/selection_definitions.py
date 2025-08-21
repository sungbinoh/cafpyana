import numpy as np
import pandas as pd


def InFV(data): # cm
    xmin = -190.
    ymin = -190.
    zmin = 10.
    xmax = 190.
    ymax =  190.
    zmax =  450.
    return (np.abs(data.x) > 10) & (np.abs(data.x) < 190) & (data.y > ymin) & (data.y < ymax) & (data.z > zmin) & (data.z < zmax)


# --- FV recommendation from TPC studies by Sungbin ---
def InFV_nohiyz(data):
    xmin = 10.
    xmax = 190.
    zmin = 10.
    zmax = 450.
    ymax_highz = 100.
    pass_xz = (np.abs(data.x) > xmin) & (np.abs(data.x) < xmax) & (data.z > zmin) & (data.z < zmax)
    pass_y = ((data.z < 250) & (np.abs(data.y) < 190.)) | ((data.z > 250) & (data.y > -190.) & (data.y < ymax_highz))
    return pass_xz & pass_y

def InFV_nohiyz_trk(data):
    xmax = 190.
    zmin = 10.
    zmax = 450.
    ymax_highz = 100.
    pass_xz = (np.abs(data.x) < xmax) & (data.z > zmin) & (data.z < zmax)
    pass_y = ((data.z < 250) & (np.abs(data.y) < 190.)) | ((data.z > 250) & (data.y > -190.) & (data.y < ymax_highz))
    return pass_xz & pass_y


# --- for event topoloy breakdown ---
def IsNu(df):
    return ~df.pdg.isna()


def IsSignal(df): # definition                                                                                                                                                                                                                                                                         
    is_fv = InFV(df.position)
    is_1mu1p0pi = (df.nmu_27MeV == 1) & (df.npi_30MeV == 0) & (df.np_50MeV == 1) & (df.npi0 == 0) & (df.mu.genE < 1.2) # & (df.np_20MeV == 1) : add with stubs
    return is_fv & is_1mu1p0pi


def Is1muNp0pi(df): # definition                                                                                                                                                                                                                                                                         
    is_fv = InFV(df.position)
    is_1mu1p0pi = (df.nmu_27MeV == 1) & (df.npi_30MeV == 0) & (df.np_50MeV > 1) & (df.npi0 == 0) #& (df.mu.genE > 0.25) # & (df.np_20MeV == 1) : add with stubs
    return is_fv & is_1mu1p0pi


def Is1muNcpi(df): # definition                                                                                                                                                                                                                                                                         
    is_fv = InFV(df.position)
    is_1mu1p0pi = (df.nmu_27MeV == 1) & (df.npi_30MeV > 0) & (df.npi0 == 0) #& (df.mu.genE > 0.25) # & (df.np_20MeV == 1) : add with stubs
    return is_fv & is_1mu1p0pi


def get_int_category(df):
    is_notnu = ~IsNu(df)
    is_nu_outfv = IsNu(df) & ~InFV(df.position)
    is_signal = IsSignal(df)
    is_1muNp0pi = Is1muNp0pi(df)
    is_1muNcpi = Is1muNcpi(df)

    # assert there's no overlap between the categories, just in case i messed something up...
    assert (is_1muNp0pi & is_1muNcpi).sum() == 0
    assert (is_1muNp0pi & is_signal).sum() == 0
    assert (is_1muNcpi & is_signal).sum() == 0
    assert (is_1muNp0pi & is_nu_outfv).sum() == 0
    assert (is_1muNcpi & is_nu_outfv).sum() == 0
    assert (is_signal & is_nu_outfv).sum() == 0
    is_other_nu_infv = IsNu(df) & InFV(df.position) & ~IsSignal(df) & ~Is1muNp0pi(df) & ~Is1muNcpi(df)

    nuint_categ = pd.Series(8, index=df.index)
    nuint_categ[is_notnu] = -1  # not nu
    nuint_categ[is_nu_outfv] = 0  # nu out of FV
    nuint_categ[is_signal] = 1    # nu in FV, signal
    nuint_categ[is_1muNp0pi] = 2  # 1mu, 0cpi, Np, 0pi0
    nuint_categ[is_1muNcpi] = 3  # 1mu, Ncpi, 0pi0
    nuint_categ[is_other_nu_infv] = 4  # nu in FV, not signal

    return nuint_categ


# --- for plotting ---
# signal / backgroundtopology breakdown
# the signal mode code MUST be the first item in the list for all the code below to work
mode_list = [1, 2, 3, 4, 0, -1]
mode_labels = ["Signal", r"$1\mu N(>1)p 0\pi^{\pm}$", r"$1\mu Ncpi$", "Other FV Nu", "Non FV Nu", "Not Nu"]

colors = ["mediumslateblue",
              "darkslateblue",
              "darkgreen",
              "crimson",
              "sienna",
              '#7f7f7f']


# --- GENIE interaction mode breakdown ---
genie_mode_list = [0, 10, 1, 2, 3]
genie_mode_labels = [r'$\nu_{\mu}$ CC QE', r'$\nu_{\mu}$ CC MEC', r'$\nu_{\mu}$ CC RES', r'$\nu_{\mu}$ CC SIS/DIS', r'$\nu_{\mu}$ CC COH'] #, r'not $\nu$']
# [r"non FV $\nu$", "Cosmic", r"$\nu$ NC"] + 
genie_mode_colors = ["#9b5580", "#390C1E", "#2c7c94", "#D88A3B", "#BFB17C"] #, '#7f7f7f']
# ["sienna", "crimson", "darkgreen"] + 