import os
import sys
import pandas as pd

# Add the head directly to sys.path
workspace_root = os.getcwd()
sys.path.insert(0, workspace_root + "/../../")

from makedf.util import *

PROTON_MASS = 0.938272
NEUTRON_MASS = 0.939565
MUON_MASS = 0.105658
PION_MASS = 0.139570
MASS_A = 22*NEUTRON_MASS + 18*PROTON_MASS - 0.34381
BE = 0.0295
MASS_Ap = MASS_A - NEUTRON_MASS + BE

def neutrino_energy(mu_p, mu_dir, p_p, p_dir):
    mu_E = mag2d(mu_p, MUON_MASS)
    p_E = mag2d(p_p, PROTON_MASS)

    dpT = transverse_kinematics(mu_p, mu_dir, p_p, p_dir)['del_Tp']
    ET = np.sqrt(dpT**2 + MASS_Ap**2) - MASS_Ap
    
    return mu_E + p_E - PROTON_MASS + ET + BE

def transverse_kinematics(mu_p, mu_dir, p_p, p_dir):
    mu_E = mag2d(mu_p, MUON_MASS)
    p_E = mag2d(p_p, PROTON_MASS)

    mu_p_x = mu_p * mu_dir.x
    mu_p_y = mu_p * mu_dir.y
    mu_p_z = mu_p * mu_dir.z
    mu_phi_x = mu_p_x/mag2d(mu_p_x, mu_p_y)
    mu_phi_y = mu_p_y/mag2d(mu_p_x, mu_p_y)

    p_p_x = p_p * p_dir.x
    p_p_y = p_p * p_dir.y
    p_p_z = p_p * p_dir.z
    p_phi_x = p_p_x/mag2d(p_p_x, p_p_y)
    p_phi_y = p_p_y/mag2d(p_p_x, p_p_y)

    mu_Tp_x = mu_phi_y*mu_p_x - mu_phi_x*mu_p_y
    mu_Tp_y = mu_phi_x*mu_p_x - mu_phi_y*mu_p_y
    mu_Tp = mag2d(mu_Tp_x, mu_Tp_y)


    p_Tp_x = mu_phi_y*p_p_x - mu_phi_x*p_p_y
    p_Tp_y = mu_phi_x*p_p_x - mu_phi_y*p_p_y
    p_Tp = mag2d(p_Tp_x, p_Tp_y)



    del_Tp_x = mu_Tp_x + p_Tp_x
    del_Tp_y = mu_Tp_y + p_Tp_y
    del_Tp = mag2d(del_Tp_x, del_Tp_y)


    del_alpha = np.arccos(-(mu_Tp_x*del_Tp_x + mu_Tp_y*del_Tp_y)/(mu_Tp*del_Tp))
    del_phi = np.arccos(-(mu_Tp_x*p_Tp_x + mu_Tp_y*p_Tp_y)/(mu_Tp*p_Tp))

    mu_E = mag2d(mu_p, MUON_MASS)

    p_E = mag2d(p_p, PROTON_MASS)

    R = MASS_A + mu_p_z + p_p_z - mu_E - p_E
    del_Lp = 0.5*R - mag2d(MASS_Ap, del_Tp)**2/(2*R)

    del_p = mag2d(del_Tp, del_Lp)


    return pd.Series({'del_p' : del_p, 
                      'del_Tp' : del_Tp, 
                      'del_phi' : del_phi, 
                      'del_alpha' : del_alpha, 
                      'mu_E' : mu_E, 
                      'p_E' : p_E})
