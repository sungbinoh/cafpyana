from time import time
t0 = time()
import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys,os
import argparse
from tqdm import tqdm

#CAFPYANA working directory
CAFPYANA_WD = '/exp/sbnd/app/users/brindenc/develop/cafpyana'
os.environ['CAFPYANA_WD'] = CAFPYANA_WD

cafpyana_wd = os.environ.get('CAFPYANA_WD')
if cafpyana_wd and cafpyana_wd not in sys.path:
    sys.path.insert(0, cafpyana_wd)
    sys.path.insert(0, cafpyana_wd + '/pyanalib')

#My imports 
SBNDANA_DIR = '/exp/sbnd/app/users/brindenc/develop/cafpyana/analysis_village/numuincl/sbnd'
sys.path.append(SBNDANA_DIR)
plt.style.use(f'{SBNDANA_DIR}/plotlibrary/numu2025.mplstyle')

from sbnd.cafclasses.slice import CAFSlice
from sbnd.cafclasses.interaction import CAFInteraction
from sbnd.cafclasses.nu import NU
from sbnd.constants import *
from sbnd.numu.numu_constants import *
from sbnd.detector.definitions import * # >= v2.2
t1 = time()
print(f'-- Time taken to import: {t1-t0:.2f}s')

DATA_DIR = '/exp/sbnd/data/users/brindenc/analyze_sbnd/numu/v10_06_00_validation/both'
PLOT_DIR = f'{DATA_DIR}/plots'
#parser = argparse.ArgumentParser()
# parser.add_argument('--fname','-f', type=str, required=True, help='name of the file to process (directory not included, its in the DATA_DIR in the script)')
# args = parser.parse_args()
# FNAME = args.fname
MC_FNAME = 'mc_full.df'
OFFBEAM_FNAME = 'data_offbeam.df'
DATA_FNAME = 'data_dev.df'
# Use postprocessed files if they exist
if os.path.exists(f'{DATA_DIR}/{MC_FNAME.replace(".df","_postprocess.df")}'):
  MC_FNAME = f'{MC_FNAME.replace(".df","_postprocess.df")}'
if os.path.exists(f'{DATA_DIR}/{OFFBEAM_FNAME.replace(".df","_postprocess.df")}'):
  OFFBEAM_FNAME = f'{OFFBEAM_FNAME.replace(".df","_postprocess.df")}'
if os.path.exists(f'{DATA_DIR}/{DATA_FNAME.replace(".df","_postprocess.df")}'):
  DATA_FNAME = f'{DATA_FNAME.replace(".df","_postprocess.df")}'
APPLY_CUTS = False
LABEL = 'SBND Internal'
PROCESS_PANDORA = True
PROCESS_SPINE = True
PROCESS_MCNU = True

print(f'MC file: {DATA_DIR}/{MC_FNAME}')
print(f'Offbeam file: {DATA_DIR}/{OFFBEAM_FNAME}')
print(f'Data file: {DATA_DIR}/{DATA_FNAME}')

import h5py
mc_pand_keys = []
offbeam_pand_keys = []
data_pand_keys = []
mc_inter_keys = []
offbeam_inter_keys = []
data_inter_keys = []
mcnu_keys = []
hdr_keys = []
offbeam_hdr_keys = []
data_hdr_keys = []
# Check the keys in the H5 Files
for fname in [MC_FNAME,OFFBEAM_FNAME,DATA_FNAME]:
  print(f'File: {DATA_DIR}/{fname}')
  with h5py.File(f'{DATA_DIR}/{fname}', 'r') as f:
    # Store all keys at the root level
    for key in f.keys():
      if 'pand' in key:
        if fname == MC_FNAME:
          mc_pand_keys.append(key)
        elif fname == OFFBEAM_FNAME:
          offbeam_pand_keys.append(key)
        elif fname == DATA_FNAME:
          data_pand_keys.append(key)
      elif 'evt' in key and 'pand' not in key:
        if fname == MC_FNAME:
          mc_inter_keys.append(key)
        elif fname == OFFBEAM_FNAME:
          offbeam_inter_keys.append(key)
        elif fname == DATA_FNAME:
          data_inter_keys.append(key)
      elif 'mcnu' in key:
        if fname == MC_FNAME:
          mcnu_keys.append(key)
      elif 'hdr' in key:
        if fname == MC_FNAME:
          hdr_keys.append(key)
        elif fname == OFFBEAM_FNAME:
          offbeam_hdr_keys.append(key)
        elif fname == DATA_FNAME:
          data_hdr_keys.append(key)

LIVETIME_DATA = 9.51e5 # From medulla
POT_MC = 0
for key in hdr_keys:
  POT_MC += pd.read_hdf(f'{DATA_DIR}/{MC_FNAME}',key=key).pot.sum()
LIVETIME_DATAOFFBEAM = 0
for key in offbeam_hdr_keys:
  LIVETIME_DATAOFFBEAM += pd.read_hdf(f'{DATA_DIR}/{OFFBEAM_FNAME}',key=key).noffbeambnb.sum()
POT_DATA = 0
#LIVETIME_DATA = 0 #FIXME: How do i find this?
print('ENSURE THAT LIVETIME_DATA IS CORRECT!!')
for key in data_hdr_keys:
  _hdr_df = pd.read_hdf(f'{DATA_DIR}/{DATA_FNAME}',key=key)
  POT_DATA += _hdr_df.pot.sum()
  #LIVETIME_DATA += _hdr_df.noffbeambnb.sum()

print(f'data POT: {POT_DATA:.2e}, mc POT: {POT_MC:.2e}, offbeam livetime: {LIVETIME_DATAOFFBEAM:.2e}, data livetime: {LIVETIME_DATA:.2e}')
t2 = time()
print(f'-- Time taken to get POT and keys: {t2-t1:.2f}s')
#MCnu
if PROCESS_MCNU:
  print('*'*35+' MCnu '+'*'*35)
  df_lengths = []
  for i,key in tqdm(enumerate(mcnu_keys),total=len(mcnu_keys),desc='Loading MCnu'):
    if i == 0:
      mcnu = NU.load(f'{DATA_DIR}/{MC_FNAME}',key=key)
      df_lengths.append(len(mcnu.data))
    else:
      _mcnu = NU.load(f'{DATA_DIR}/{MC_FNAME}',key=key)
      df_lengths.append(len(_mcnu.data))
      mcnu.combine(_mcnu)
  print(np.array(df_lengths).sum(),len(mcnu.data))

  # Post process
  mcnu.scale_to_pot(POT_DATA,sample_pot=POT_MC)

  mcnu.add_fv()
  mcnu.add_av()

  mcnu.cut_muon(cut=APPLY_CUTS,min_ke=0.1)
  mcnu.cut_fv(cut=APPLY_CUTS)
  mcnu.cut_cosmic(cut=APPLY_CUTS)
  mcnu.cut_cont(cut=APPLY_CUTS)

  with open('mcnu_keys.txt','w') as f:
    for k in mcnu.data.keys():
      f.write(f'{k}\n')
t3 = time()
print(f'-- Time taken to post process MCnu: {t3-t2:.2f}s')
#Pandora
if PROCESS_PANDORA:
  print('*'*35+' Pandora '+'*'*35)
  #TODO: Process the Pandora data
  #slc_data = None 
  df_lengths = []
  for i,key in tqdm(enumerate(mc_pand_keys),total=len(mc_pand_keys),desc='Loading Pandora'):
    if i == 0:
      slc = CAFSlice.load(f'{DATA_DIR}/{MC_FNAME}', key=key, pot=POT_MC)
      slc.scale_to_pot(POT_DATA,sample_pot=POT_MC)
      df_lengths.append(len(slc.data))
    else:
      _slc = CAFSlice.load(f'{DATA_DIR}/{MC_FNAME}', key=key, pot=POT_MC)
      _slc.scale_to_pot(POT_DATA,sample_pot=POT_MC)
      df_lengths.append(len(_slc.data))
      slc.combine(_slc)
  print(np.array(df_lengths).sum(),len(slc.data))
  #Add offbeam data
  offbeam_df_lengths = []
  for i,key in tqdm(enumerate(offbeam_pand_keys),total=len(offbeam_pand_keys),desc='Loading Pandora offbeam'):
    _slc_intime = CAFSlice.load(f'{DATA_DIR}/{OFFBEAM_FNAME}', key=key, livetime=LIVETIME_DATAOFFBEAM)
    _slc_intime.scale_to_livetime(LIVETIME_DATA,sample_livetime=LIVETIME_DATAOFFBEAM)
    df_lengths.append(len(_slc_intime.data))
    slc.combine(_slc_intime,duplicate_ok=True)
  print(np.array(offbeam_df_lengths+df_lengths).sum(),len(slc.data))
  #Make slc_data
  data_df_lengths = []
  for i,key in tqdm(enumerate(data_pand_keys),total=len(data_pand_keys),desc='Loading Pandora data'):
    if i == 0:
      slc_data = CAFSlice.load(f'{DATA_DIR}/{DATA_FNAME}', key=key)
      data_df_lengths.append(len(slc_data.data))
    else:
      _slc_data = CAFSlice.load(f'{DATA_DIR}/{DATA_FNAME}', key=key)
      data_df_lengths.append(len(_slc_data.data))
      slc_data.combine(_slc_data)
  print(np.array(data_df_lengths).sum(),len(slc_data.data))

  #Post process
  #Set pandora containment and event type
  mcnu = slc.set_mcnu_containment(mcnu)
  mcnu.add_event_type('pandora')

  if slc_data is not None:
    slc_list = [slc,slc_data]
  else:
    slc_list = [slc]
  for s in slc_list:
    s.clean(dummy_vals=[-9999,-999,999,9999,-5])

    s.add_has_muon()
    s.add_in_av()
    s.add_in_fv()
    s.cut_cosmic(cut=APPLY_CUTS,fmatch_score=320,nu_score=0.5,use_opt0=True,use_isclearcosmic=True)
    s.cut_fv(cut=APPLY_CUTS)
    s.cut_muon(cut=APPLY_CUTS,min_ke=0.1)
    s.cut_is_cont(cut=APPLY_CUTS)

  #I know there is an offbeam sample combined. There is a nan check in the function that assigns these as cosmics
  slc.add_event_type()

  with open('slc_keys.txt','w') as f:
    for k in slc.data.keys():
      f.write(f'{k}\n')

  PAND_CUTS = ['cosmic','fv','muon']
  PAND_CUTS_CONT = PAND_CUTS + ['cont']
  pur,eff,f1 = slc.get_pur_eff_f1(mcnu,PAND_CUTS_CONT,categories=[0,1])
  print('Pandora cuts:')
  print(PAND_CUTS_CONT)
  print('Pandora pur, eff, f1:')
  print(pur,eff,f1)
else:
  slc_data = None
  slc = None
t4 = time()
print(f'-- Time taken to post process Pandora: {t4-t3:.2f}s')
#SPINE
if PROCESS_SPINE:
  print('*'*35+' SPINE '+'*'*35)
  #TODO: Process the SPINE data
  #inter_data = None
  #Add MC information
  df_lengths = []
  for i,key in tqdm(enumerate(mc_inter_keys),total=len(mc_inter_keys),desc='Loading SPINE'):
    if i == 0:
      inter = CAFInteraction.load(f'{DATA_DIR}/{MC_FNAME}', key=key, pot=POT_MC)
      inter.scale_to_pot(POT_DATA,sample_pot=POT_MC)
      df_lengths.append(len(inter.data))
    else:
      _inter = CAFInteraction.load(f'{DATA_DIR}/{MC_FNAME}', key=key, pot=POT_MC)
      _inter.scale_to_pot(POT_DATA,sample_pot=POT_MC)
      df_lengths.append(len(_inter.data))
      inter.combine(_inter)
  print(np.array(df_lengths).sum(),len(inter.data))

  #Add offbeam data
  offbeam_df_lengths = []
  for i,key in tqdm(enumerate(offbeam_inter_keys),total=len(offbeam_inter_keys),desc='Loading SPINE'):
    _inter = CAFInteraction.load(f'{DATA_DIR}/{OFFBEAM_FNAME}', key=key, livetime=LIVETIME_DATAOFFBEAM)
    _inter.scale_to_livetime(LIVETIME_DATA,sample_livetime=LIVETIME_DATAOFFBEAM)
    offbeam_df_lengths.append(len(_inter.data))
    inter.combine(_inter,duplicate_ok=True)
  print(np.array(offbeam_df_lengths+df_lengths).sum(),len(inter.data))

  #Make inter_data
  data_df_lengths = []
  for i,key in tqdm(enumerate(data_inter_keys),total=len(data_inter_keys),desc='Loading SPINE data'):
    if i == 0:
      inter_data = CAFInteraction.load(f'{DATA_DIR}/{DATA_FNAME}', key=key)
      data_df_lengths.append(len(inter_data.data))
    else:
      _inter_data = CAFInteraction.load(f'{DATA_DIR}/{DATA_FNAME}', key=key)
      data_df_lengths.append(len(_inter_data.data))
      inter_data.combine(_inter_data)
  print(np.array(data_df_lengths).sum(),len(inter_data.data))
  print(len(inter_data.data),np.array(data_df_lengths).sum(),len(inter.data),np.array(offbeam_df_lengths+df_lengths).sum())

  #Post process
  #Set pandora containment and event type
  mcnu = inter.set_mcnu_containment(mcnu)
  mcnu.add_event_type('spine')

  if inter_data is not None:
    inter_list = [inter,inter_data]
  else:
    inter_list = [inter]
  for iinter in inter_list:
    iinter.clean(dummy_vals=[-9999,-999,999,9999,-5,np.inf,-np.inf])
    iinter.add_in_fv()
    iinter.add_in_av()

    iinter.cut_cosmic(cut=APPLY_CUTS)
    iinter.cut_fv(cut=APPLY_CUTS)
    iinter.cut_muon(cut=APPLY_CUTS,min_ke=0.1)
    iinter.cut_is_cont(cut=APPLY_CUTS)
    iinter.cut_start_dedx(cut=APPLY_CUTS,dedx=4.17)
    iinter.cut_cosmic_score(cut=APPLY_CUTS,score=102.35)

  inter.add_event_type()

  with open('inter_keys.txt','w') as f:
    for k in inter.data.keys():
      f.write(f'{k}\n')

  SPINE_CUTS = ['cosmic','fv','muon','start_dedx','cosmic_score']
  SPINE_CUTS_CONT = ['cosmic','fv','muon','cont']
  SPINE_CUTS_CONT_NOFM = ['fv','muon','cont']
  pur,eff,f1 = inter.get_pur_eff_f1(mcnu,SPINE_CUTS,categories=[0,1])
  print('SPINE cuts:')
  print(SPINE_CUTS)
  print('SPINE pur, eff, f1:')
  print(pur,eff,f1)
  pur,eff,f1 = inter.get_pur_eff_f1(mcnu,SPINE_CUTS_CONT,categories=[0])
  print('SPINE contained cuts:')
  print(SPINE_CUTS_CONT)
  print('SPINE contained pur, eff, f1:')
  print(pur,eff,f1)
  pur,eff,f1 = inter.get_pur_eff_f1(mcnu,SPINE_CUTS_CONT_NOFM,categories=[0])
  print('SPINE contained no FM cuts:')
  print(SPINE_CUTS_CONT_NOFM)
  print('SPINE contained no FM pur, eff, f1:')
  print(pur,eff,f1)
else:
  inter_data = None
  inter = None
t5 = time()
print(f'-- Time taken to post process SPINE: {t5-t4:.2f}s')

#Save data to h5 files
if PROCESS_MCNU:
  mcnu.data.to_hdf(f"{DATA_DIR}/{MC_FNAME.replace('.df','_postprocess.df')}",key='mcnu')
if PROCESS_PANDORA:
  slc.data.to_hdf(f"{DATA_DIR}/{MC_FNAME.replace('.df','_postprocess.df')}",key='pandora')
if PROCESS_SPINE:
  inter.data.to_hdf(f"{DATA_DIR}/{MC_FNAME.replace('.df','_postprocess.df')}",key='spine')
if slc_data is not None:
  slc_data.data.to_hdf(f"{DATA_DIR}/{DATA_FNAME.replace('.df','_postprocess.df')}",key='pandora')
if inter_data is not None:
  inter_data.data.to_hdf(f"{DATA_DIR}/{DATA_FNAME.replace('.df','_postprocess.df')}",key='spine')
t6 = time()
print(f'-- Time taken to save data to h5 files: {t6-t5:.2f}s')
print(f'-- Total time taken: {t6-t0:.2f}s')
