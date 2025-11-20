from time import time
t0 = time()
import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys,os
import argparse
import gc
from tqdm import tqdm
from naming import * #Import cut names

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

DATA_DIR = '/pnfs/sbnd/persistent/users/brindenc/numuincl_analysis/v10_06_00/pandora'
#parser = argparse.ArgumentParser()
# parser.add_argument('--fname','-f', type=str, required=True, help='name of the file to process (directory not included, its in the DATA_DIR in the script)')
# args = parser.parse_args()
# FNAME = args.fname
#MC_FNAME = 'mc_med_syst2.df'
#MC_FNAME = 'mc_med2.df'
#MC_FNAME = 'pds_var/cv/mc_small_2.df'
#MC_FNAME = 'pds_var/PMTGainFluct/mc_small_2.df'
#MC_FNAME = 'pds_var/PMTHighNoise/mc_small_2.df'
#MC_FNAME = 'pds_var/PMTLowEff/mc_small_2.df'
MC_FNAME = 'mc_syst/mc_tinypand_fullsyst_cut3.df'
OFFBEAM_FNAME = 'data_offbeam2.df'
DATA_FNAME = 'data_dev2.df'
APPLY_CUTS = False
PROCESS_PANDORA = True
PROCESS_SPINE = False
PROCESS_MCNU = True
PROCESS_WEIGHTS = False
MAX_OFFBEAM_FILES = 2  # Limit number of offbeam files to load (set to None to load all)
MAX_MC_FILES = None # Limit number of MC files to load (set to None to load all)
if APPLY_CUTS:
  suffix = 'postprocess_cut'
else:
  suffix = 'postprocess'

if PROCESS_PANDORA and not PROCESS_SPINE:
  suffix += '_pandora'
elif PROCESS_SPINE and not PROCESS_PANDORA:
  suffix += '_spine'

assert os.path.exists(f'{DATA_DIR}/{MC_FNAME}'), f'MC file {DATA_DIR}/{MC_FNAME} does not exist'
assert os.path.exists(f'{DATA_DIR}/{OFFBEAM_FNAME}'), f'Offbeam file {DATA_DIR}/{OFFBEAM_FNAME} does not exist'
assert os.path.exists(f'{DATA_DIR}/{DATA_FNAME}'), f'Data file {DATA_DIR}/{DATA_FNAME} does not exist'

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
      print(f'- Key: {key}')
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

# Limit offbeam files if specified
if MAX_OFFBEAM_FILES is not None:
  if len(offbeam_pand_keys) > MAX_OFFBEAM_FILES:
    print(f'Limiting offbeam Pandora files from {len(offbeam_pand_keys)} to {MAX_OFFBEAM_FILES}')
    offbeam_pand_keys = offbeam_pand_keys[:MAX_OFFBEAM_FILES]
  if len(offbeam_inter_keys) > MAX_OFFBEAM_FILES:
    print(f'Limiting offbeam SPINE files from {len(offbeam_inter_keys)} to {MAX_OFFBEAM_FILES}')
    offbeam_inter_keys = offbeam_inter_keys[:MAX_OFFBEAM_FILES]
  if len(offbeam_hdr_keys) > MAX_OFFBEAM_FILES:
    print(f'Limiting offbeam header files from {len(offbeam_hdr_keys)} to {MAX_OFFBEAM_FILES}')
    offbeam_hdr_keys = offbeam_hdr_keys[:MAX_OFFBEAM_FILES]

if MAX_MC_FILES is not None:
  if len(mc_pand_keys) > MAX_MC_FILES:
    print(f'Limiting MC Pandora files from {len(mc_pand_keys)} to {MAX_MC_FILES}')
    mc_pand_keys = mc_pand_keys[:MAX_MC_FILES]
  if len(mc_inter_keys) > MAX_MC_FILES:
    print(f'Limiting MC SPINE files from {len(mc_inter_keys)} to {MAX_MC_FILES}')
    mc_inter_keys = mc_inter_keys[:MAX_MC_FILES]
  if len(mcnu_keys) > MAX_MC_FILES:
    print(f'Limiting MC MCnu files from {len(mcnu_keys)} to {MAX_MC_FILES}')
    mcnu_keys = mcnu_keys[:MAX_MC_FILES]
  if len(hdr_keys) > MAX_MC_FILES:
    print(f'Limiting MC header files from {len(hdr_keys)} to {MAX_MC_FILES}')
    hdr_keys = hdr_keys[:MAX_MC_FILES]

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
if len(data_hdr_keys) == 0:
  raise ValueError('No data header keys found')
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
      del _mcnu  # Free memory after combining
      gc.collect()
  print(np.array(df_lengths).sum(),len(mcnu.data))

  # Post process
  mcnu.scale_to_pot(POT_DATA,sample_pot=POT_MC)

  mcnu.add_fv()
  mcnu.add_av()

  #Never cut on truth quantities
  mcnu.cut_muon(cut=False,min_ke=0.1)
  mcnu.cut_fv(cut=False)
  mcnu.cut_cosmic(cut=False)
  mcnu.cut_cont(cut=False)

  with open('mcnu_keys.txt','w') as f:
    for k in mcnu.data.keys():
      f.write(f'{k}\n')
t3 = time()
print(f'-- Time taken to post process MCnu: {t3-t2:.2f}s')
#Pandora
if PROCESS_PANDORA:
  print('*'*35+' Pandora '+'*'*35)
  #TODO: Process the Pandora data
  slc_data = None 
  df_lengths = []
  # Load all MC files without scaling first (memory efficient)
  for i,key in tqdm(enumerate(mc_pand_keys),total=len(mc_pand_keys),desc='Loading Pandora'):
    if i == 0:
      slc = CAFSlice.load(f'{DATA_DIR}/{MC_FNAME}', key=key, pot=POT_MC)
      df_lengths.append(len(slc.data))
    else:
      _slc = CAFSlice.load(f'{DATA_DIR}/{MC_FNAME}', key=key, pot=POT_MC)
      df_lengths.append(len(_slc.data))
      slc.combine(_slc)
      del _slc  # Free memory after combining
      gc.collect()  # Force garbage collection
  # Scale MC once after combining (memory efficient)
  slc.scale_to_pot(POT_DATA,sample_pot=POT_MC)
  gc.collect()  # Clean up after scaling
  #print(np.array(df_lengths).sum(),len(slc.data))
  #Add offbeam data - process in batches to avoid memory issues
  # Process offbeam files in smaller batches, combining incrementally
  BATCH_SIZE = 3  # Process 3 files at a time
  slc_offbeam = None
  offbeam_df_lengths = []
  n_batches = (len(offbeam_pand_keys) + BATCH_SIZE - 1) // BATCH_SIZE
  for batch_idx, batch_start in enumerate(range(0, len(offbeam_pand_keys), BATCH_SIZE), 1):
    batch_keys = offbeam_pand_keys[batch_start:batch_start+BATCH_SIZE]
    slc_batch = None
    for i,key in enumerate(batch_keys):
      if i == 0:
        slc_batch = CAFSlice.load(f'{DATA_DIR}/{OFFBEAM_FNAME}', key=key, livetime=LIVETIME_DATAOFFBEAM)
        offbeam_df_lengths.append(len(slc_batch.data))
      else:
        _slc_intime = CAFSlice.load(f'{DATA_DIR}/{OFFBEAM_FNAME}', key=key, livetime=LIVETIME_DATAOFFBEAM)
        offbeam_df_lengths.append(len(_slc_intime.data))
        slc_batch.combine(_slc_intime)
        del _slc_intime
        gc.collect()
    # Combine batch with main offbeam accumulator
    if slc_offbeam is None:
      slc_offbeam = slc_batch
    else:
      slc_offbeam.combine(slc_batch)
      del slc_batch
      gc.collect()
  # Scale offbeam once after combining (memory efficient)
  if slc_offbeam is not None:
    slc_offbeam.scale_to_livetime(LIVETIME_DATA,sample_livetime=LIVETIME_DATAOFFBEAM)
    slc.combine(slc_offbeam,duplicate_ok=True)
    del slc_offbeam  # Free memory
    gc.collect()
  #print(np.array(offbeam_df_lengths+df_lengths).sum(),len(slc.data))
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
      del _slc_data  # Free memory after combining
      gc.collect()
  #print(np.array(data_df_lengths).sum(),len(slc_data.data))
  #Post process
  #Set pandora containment and event type
  print(f'Setting mcnu containment for {len(mcnu.data)} events')
  mcnu = slc.set_mcnu_containment(mcnu)
  # Force garbage collection after memory-intensive operation
  gc.collect()
  print('Adding event type to mcnu')
  mcnu.add_event_type('pandora')
  print('Added event type to mcnu')

  if slc_data is not None:
    slc_list = [slc,slc_data]
  else:
    slc_list = [slc]
  print(f'Post processing {len(slc_list)} slc objects')
  for i,s in enumerate(slc_list):
    s.clean(dummy_vals=[-9999,-999,999,9999,-5])

    s.add_has_muon()
    s.add_in_av()
    s.add_in_fv()
    print(f'Adding event type for {len(s.data)} events')
    #Make slc signal to store the true signal information
    if i == 0:
      #I know there is an offbeam sample combined. There is a nan check in the function that assigns these as cosmics
      s.add_event_type()
      #s.add_stat_unc()
      slc_signal = s.copy()
      #Keep only true signal events
      is_signal = np.isin(s.data.truth.event_type,[0,1])
      slc_signal.data = s.data[is_signal]

    print(f'Adding cuts')
    s.cut_flashmatch(cut=APPLY_CUTS)
    s.cut_fv(cut=APPLY_CUTS)
    s.cut_muon(cut=APPLY_CUTS,min_ke=0.1)
    s.cut_cosmic(cut=APPLY_CUTS,fmatch_score=320,nu_score=0.5,use_opt0=True,use_isclearcosmic=False)
    s.cut_lowz(cut=APPLY_CUTS,z_max=6,include_start=True)
    s.cut_is_cont(cut=False) #Don't apply containment cut


  with open('slc_keys.txt','w') as f:
    for k in slc.data.keys():
      f.write(f'{k}\n')

  #PAND_CUTS = ['cosmic','fv','muon']
  #PAND_CUTS_CONT = PAND_CUTS + ['cont']
  from naming import PAND_CUTS_CONT, PAND_CUTS
  pur,eff,f1,_,_,_ = slc.get_pur_eff_f1(mcnu,PAND_CUTS_CONT,categories=[0,1])
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
  inter_data = None
  #Add MC information - load all without scaling first (memory efficient)
  df_lengths = []
  for i,key in tqdm(enumerate(mc_inter_keys),total=len(mc_inter_keys),desc='Loading SPINE'):
    if i == 0:
      inter = CAFInteraction.load(f'{DATA_DIR}/{MC_FNAME}', key=key, pot=POT_MC)
      df_lengths.append(len(inter.data))
    else:
      _inter = CAFInteraction.load(f'{DATA_DIR}/{MC_FNAME}', key=key, pot=POT_MC)
      df_lengths.append(len(_inter.data))
      inter.combine(_inter)
      del _inter  # Free memory after combining
      gc.collect()
  # Scale MC once after combining (memory efficient)
  inter.scale_to_pot(POT_DATA,sample_pot=POT_MC)
  gc.collect()
  #print(np.array(df_lengths).sum(),len(inter.data))

  #Add offbeam data - process in batches to avoid memory issues
  # Process offbeam files in smaller batches, combining incrementally
  BATCH_SIZE = 3  # Process 3 files at a time
  inter_offbeam = None
  offbeam_df_lengths = []
  n_batches = (len(offbeam_inter_keys) + BATCH_SIZE - 1) // BATCH_SIZE
  for batch_idx, batch_start in enumerate(range(0, len(offbeam_inter_keys), BATCH_SIZE), 1):
    batch_keys = offbeam_inter_keys[batch_start:batch_start+BATCH_SIZE]
    inter_batch = None
    for i,key in enumerate(batch_keys):
      if i == 0:
        inter_batch = CAFInteraction.load(f'{DATA_DIR}/{OFFBEAM_FNAME}', key=key, livetime=LIVETIME_DATAOFFBEAM)
        offbeam_df_lengths.append(len(inter_batch.data))
      else:
        _inter = CAFInteraction.load(f'{DATA_DIR}/{OFFBEAM_FNAME}', key=key, livetime=LIVETIME_DATAOFFBEAM)
        offbeam_df_lengths.append(len(_inter.data))
        inter_batch.combine(_inter)
        del _inter
        gc.collect()
    # Combine batch with main offbeam accumulator
    if inter_offbeam is None:
      inter_offbeam = inter_batch
    else:
      inter_offbeam.combine(inter_batch)
      del inter_batch
      gc.collect()
  # Scale offbeam once after combining (memory efficient)
  if inter_offbeam is not None:
    inter_offbeam.scale_to_livetime(LIVETIME_DATA,sample_livetime=LIVETIME_DATAOFFBEAM)
    inter.combine(inter_offbeam,duplicate_ok=True)
    del inter_offbeam  # Free memory
    gc.collect()
  #print(np.array(offbeam_df_lengths+df_lengths).sum(),len(inter.data))

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
      del _inter_data  # Free memory after combining
      gc.collect()
  #print(np.array(data_df_lengths).sum(),len(inter_data.data))
  #print(len(inter_data.data),np.array(data_df_lengths).sum(),len(inter.data),np.array(offbeam_df_lengths+df_lengths).sum())

  #Post process
  #Set pandora containment and event type
  print(f'Setting mcnu containment for {len(mcnu.data)} events')
  mcnu = inter.set_mcnu_containment(mcnu)
  print('Adding event type to mcnu')
  mcnu.add_event_type('spine')

  if inter_data is not None:
    inter_list = [inter,inter_data]
  else:
    inter_list = [inter]
  for i,iinter in enumerate(inter_list):
    iinter.clean(dummy_vals=[-9999,-999,999,9999,-5,np.inf,-np.inf])
    iinter.add_in_fv()
    iinter.add_in_av()

    if i == 0:
      #I know there is an offbeam sample combined. There is a nan check in the function that assigns these as cosmics
      iinter.add_event_type()
      #iinter.add_stat_unc()
      inter_signal = iinter.copy()
      #Keep only true signal events
      is_signal = np.isin(iinter.data.truth.event_type,[0,1])
      inter_signal.data = iinter.data[is_signal]

    iinter.cut_time_contained(cut=APPLY_CUTS)
    iinter.cut_cosmic(cut=APPLY_CUTS)
    iinter.cut_fv(cut=APPLY_CUTS)
    iinter.cut_muon(cut=APPLY_CUTS,min_ke=0.1)
    iinter.cut_cosmic_score(cut=APPLY_CUTS,score=102.35)
    iinter.cut_lowz(cut=APPLY_CUTS,z_max=6,include_start=True)
    iinter.cut_start_dedx(cut=APPLY_CUTS,dedx=4.17)
    iinter.cut_is_cont(cut=False) #Don't apply containment cut

  with open('inter_keys.txt','w') as f:
    for k in inter.data.keys():
      f.write(f'{k}\n')

  from naming import SPINE_CUTS, SPINE_CUTS_CONT, SPINE_CUTS_CONT_NOFM
  pur,eff,f1,_,_,_ = inter.get_pur_eff_f1(mcnu,SPINE_CUTS,categories=[0,1])
  print('SPINE cuts:')
  print(SPINE_CUTS)
  print('SPINE pur, eff, f1:')
  print(pur,eff,f1)
  pur,eff,f1,_,_,_ = inter.get_pur_eff_f1(mcnu,SPINE_CUTS_CONT,categories=[0])
  print('SPINE contained cuts:')
  print(SPINE_CUTS_CONT)
  print('SPINE contained pur, eff, f1:')
  print(pur,eff,f1)
  # pur,eff,f1,_,_,_ = inter.get_pur_eff_f1(mcnu,SPINE_CUTS_CONT_NOFM,categories=[0])
  # print('SPINE contained no FM cuts:')
  # print(SPINE_CUTS_CONT_NOFM)
  # print('SPINE contained no FM pur, eff, f1:')
  # print(pur,eff,f1)
else:
  inter_data = None
  inter = None
t5 = time()
print(f'-- Time taken to post process SPINE: {t5-t4:.2f}s')

#Add universe weights back before saving
if PROCESS_WEIGHTS:
  print('*'*35+' Adding Universe Weights '+'*'*35)
  if PROCESS_PANDORA:
    print(f'Adding universe weights to Pandora from {len(mc_pand_keys)} keys...')
    slc.add_universe_weights(f'{DATA_DIR}/{MC_FNAME}', keys=mc_pand_keys, duplicate_ok=False)
    gc.collect()
    if 'slc_signal' in locals():
      print(f'Adding universe weights to Pandora signal from {len(mc_pand_keys)} keys...')
      slc_signal.add_universe_weights(f'{DATA_DIR}/{MC_FNAME}', keys=mc_pand_keys, duplicate_ok=False)
      gc.collect()
  if PROCESS_SPINE:
    print(f'Adding universe weights to SPINE from {len(mc_inter_keys)} keys...')
    inter.add_universe_weights(f'{DATA_DIR}/{MC_FNAME}', keys=mc_inter_keys, duplicate_ok=False)
    gc.collect()
    if 'inter_signal' in locals():
      print(f'Adding universe weights to SPINE signal from {len(mc_inter_keys)} keys...')
      inter_signal.add_universe_weights(f'{DATA_DIR}/{MC_FNAME}', keys=mc_inter_keys, duplicate_ok=False)
      gc.collect()
  print('Finished adding universe weights')

#Save data to h5 files
if PROCESS_MCNU:
  mcnu.data.to_hdf(f"{DATA_DIR}/{MC_FNAME.replace('.df',f'_{suffix}.df')}",key='mcnu')
if PROCESS_PANDORA:
  slc.data.to_hdf(f"{DATA_DIR}/{MC_FNAME.replace('.df',f'_{suffix}.df')}",key='pandora')
  if 'slc_signal' in locals():
    slc_signal.data.to_hdf(f"{DATA_DIR}/{MC_FNAME.replace('.df',f'_{suffix}.df')}",key='pandora_signal')
if PROCESS_SPINE:
  inter.data.to_hdf(f"{DATA_DIR}/{MC_FNAME.replace('.df',f'_{suffix}.df')}",key='spine')
  if 'inter_signal' in locals():
    inter_signal.data.to_hdf(f"{DATA_DIR}/{MC_FNAME.replace('.df',f'_{suffix}.df')}",key='spine_signal')
print('Not saving data files, uncomment these lines to do that')
if slc_data is not None:
  slc_data.data.to_hdf(f"{DATA_DIR}/{DATA_FNAME.replace('.df',f'_{suffix}.df')}",key='pandora')
if inter_data is not None:
  inter_data.data.to_hdf(f"{DATA_DIR}/{DATA_FNAME.replace('.df',f'_{suffix}.df')}",key='spine')
t6 = time()
print(f'-- Time taken to save data to h5 files: {t6-t5:.2f}s')
print(f'-- Total time taken: {t6-t0:.2f}s')
