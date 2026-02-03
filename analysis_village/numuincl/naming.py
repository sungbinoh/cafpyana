PRELIM_LABEL = 'SBND Preliminary Simulation'
INTERNAL_LABEL = 'SBND Internal Simulation'

#PANDORA_QUALIFIER = r'Pandora TPC reconstruction ($\chi^2$ cut, by TPC)'
#PANDORA_QUALIFIER = r'Algorithm-based TPC reconstruction, local $\chi^2$ optical cut'
PANDORA_QUALIFIER = r'Flash Match ($\chi^2$ cut, by TPC) + hit-based TPC reconstruction'
PANDORA_QUALIFIER = f'Pandora '
PANDORA_QUALIFIER_PRELIM_LABEL = PANDORA_QUALIFIER + PRELIM_LABEL
PANDORA_QUALIFIER_INTERNAL_LABEL = PANDORA_QUALIFIER + INTERNAL_LABEL
PANDORA_QUALIFIER_INTERNAL_LABEL = PANDORA_QUALIFIER + INTERNAL_LABEL

#SPINE_QUALIFIER = r'SPINE TPC reconstruction (best match, whole detector)'
SPINE_QUALIFIER = 'Flash Match (best match, whole detector) + SPINE (ML) TPC reconstruction'
SPINE_QUALIFIER = f'SPINE '
SPINE_QUALIFIER_PRELIM_LABEL = SPINE_QUALIFIER + PRELIM_LABEL
SPINE_QUALIFIER_INTERNAL_LABEL = SPINE_QUALIFIER + INTERNAL_LABEL
SPINE_QUALIFIER_INTERNAL_LABEL = SPINE_QUALIFIER + INTERNAL_LABEL
# Cut lists
PAND_CUTS = ['flashpe','flashmatch','fv','muon','cosmic','lowz']
PAND_CUT_LABELS = ['Flash PE > 2000','Flash Match\n(best BCFM)','Fiducial','Has Muon','Flash Score\n(score cut)','Low Z']

PAND_CUTS_CONT = ['flashpe','flashmatch','fv','muon','cont','cosmic']
PAND_CUT_LABELS_CONT = ['Flash PE > 2000','Flash Match (by TPC)','Fiducial','Has Muon','Muon Contained','Flash Score']

SPINE_CUTS = ['time_contained','cosmic','fv','muon','cosmic_score','lowz','start_dedx']
SPINE_CUT_LABELS = ['Time Contained','Flash Match\n(best match, whole detector)','Fiducial','Has Muon','Flash Score','Low Z','Start dE/dx']

SPINE_CUTS_CONT = ['time_contained','cosmic','fv','muon','cont']
SPINE_CUT_LABELS_CONT = ['Time Contained','Flash Match\n(best match, whole detector)','Fiducial','Has Muon','Muon Contained']

SPINE_CUTS_CONT_NOFM = ['time_contained','fv','muon','cont']
SPINE_CUT_LABELS_CONT_NOFM = ['Time Contained','Fiducial','Has Muon','Contained']