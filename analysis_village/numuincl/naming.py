PRELIM_LABEL = 'SBND Preliminary Simulation'
INTERNAL_LABEL = 'SBND Internal Simulation'

#PANDORA_QUALIFIER = r'Pandora TPC reconstruction ($\chi^2$ cut, by TPC)'
#PANDORA_QUALIFIER = r'Algorithm-based TPC reconstruction, local $\chi^2$ optical cut'
PANDORA_QUALIFIER = r'Flash Match ($\chi^2$ cut, by TPC) + hit-based TPC reconstruction'
PANDORA_QUALIFIER_PRELIM_LABEL = PANDORA_QUALIFIER + '\n' + PRELIM_LABEL
PANDORA_QUALIFIER_INTERNAL_LABEL = PANDORA_QUALIFIER + '\n' + INTERNAL_LABEL

#SPINE_QUALIFIER = r'SPINE TPC reconstruction (best match, whole detector)'
SPINE_QUALIFIER = 'Flash Match (best match, whole detector) + SPINE (ML) TPC reconstruction'
SPINE_QUALIFIER_PRELIM_LABEL = SPINE_QUALIFIER + '\n' + PRELIM_LABEL
SPINE_QUALIFIER_INTERNAL_LABEL = SPINE_QUALIFIER + '\n' + INTERNAL_LABEL

# Cut lists
PAND_CUTS = ['cosmic','fv','muon']
PAND_CUT_LABELS = ['Flash Match\n($\chi^2$ cut, by TPC)','Fiducial','Has Muon']

PAND_CUTS_CONT = PAND_CUTS + ['cont']
PAND_CUT_LABELS_CONT = PAND_CUT_LABELS + ['Muon Contained']

SPINE_CUTS = ['cosmic','fv','muon','start_dedx','cosmic_score']
SPINE_CUT_LABELS = ['Flash Match\n(best match, whole detector)','Fiducial','Has Muon','Start dE/dx','Flash Score']

SPINE_CUTS_CONT = ['cosmic','fv','muon','cont']
SPINE_CUT_LABELS_CONT = ['Flash Match\n(best match, whole detector)','Fiducial','Has Muon','Muon Contained']

SPINE_CUTS_CONT_NOFM = ['fv','muon','cont']
SPINE_CUT_LABELS_CONT_NOFM = ['Fiducial','Has Muon','Contained']