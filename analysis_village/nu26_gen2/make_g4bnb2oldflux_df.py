from makedf.makedf import *
from pyanalib.pandas_helpers import *

INPUT_G4BNB = "/exp/sbnd/app/users/sungbino/framework_dev/g4bnb2old/nuMom_immParent_G4BNB.root"
INPUT_OLD = "/exp/sbnd/app/users/sungbino/framework_dev/g4bnb2old/nuMom_immParent_oldFiles.root"

def make_g4bnb2oldflux_df(f):
    this_mcdf = loadbranches(f["recTree"], branches)

    print(this_mcdf)
    
