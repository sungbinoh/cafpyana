from pyanalib.pandas_helpers import loadbranches
from analysis_village.SinglePhotonWC.make_wc_candidates import make_wc_candidate_df
from analysis_village.SinglePhotonWC.make_wc_events import make_wc_event_df


def make_hdrdf(f):
    """rec.hdr pot and noffbeambnb (as makedf.make_hdrdf, without its CVMFS imports)."""
    return loadbranches(f["recTree"], ["rec.hdr.run", "rec.hdr.subrun", "rec.hdr.evt", "rec.hdr.pot",
                                       "rec.hdr.noffbeambnb", "rec.hdr.nbnbinfo", "rec.hdr.first_in_subrun"]).rec.hdr


# Data (on- or off-beam): the sp_mc.py tables without the truth table, plus the
# header (pot, noffbeambnb) for the POT / off-beam gate normalisation.
DFS = [make_wc_candidate_df, make_wc_event_df, make_hdrdf]
NAMES = ["wc", "wcevt", "hdr"]
