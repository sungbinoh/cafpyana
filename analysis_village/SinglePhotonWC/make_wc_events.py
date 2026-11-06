"""One row per CAF event, index entry: Wire-Cell bookkeeping, so events
without a candidate (no row in the wc table) can be accounted for.

  valid        1: Wire-Cell output exists for the event. 0: none (a step-2
               job that failed every try, a partial merge, or an event lost
               upstream); a production failure, not a reconstruction outcome
  n_nu         candidates (T_tagger rows)
  nu_n_flashes, nu_n_flashes_in_window, nu_n_bundles   from Trun
  n_bundle_reason_<r>   T_bundle rows (flash-cluster bundles) per reason:
               0 candidate (main), 1 candidate (demoted main), 2 all
               cosmic-vetoed, 3 below the 15 cm length floor, 4 STM-only,
               5 nothing eligible, 6 flash-group dedup

valid = 1 with n_nu = 0 is a processed event Wire-Cell found no candidate in;
the reason counts say which step dropped its bundles (none: no bundle).
Joins to mcnu on (file, __ntuple, entry).
"""

import pandas as pd

from makedf import wirecell as wc

BUNDLE_REASONS = range(7)
RUN_COLUMNS = ["nu_n_flashes", "nu_n_flashes_in_window", "nu_n_bundles"]


def make_wc_event_df(f):
    ev = wc.make_wc_evtdf(f)[["run", "subrun", "evt", "valid", "n_nu"]]
    run = wc.make_wc_rundf(f).reset_index("row", drop=True)[RUN_COLUMNS]
    reason = wc.make_wc_bundledf(f, ["reason"]).reason
    n = (pd.crosstab(reason.index.get_level_values("entry"), reason)
         .reindex(index=ev.index, columns=BUNDLE_REASONS, fill_value=0)
         .add_prefix("n_bundle_reason_"))
    n.columns.name = None
    # valid = 0 events have no Trun row: their counts are 0
    return ev.join(run.reindex(ev.index).fillna(0).astype(int)).join(n)
