import pandas as pd
import numpy as np
import sqlite3
import uproot

def _etau_iov(run):
    return __iov(run, SBND_etau_cal_iovdf)

def __iov(run, df):
    return pd.cut(run, list(df.run) + [np.inf], labels=df.iov, right=False)


def get_etau(hit_df):
    iov = _etau_iov(hit_df.run)
    etaudf = pd.DataFrame({"iov": iov})
    etau = etaudf.merge(SBND_etau_cal_df, on=["iov"], how="left", validate="many_to_one")
    etau.index = hit_df.index
    etau_out = np.where(hit_df.tpc == 0, etau.etau_E, etau.etau_W)
    etau_out = pd.Series(etau_out).fillna(np.inf).to_numpy()
    return etau_out

##############################
# SBND calo db files
##############################
sbnd_data_v = "v01_42_00"
SBND_etau_cal_f = "/cvmfs/sbnd.opensciencegrid.org/products/sbnd/sbnd_data/" + sbnd_data_v + "/CalibrationDatabase/tpc_elifetime.db"
SBND_etau_cal_db = "tpc_elifetime_data"
SBND_etau_cal_iov = "tpc_elifetime_iovs"

conn = sqlite3.connect(SBND_etau_cal_f)
cursor = conn.cursor()
cursor.execute("SELECT * FROM %s" % SBND_etau_cal_db)
rows = cursor.fetchall()
data = list(zip(*rows))
SBND_etau_cal_df = pd.DataFrame({
  "iov": data[0],
  "etau_E": data[5],
  "etau_W": data[8],  
})

cursor.execute("SELECT * FROM %s WHERE ACTIVE=1" % SBND_etau_cal_iov)
rows = cursor.fetchall()
data = list(zip(*rows))
SBND_etau_cal_iovdf = pd.DataFrame({
  "iov": data[0],
  "begin_time": data[1],
})
SBND_etau_cal_iovdf["run"] = SBND_etau_cal_iovdf.begin_time % 1000000000
SBND_etau_cal_iovdf.sort_values(by="run", inplace=True)
conn.close()
