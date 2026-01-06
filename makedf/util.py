import numpy as np
import sys

def mag(x, y, z):
    return np.sqrt(x**2 + y**2 + z**2)

def magdf(df):
    return mag(df.x, df.y, df.z)

def mag2d(x, y):
    return np.sqrt(x**2 + y**2)

def dmagdf(df1, df2):
    return mag(df1.x - df2.x, df1.y - df2.y, df1.z - df2.z)

def dotdf(df1, df2):
    return df1.x*df2.x + df1.y*df2.y + df1.z*df2.z 

def unitdf(df):
    return df.divide(magdf(df), axis=0)

def InAV(df,det="SBND"):
    if det == "SBND":
        xmin = -195
        ymin = -195
        zmin = 5
        xmax = 195
        ymax = 195
        zmax = 495
    return (df.x > xmin) & (df.x < xmax) & (df.y > ymin) & (df.y < ymax) & (df.z > zmin) & (df.z < zmax)

def InFV(df, inzback, inx=10, iny=10, inzfront=10, det="ICARUS"):
    if det == "ICARUS":
        xmin_C0 = -358.49
        xmax_C0 = -61.94

        xmax_C1 = -xmin_C0
        xmin_C1 = -xmax_C0

        ymin = -181.85999999999999
        ymax = 134.96

        zmin = -894.950652270838
        zmax = 894.950652270838

        xmin_C0 = xmin_C0 + inx
        xmax_C0 = xmax_C0 - inx
        xmin_C1 = xmin_C1 + inx
        xmax_C1 = xmax_C1 - inx

        ymin = ymin + iny
        ymax = ymax - iny

        zmin = zmin + inzfront
        zmax = zmax - inzback

        return (((df.x < xmax_C0) & (df.x > xmin_C0)) | ((df.x < xmax_C1) & (df.x > xmin_C1))) &\
            (df.y < ymax) & (df.y > ymin) & (df.z < zmax) & (df.z > zmin)
    
    elif det == "SBND":
        xmin = -190
        ymin = -190
        zmin = 10
        xmax = 190
        ymax =  190
        zmax =  450.
        return (df.x > xmin) & (df.x < xmax) & (df.y > ymin) & (df.y < ymax) & (df.z > zmin) & (df.z < zmax)

    elif det == "SBND_nohighyz":
        xmax = 190.
        zmin = 10.
        zmax = 450.
        ymax_highz = 100.
        pass_xz = (np.abs(df.x) < xmax) & (df.z > zmin) & (df.z < zmax)
        pass_y = ((df.z < 250) & (np.abs(df.y) < 190.)) | ((df.z > 250) & (df.y > -190.) & (df.y < ymax_highz))
        return pass_xz & pass_y

    else:
        raise NameError("DETECTOR not valid, should be SBND or ICARUS")

def TrkInFV(df):
    return InFV(df, 15.)

def SlcInFV(df):
    return InFV(df, 100.)


