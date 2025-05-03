import numpy as np

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
        xmin = -199.15 + 10
        ymin = -200. + 10
        zmin = 0.0 + 10
        xmax = 199.15 - 10
        ymax =  200. - 10
        zmax =  500. - 50
        return (df.x > xmin) & (df.x < xmax) & (df.y > ymin) & (df.y < ymax) & (df.z > zmin) & (df.z < zmax)
    
    else:
        raise NameError("DETECTOR not valid, should be SBND or ICARUS")

def TrkInFV(df):
    return InFV(df, 15.)

def SlcInFV(df):
    return InFV(df, 100.)


