import numpy as np
import pandas as pd

from analysis_village.mcs.modules.TMatrixDSym import *
from analysis_village.mcs.modules.TrackTrajectory import *
from analysis_village.mcs.modules.RangeE import *

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)s:%(asctime)s:%(name)s:%(message)s')
file_handler = logging.FileHandler('MCSFit.log')
file_handler.setLevel(logging.INFO)
file_handler.setFormatter(formatter)

SEG_LEN = 14
MIN_NSEG = 3
MIN_NHITS_PER_SEG = 5
LAR_RADL_INV = 1./14.0;

# track E calculation
PID = 13
MOM_DEP_CONST = True
E_LOSS_MODE = 0 # Default is MPV Landau. Choose 1 for MIP (constant); 2 for Bethe-Bloche
E_LOSS_NSTEPS = 10

# fit configurations
P_MIN = 0.1
P_MAX = 10.1
P_STEP = 0.01
ANG_RESOL = 5 # [0MOM_DEP_CONSTNG_RESOL 0., 5.5, -3.5, 0.]
CUT_SEG = 0

# get segdf from fitMcsPrepare
def fitMcs(segdf, pid, mom_dep_const, ang_resol=ANG_RESOL, cut_seg=CUT_SEG):
    # TODO: truncate segments near track end to simulate exiting tracks
    if cut_seg > 0:
        if len(segdf) > (cut_seg+1): 
            segdf = segdf[:-cut_seg]
        else:
            logging.info("not enough segments to truncate at track end, skipping traj")
            return None

    fwdResult = doLikelihoodScan(segdf.dtheta, segdf.segradlengths, segdf.cumlenfwd, fwd_fit=True, mom_dep_const=mom_dep_const, pid=pid, ang_resol=ang_resol)
    bwdResult = doLikelihoodScan(segdf.dtheta, segdf.segradlengths, segdf.cumlenbwd, fwd_fit=False, mom_dep_const=mom_dep_const, pid=pid, ang_resol=ang_resol)
    return fwdResult, bwdResult


def fitMcsPrepare(traj):
    breakpoints, segradlengths, cumseglens, resranges = breakTrajInSegments(traj)
    #print("[TrajectoryMCSFitter::fitMcsPrepare] start")
    #print(breakpoints)
    #print(segradlengths)
    #print(cumseglens)
    #print(resranges)

    if (len(segradlengths)<2):
        logging.info("not enough segments for MCS fit, skipping traj")
        return None
    
    dtheta = []
    dtheta_xz = []
    dtheta_yz = []
    dtheta_xz_rot = []
    dtheta_yz_rot = []
    dtheta_xz_prime = []
    dtheta_yz_prime = []
    for p in range(len(segradlengths)):
        pcdir1 = linearRegression(traj, breakpoints[p], breakpoints[p+1]) # normalized to 1
        if pcdir1 is None:
            return None
        if (p > 0):
            if ((segradlengths[p] < -100) | (segradlengths[p-1] < -100)):
                dtheta.append(-999)
                dtheta_xz.append(-999)
                dtheta_yz.append(-999)
                dtheta_xz_rot.append(-999)
                dtheta_yz_rot.append(-999)
                dtheta_xz_prime.append(-999)
                dtheta_yz_prime.append(-999)
                logging.info("bad breakpoint")

            # TODO: cut on hits
            # elif ((~breakpointsgood[p]) | (~breakpointsgood[p-1])):
            #     dtheta.append(-999)
            #     logging.info("bad breakpoint")

            else:
                cosval = np.dot(pcdir0, pcdir1)
                dt = 1e3*np.arccos(cosval) # should we try to use expansion for small angles?
                dtheta.append(dt)

                # Project to x–z plane (x = index 0, z = index 2)
                v0_xz = pcdir0[[0, 2]]
                v1_xz = pcdir1[[0, 2]]
                v0_xz = v0_xz / np.linalg.norm(v0_xz)
                v1_xz = v1_xz / np.linalg.norm(v1_xz)
                angle_xz_signed = 1e3*(np.arctan2(v0_xz[1], v0_xz[0]) - np.arctan2(v1_xz[1], v1_xz[0]))
                dtheta_xz.append(angle_xz_signed)
                
                v0_yz = pcdir0[[1, 2]]
                v1_yz = pcdir1[[1, 2]]
                v0_yz = v0_yz / np.linalg.norm(v0_yz)
                v1_yz = v1_yz / np.linalg.norm(v1_yz)
                angle_yz_signed = 1e3*(np.arctan2(v0_yz[1], v0_yz[0]) - np.arctan2(v1_yz[1], v1_yz[0]))
                dtheta_yz.append(angle_yz_signed)

                pcz = np.array([0, 0, 1])
                cos_pcdir0_z = np.dot(pcdir0, pcz)
                pcdir0_X_z = np.cross(pcdir0, pcz)
                pcdir0_X_z = pcdir0_X_z / np.linalg.norm(pcdir0_X_z)
                rotated_pcdir1 = rotate_vector(pcdir1, pcdir0_X_z, cos_pcdir0_z)
                angle_xz_rot = 1e3*np.arctan2(rotated_pcdir1[1], rotated_pcdir1[2])
                angle_yz_rot = 1e3*np.arctan2(rotated_pcdir1[0], rotated_pcdir1[2])
                dtheta_xz_rot.append(angle_xz_rot)
                dtheta_yz_rot.append(angle_yz_rot)

                # Follow angle definitions in https://arxiv.org/pdf/2605.03048
                ## z' is muon dir of the previous segment
                this_zprime = pcdir0 / np.linalg.norm(pcdir0)
                ## y' = z' X x & x' = -z' X y'
                this_xvec = np.array([1, 0, 0])
                this_yprime = np.cross(this_zprime, this_xvec)
                this_yprime /= np.linalg.norm(this_yprime)
                this_xprime = np.cross(-1. * this_zprime, this_yprime)
                ## now, collect pcdir1's compoments in x',y',z'
                pcdir1_xprime = np.dot(this_xprime, pcdir1)
                pcdir1_yprime = np.dot(this_yprime, pcdir1)
                pcdir1_zprime = np.dot(this_zprime, pcdir1)
                theta_xz_prime = 1e3*np.arctan2(pcdir1_xprime, pcdir1_zprime)
                theta_yz_prime = 1e3*np.arctan2(pcdir1_yprime, pcdir1_zprime)
                dtheta_xz_prime.append(theta_xz_prime)
                dtheta_yz_prime.append(theta_yz_prime)

        pcdir0 = pcdir1
    
    cumseglens = np.array(cumseglens)
    cumlenfwd = cumseglens[:-2]
    cumlenbwd = (cumseglens[-1]-cumseglens[::-1])[:-2]

    resranges = np.array(resranges)
    #print(resranges)
    resranges_fwd = resranges[:-2] # R.R. at the begining of seg0
    resranges_fwd_brk_p = resranges[1:-1] # R.R. at breaking points of the two segments (seg0 ~ pcdir0, seg1 ~ pcdir1)
    resranges_fwd_next_seg_end = resranges[2:] # R.R. at the end of the seg 1
    resranges_bwd = (resranges[0]-resranges[::-1])[:-2]
    resranges_fwd_P = muonRangeP(resranges_fwd)/1e3 # [GeV]
    resranges_mid_seg0 = (resranges_fwd + resranges_fwd_brk_p) / 2. # R.R. on midpoint of seg0
    resranges_mid_seg1 = (resranges_fwd_brk_p + resranges_fwd_next_seg_end) / 2. # R.R. on midpoint of seg1
    resranges_mid_seg0_P = muonRangeP(resranges_mid_seg0)/1e3 # [GeV] # P_range for midpoint in seg0
    resranges_mid_seg1_P = muonRangeP(resranges_mid_seg1)/1e3 # [GeV] # P_range for midpoint in seg1
    resranges_mid_mean_P = (resranges_mid_seg0_P + resranges_mid_seg1_P) / 2. # [GeV], mean of the P_range at the two midpoints
    segdf = pd.DataFrame( {"breakpoint": range(len(dtheta))
                           , "dtheta": dtheta
                           , "dtheta_xz": dtheta_xz
                           , "dtheta_yz": dtheta_yz
                           , "angle_xz_rot": dtheta_xz_rot
                           , "angle_yz_rot": dtheta_yz_rot
                           , "dtheta_xz_prime": dtheta_xz_prime
                           , "dtheta_yz_prime": dtheta_yz_prime
                           , "segradlengths": segradlengths[:-1]
                           , "cumlenfwd": cumlenfwd
                           , "cumlenbwd": cumlenbwd
                           , "resranges_fwd": resranges_fwd
                           , "resranges_fwd_brk_p": resranges_fwd_brk_p
                           , "resranges_fwd_next_seg_end": resranges_fwd_next_seg_end
                           , "resranges_bwd": resranges_bwd
                           , "resranges_fwd_P": resranges_fwd_P
                           , "resranges_mid_seg0_P": resranges_mid_seg0_P
                           , "resranges_mid_seg1_P": resranges_mid_seg1_P
                           , "resranges_mid_mean_P": resranges_mid_mean_P
                           } )

    #print(segdf)
    return segdf


def breakTrajInSegments(traj):
    trajlen = traj.length
    nseg = max(MIN_NSEG, int(trajlen)//SEG_LEN)
    seglen = trajlen/nseg

    cumseglens = [0] #first segment has zero cumulative length from previous segments
    thislen = 0.

    next_valid = traj.first_point()
    breakpoints = [next_valid]
    segradlengths = []
    resranges = [traj.rr_at_point(next_valid)]
    pos0 = traj.location_at_point(next_valid)

    next_valid = traj.next_point(next_valid)

    #print("[TrajectoryMCSFitter::breakTrajInSegments] next_valid")
    #print(breakpoints)
    #print(next_valid)

    npoints = 0
    while next_valid is not None:
        pos1 = traj.location_at_point(next_valid)
        thislen += np.linalg.norm(pos1-pos0)
        pos0 = pos1
        npoints += 1
        if (thislen >= seglen):
            breakpoints.append(next_valid)
            if (npoints >= MIN_NHITS_PER_SEG):
                segradlengths.append(thislen*LAR_RADL_INV)
            else:
                logging.info("not enough hits in segment")
                segradlengths.append(-999.)
            cumseglens.append(cumseglens[-1]+thislen)
            resranges.append(traj.rr_at_point(next_valid))
            thislen = 0.
            npoints = 0.
        
        next_valid = traj.next_point(next_valid)

    if (thislen > 0): # end of trarck
        last_valid = traj.last_point()
        
        end_boundary = last_valid[:-1] + (last_valid[-1] + 1,)
        breakpoints.append(end_boundary)
        
        segradlengths.append(thislen*LAR_RADL_INV)
        cumseglens.append(cumseglens[-1]+thislen)
        resranges.append(traj.rr_at_point(last_valid))

    return breakpoints, segradlengths, cumseglens, resranges


def doLikelihoodScan(dtheta, segradlengths, cumlen, fwd_fit=True, mom_dep_const=MOM_DEP_CONST, pid=PID, ang_resol=ANG_RESOL, factors_list=None):
    best_idx  = -1
    best_logL = float("inf")
    best_p    = -1.0
    vlogL = []
    p_test = P_MIN
    while p_test <= P_MAX:
        logL = mcsLikelihood(p_test, dtheta, segradlengths, cumlen, fwd_fit, mom_dep_const, pid, ang_resol, factors_list)
        if logL < best_logL:
            best_p    = p_test
            best_logL = logL
            best_idx  = len(vlogL)
        vlogL.append(logL)
        p_test += P_STEP
        
    lunc = -1.0
    if best_idx > 0:
        for j in range(best_idx-1, -1, -1):
            dLL = vlogL[j] - vlogL[best_idx]
            if dLL < 0.5:
                lunc = (best_idx - j) * P_STEP
            else:
                break
    runc = -1.0
    if best_idx < len(vlogL) - 1:
        for j in range(best_idx+1, len(vlogL)):
            dLL = vlogL[j] - vlogL[best_idx]
            if dLL < 0.5:
                runc = (j - best_idx) * P_STEP
            else:
                break
                
    return {"p": best_p, "pUnc": max(lunc, runc), "logL": best_logL}

def linearRegression(traj, first_p, last_p): 
    npoints = 0
    points = []
    
    # Extract the shared tuple prefix and the integer indices
    prefix = first_p[:-1]
    start_idx = first_p[-1]
    end_idx = last_p[-1]

    # First loop: Collect points
    #for hit_idx in range(start_idx, end_idx, -1):
    for hit_idx in range(start_idx, end_idx):
        idx = prefix + (hit_idx,)  # Reconstruct tuple
        points.append(traj.location_at_point(idx))
        npoints += 1
        
    if npoints == 0:
        return None
        
    sumpoints = np.array(points).sum(axis=0)
    avgpos = sumpoints/npoints
    norm = 1/npoints

    m = np.zeros((3,3))
    
    # Second loop: Build covariance matrix
    #for hit_idx in range(start_idx, end_idx, -1):
    for hit_idx in range(start_idx, end_idx):
        idx = prefix + (hit_idx,)  # Reconstruct tuple
        p = traj.location_at_point(idx)
        xxw0 = p[0]-avgpos[0]
        yyw0 = p[1]-avgpos[1]
        zzw0 = p[2]-avgpos[2]
        m[0, 0] += xxw0*xxw0*norm
        m[0, 1] += xxw0*yyw0*norm
        m[0, 2] += xxw0*zzw0*norm
        m[1, 0] += yyw0*xxw0*norm
        m[1, 1] += yyw0*yyw0*norm
        m[1, 2] += yyw0*zzw0*norm
        m[2, 0] += zzw0*xxw0*norm
        m[2, 1] += zzw0*yyw0*norm
        m[2, 2] += zzw0*zzw0*norm
        
    me = TMatrixDSymEigen(m)
    eigenval = me.get_eigen_values()
    eigenvec = me.get_eigen_vectors().T

    maxevalidx = np.argmax(eigenval)
    pcdir = np.array([eigenvec[0, maxevalidx], eigenvec[1, maxevalidx], eigenvec[2, maxevalidx]])

    # Notice we can still use the original first_p tuple here!
    first_dir = traj.direction_at_point(first_p)
                                                                                                                                                                                                                                                                                 
    if first_dir is None:
        return None
    if (np.dot(first_dir,pcdir)<0.):
        pcdir*=-1.
        
    return pcdir

def mcsLikelihood(p, theta_list, segradlengths_list, cumlen, fwd_fit=True, mom_dep_const=MOM_DEP_CONST, pid=PID, ang_resol=ANG_RESOL, factors_list=None):
    beg = 0 if fwd_fit else len(theta_list) - 1
    end = len(theta_list) if fwd_fit else -1
    incr = 1 if fwd_fit else -1

    m = mass(pid)
    p2, m2 = p*p, m*m
    Etot = np.sqrt(p2 + m2)
    Eij2 = 0.0

    FIXED_TERM = 0.5 * np.log(2.0 * np.pi) 
    result = 0.0

    for i in range(beg, end, incr):
        # if segradlengths_list[i] < 0:
        #     continue
        # TODO: need to specify as -999 for 2D fit
        if theta_list[i] < 0:
            continue

        if E_LOSS_MODE == 1: # ELoss mode: MIP (constant)
            kcal = 0.002105
            Eij = Etot - kcal * cumlen[i] # energy at this segment
            Eij2 = Eij * Eij
        else:
            Eij = getE(Etot, cumlen[i], m) # Non constant energy loss distribution
            Eij2 = Eij * Eij

        if Eij2 <= m2:
            result = float("inf")
            break

        pij = np.sqrt(Eij2 - m2) 
        tH0 = highlandFormula(pij, segradlengths_list[i], m, mom_dep_const)

        if factors_list is not None:
            tH0 *= factors_list[i]
        else: # 3D
            dim_factor = 2

        rms = np.sqrt(dim_factor * (tH0*tH0 + ang_resol*ang_resol))
        if rms == 0.0:
            logging.error("RMS cannot be zero")
            return float("inf")
        arg = theta_list[i] / rms
        chi2 = (np.log(rms) + 0.5*arg*arg + FIXED_TERM)
        result += chi2
    return result

def hlPrediction(p, theta_list, segradlengths_list, cumlen, fwd_fit=True, mom_dep_const=MOM_DEP_CONST, pid=PID, ang_resol=ANG_RESOL, factors_list=None):
    beg = 0 if fwd_fit else len(theta_list) - 1
    end = len(theta_list) if fwd_fit else -1
    incr = 1 if fwd_fit else -1

    m = mass(pid)
    p2, m2 = p*p, m*m
    Etot = np.sqrt(p2 + m2)
    Eij2 = 0.0

    FIXED_TERM = 0.5 * np.log(2.0 * np.pi) 
    result = 0.0

    tH0_list = []
    rms_list = []
    pij_list = []
    for i in range(beg, end, incr):
        # if segradlengths_list[i] < 0:
        #     continue
        # TODO: need to specify as -999 for 2D fit
        if theta_list[i] < 0:
            continue

        if E_LOSS_MODE == 1: # ELoss mode: MIP (constant)
            kcal = 0.002105
            Eij = Etot - kcal * cumlen[i] # energy at this segment
            Eij2 = Eij * Eij
        else:
            Eij = getE(Etot, cumlen[i], m) # Non constant energy loss distribution
            Eij2 = Eij * Eij

        if Eij2 <= m2:
            result = float("inf")
            break

        pij = np.sqrt(Eij2 - m2) 
        tH0 = highlandFormula(pij, segradlengths_list[i], m, mom_dep_const)

        if factors_list is not None:
            tH0 *= factors_list[i]
        else: # 3D
            dim_factor = 2

        rms = np.sqrt(dim_factor * (tH0*tH0 + ang_resol*ang_resol))
        if rms == 0.0:
            logging.error("RMS cannot be zero")
            return float("inf")
        arg = theta_list[i] / rms
        chi2 = (np.log(rms) + 0.5*arg*arg + FIXED_TERM)
        result += chi2

        tH0_list.append(tH0)
        rms_list.append(rms)
        pij_list.append(pij)

    ret = {"chi2": result, "tH0": tH0_list, "rms": rms_list, "pij": pij_list}
    return ret



def momentumDependentConstant(pij):
    # https://arxiv.org/abs/1703.06187
    a = 0.1049
    c = 11.0038
    return ((a/(pij*pij)) + c)


def highlandFormula(pij, segradlength, m, mom_dep_const=False):
    TUNED_HL_TERM = 11.0038 # https://arxiv.org/abs/1703.06187
    HL_TERM2 = 0.038
    m2 = m*m
    beta = np.sqrt(1.0 - ((m2) / (pij*pij + m2)))
    tH0 = ((momentumDependentConstant(pij) if mom_dep_const else TUNED_HL_TERM) / (pij*beta))*(1.0 + HL_TERM2*np.log(segradlength))*np.sqrt(segradlength)
    return tH0


def energyLossLandau(mass2, e2, x):
    # eq. (33.11) in http://pdg.lbl.gov/2016/reviews/rpp2016-rev-passage-particles-matter.pdf (except density correction is ignored)
    if x <= 0.0:
        return 0.0
    
    Iinv2 = 1.0 / (188.0E-6 * 188.0E-6)
    matConst = 1.4 * 18.0 / 40.0 # density*Z/A
    me = 0.511 # MeV
    kappa = 0.307075
    j = 0.200
    
    beta2 = (e2 - mass2) / e2
    gamma2 = 1.0 / (1.0 - beta2)
    epsilon = 0.5 * kappa * x * matConst / beta2
    
    return 0.001 * epsilon * (np.log(2.0 * me * beta2 * gamma2 * epsilon * Iinv2) + j - beta2)


def energyLossBetheBloche(mass, e2):
    # stolen, mostly, from GFMaterialEffects.
    Iinv = 1.0 / 188.0E-6
    matConst = 1.4 * 18.0 / 40.0 # density*Z/A
    me = 0.511
    kappa = 0.307075
    
    beta2 = (e2 - mass * mass) / e2
    gamma2 = 1.0 / (1.0 - beta2)
    massRatio = me / mass
    argument = (2.0 * me * gamma2 * beta2 * Iinv) * (1 + 2 * (gamma2)**0.5 * massRatio + massRatio**2)**0.5
    
    dedx = kappa * matConst / beta2
    
    if mass == 0.0:
        return 0.0
    if argument <= np.exp(beta2):
        dedx = 0.0
    else:
        dedx *= (np.log(argument) - beta2) * 1.0E-3 # Bethe-Bloch, converted to GeV/cm
        if dedx < 0.0:
            dedx = 0.0
    
    return dedx


# def dEdx_Bethe_Bloch(KE, mass):
#     gamma = (KE/mass)+1.0
#     beta = np.sqrt(1-(1.0/(gamma*gamma)))
#     Wmax = (2.0 * Me * pow(beta * gamma, 2)) / (1.0 + 2.0 * Me * (gamma / mass) + pow((Me / mass),2))
#     delta = density_correction(beta, gamma)
#     # density correction
#     f = rho * K * (18.0 / A) * pow(1. / beta, 2)
#     a0 = 0.5 * np.log(2.0 * Me * pow(beta * gamma, 2) * Wmax / (I * I))
#     this_dEdx = f * ( a0 - pow(beta, 2) - delta / 2.0) # [MeV/cm]

#     return this_dEdx;


def getE(initial_E, length_travelled, m, eLossMode=0):
    nElossSteps = 10
    step_size = length_travelled / nElossSteps
    current_E = initial_E
    m2 = m * m
    
    for i in range(nElossSteps):
        if eLossMode == 2:
            dedx = energyLossBetheBloche(m, current_E * current_E)
            current_E -= (dedx * step_size)
        else:
            # MPV of Landau energy loss distribution
            current_E -= energyLossLandau(m2, current_E * current_E, step_size)
        
        if current_E <= m:
            # print("WARNING: current_E less than mu mass. it is", current_E) # p_min hypothesis is small, non-fatal warning
            return 0.0
    
    return current_E


def mass(pid): # GeV
    if (np.abs(pid)==13):
        MASS_MUON = 0.105658
        return MASS_MUON 
    if (np.abs(pid)==211):
        MASS_PION = 0.139570
        return MASS_PION

def rotate_vector(v, axis, cos_theta):
    """
    Rotate vector `v` around `axis` by `theta_rad` radians using Rodrigues' formula.
    
    Parameters:
    - v: np.ndarray shape (3,) - vector to rotate
    - axis: np.ndarray shape (3,) - rotation axis
    - theta_rad: float - rotation angle in radians
    
    Returns:
    - rotated vector: np.ndarray shape (3,)
    """
    axis = axis / np.linalg.norm(axis)  # Ensure unit vector
    v = np.asarray(v)
    sin_theta = np.sqrt(1 - cos_theta * cos_theta)
    cross = np.cross(axis, v)
    dot = np.dot(axis, v)
    
    v_rot = v * cos_theta + cross * sin_theta + axis * dot * (1 - cos_theta)
    return v_rot
    
# def getRangeE(track_length, m, eLossMode=0):
#     nElossSteps = 10
#     step_size = track_length / nElossSteps
#     current_E = 0
#     m2 = m * m
    
#     for i in range(nElossSteps):
#         if eLossMode == 2:
#             dedx = energyLossBetheBloche(m, current_E * current_E)
#             current_E -= (dedx * step_size)
#         else:
#             # MPV of Landau energy loss distribution
#             current_E -= energyLossLandau(m2, current_E * current_E, step_size)
        
#         if current_E <= m:
#             # print("WARNING: current_E less than mu mass. it is", current_E) # p_min hypothesis is small, non-fatal warning
#             return 0.0
    
#     return current_E
