import numpy as np
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

mass_muon = 105.658*1e-3
mass_chpion = 139.57*1e-3
mass_proton = 938.272*1e-3

class MCSFIT:
    def __init__(self, angle_mode = 0, dim_factor = 2., ang_trunc = False,
                 kappa_a = 0.105, kappa_c = 11.004, ANG_RESOL=3.0, cut_seg = 2,
                 P_MIN=0.01, P_MAX=7.5, P_STEP=0.01,
                 MOM_DEP_CONST=True, mass=mass_muon, E_LOSS_MODE = 0, epslion = 0.038):
        self.density_C = 5.2146
        self.density_y0 = 0.2
        self.density_y1 = 3.0
        self.density_a = 0.19559
        self.density_k = 3.0
        self.ln10 = np.log(10)
        self.ang_trunc_sigma = 2.5
        
        self.angle_mode = angle_mode
        self.dim_factor = dim_factor
        self.ang_trunc = ang_trunc
        self.P_MIN = P_MIN
        self.P_MAX = P_MAX
        self.P_STEP = P_STEP
        self.MOM_DEP_CONST = MOM_DEP_CONST
        self.ANG_RESOL = ANG_RESOL
        self.E_LOSS_MODE = E_LOSS_MODE
        self.mass = mass
        self.kappa_a = kappa_a
        self.kappa_c = kappa_c
        self.epslion = epslion
        self.cut_seg = cut_seg

    def energyLossLandau(self, e2, x):
        if x <= 0.0:
            return 0.0

        Iinv2 = 1.0 / (188.0E-6)**2
        matConst = 1.4 * 18.0 / 40.0
        me = 0.511
        kappa = 0.307075
        j = 0.200

        mass2 = self.mass**2
        beta2 = (e2 - mass2) / e2
        gamma2 = 1.0 / (1.0 - beta2)
        epsilon = 0.5 * kappa * x * matConst / beta2

        return 1e-3*epsilon * (np.log(2.0 * me * beta2 * gamma2 * epsilon * Iinv2) + j - beta2)

    def energyLossLandauWithDensityCorrection(self, e2, x):
        if x <= 0.0:
            return 0.0

        Iinv2 = 1.0 / (188.0E-6)**2
        matConst = 1.4 * 18.0 / 40.0
        me = 0.511
        kappa = 0.307075
        j = 0.200

        mass2 = self.mass**2
        beta2 = (e2 - mass2) / e2
        gamma2 = 1.0 / (1.0 - beta2)
        epsilon = 0.5 * kappa * x * matConst / beta2

        betagamma = np.sqrt(beta2 * gamma2)
        density_y = np.log10(betagamma)

        if density_y > self.density_y1:
            delta = 2.0 * self.ln10 * density_y - self.density_C
        elif density_y < self.density_y0:
            delta = 0.0
        else:
            delta = 2.0 * self.ln10 * density_y - self.density_C + self.density_a * ((self.density_y1 - density_y) ** self.density_k)
        
        return 1e-3*epsilon * (np.log(2.0 * me * beta2 * gamma2 * epsilon * Iinv2) + j - beta2 - delta)
    
    def energyLossBetheBloch(self, e2):
        Iinv = 1.0 / 188.0E-6
        matConst = 1.4 * 18.0 / 40.0
        me = 0.511
        kappa = 0.307075

        beta2 = (e2 - self.mass * self.mass) / e2
        gamma2 = 1.0 / (1.0 - beta2)
        massRatio = me * 1e-3 / self.mass
        argument = (2.0 * me * gamma2 * beta2 * Iinv) / np.sqrt(1 + 2 * np.sqrt(gamma2) * massRatio + massRatio**2)

        dedx = kappa * matConst / beta2

        if self.mass == 0.0 or argument <= np.exp(beta2):
            return 0.0
        else:
            dedx *= (np.log(argument) - beta2) * 1e-3
            return max(dedx, 0.0)

    def energyLossBetheBlochWithDensityCorrection(self, e2):
        Iinv = 1.0 / 188.0E-6
        matConst = 1.36 * 18.0 / 39.948
        me = 0.511
        kappa = 0.307075

        beta2 = (e2 - self.mass * self.mass) / e2
        gamma2 = 1.0 / (1.0 - beta2)
        massRatio = me * 1e-3 / self.mass
        argument = (2.0 * me * gamma2 * beta2 * Iinv) / np.sqrt(1 + 2 * np.sqrt(gamma2) * massRatio + massRatio**2)

        dedx = kappa * matConst / beta2

        betagamma = np.sqrt(beta2 * gamma2)
        density_y = np.log10(betagamma)

        if density_y > self.density_y1:
            delta = 2.0 * self.ln10 * density_y - self.density_C
        elif density_y < self.density_y0:
            delta = 0.0
        else:
            delta = 2.0 * self.ln10 * density_y - self.density_C + self.density_a * ((self.density_y1 - density_y) ** self.density_k)

        density_corr = delta / 2.0

        if self.mass == 0.0 or argument <= np.exp(beta2):
            return 0.0
        else:
            dedx *= (np.log(argument) - beta2 - density_corr) *1e-3
            return max(dedx, 0.0)

    def momentumDependentConstant(self, pij):
        ## == parameterization of MicroBooNE is a bit tricky
        ## ==== Input momentum should be in [GeV] unit, and output constant has [MeV] unit...
        return (self.kappa_a / (pij**2)) + self.kappa_c

    def highlandFormula(self, pij, segradlength):
        TUNED_HL_TERM = self.kappa_a
        m2 = self.mass * self.mass
        beta = np.sqrt(1.0 - (m2 / (pij**2 + m2)))
        kappa = self.momentumDependentConstant(pij) if self.MOM_DEP_CONST else TUNED_HL_TERM
        tH0 = (kappa / (1e3 * pij * beta)) * (1.0 + self.epslion * np.log(segradlength)) * np.sqrt(segradlength)
        return tH0*1e3
    
    def getE(self, initial_E, length_travelled):
        nElossSteps = 10
        step_size = length_travelled / nElossSteps
        current_E = initial_E
        m2 = self.mass * self.mass

        for i in range(nElossSteps):
            if self.E_LOSS_MODE == 2:
                dedx = self.energyLossBetheBloch(current_E * current_E)
                current_E -= (dedx * step_size)
            elif self.E_LOSS_MODE == 0:
                # MPV of Landau energy loss distribution
                #current_E -= self.energyLossLandau(current_E * current_E, step_size)
                current_E -= self.energyLossLandauWithDensityCorrection(current_E * current_E, step_size)
            else:
                dedx = self.energyLossBetheBlochWithDensityCorrection(current_E * current_E)
                current_E -= (dedx * step_size)

            if current_E <= self.mass:
                # print("WARNING: current_E less than mu mass. it is", current_E) # p_min hypothesis is small, non-fatal warning
                return 0.0

        return current_E
    def mcsLikelihood(self, p, theta_list, segradlengths_list, cumlen, fwd_fit=True, ang_factor=1., dim_factor=1.):
        beg = 0 if fwd_fit else len(theta_list) - 1
        end = len(theta_list) if fwd_fit else -1
        incr = 1 if fwd_fit else -1

        m = self.mass
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

            if self.E_LOSS_MODE == 1: # ELoss mode: MIP (constant)
                kcal = 0.002105
                Eij = Etot - kcal * cumlen[i] # energy at this segment
                Eij2 = Eij * Eij
            else:
                Eij = self.getE(Etot, cumlen[i]) # Non constant energy loss distribution        
                Eij2 = Eij * Eij

            if Eij2 <= m2:
                #print("Eij2 <= m2") # debug
                result = float("inf")
                break

            pij = np.sqrt(Eij2 - m2)
            tH0 = self.highlandFormula(pij, segradlengths_list[i])
            tH0 *= ang_factor
            
            rms = np.sqrt(dim_factor * (tH0*tH0 + self.ANG_RESOL*self.ANG_RESOL))
            if rms == 0.0:
                logging.error("RMS cannot be zero")
                return float("inf")
            arg = theta_list[i] / rms
            chi2 = (np.log(rms) + 0.5*arg*arg + FIXED_TERM)
            
            this_truc_lim = np.sqrt(dim_factor * (tH0*tH0*self.ang_trunc_sigma*self.ang_trunc_sigma + self.ANG_RESOL*self.ANG_RESOL))
            if self.ang_trunc and theta_list[i] > this_truc_lim:
                chi2 = 0
            result += chi2
        return result

    def doLikelihoodScan(self, dtheta, angle_xz_rot, angle_yz_rot, segradlengths, cumlen, fwd_fit=True, factors_list=None):
        best_idx  = -1
        best_logL = float("inf")
        best_p    = -1.0
        vlogL = []
        p_test = self.P_MIN
        while p_test <= self.P_MAX:
            if self.angle_mode == 0: ## 3D angle with Gaussian likelihood
                logL = self.mcsLikelihood(p_test, dtheta, segradlengths, cumlen, fwd_fit, 1., self.dim_factor)
            else: ## 2D angle with Gaussian likelihood
                xz_logL = self.mcsLikelihood(p_test, angle_xz_rot, segradlengths, cumlen, fwd_fit, 1., 1.)
                yz_logL = self.mcsLikelihood(p_test, angle_yz_rot, segradlengths, cumlen, fwd_fit, 1., 1.)
                logL = xz_logL + yz_logL

            if logL < best_logL:
                best_p    = p_test
                best_logL = logL
                best_idx  = len(vlogL)
            #print('p_test: %f, logL: %f' %(p_test, logL)) ## debug
            vlogL.append(logL)
            p_test += self.P_STEP

        lunc = -1.0
        if best_idx > 0:
            for j in range(best_idx-1, -1, -1):
                dLL = vlogL[j] - vlogL[best_idx]
                if dLL < 0.5:
                    lunc = (best_idx - j) * self.P_STEP
                else:
                    break
        runc = -1.0
        if best_idx < len(vlogL) - 1:
            for j in range(best_idx+1, len(vlogL)):
                dLL = vlogL[j] - vlogL[best_idx]
                if dLL < 0.5:
                    runc = (j - best_idx) * self.P_STEP
                else:
                    break

        return {"p": best_p, "pUnc": max(lunc, runc), "logL": best_logL}
    def fitMcs(self, segdf):
        # TODO: truncate segments near track end to simulate exiting tracks     
        if self.cut_seg > 0:
            if len(segdf) > (self.cut_seg+1):
                segdf = segdf[:-self.cut_seg]
            else:
                logging.info("not enough segments to truncate at track end, skipping traj")
                return None

        fwdResult = self.doLikelihoodScan(segdf.dtheta, segdf.angle_xz_rot, segdf.angle_yz_rot, segdf.segradlengths, segdf.cumlenfwd, fwd_fit=True)
        bwdResult = self.doLikelihoodScan(segdf.dtheta, segdf.angle_xz_rot, segdf.angle_yz_rot, segdf.segradlengths, segdf.cumlenbwd, fwd_fit=False)
        return fwdResult, bwdResult, len(segdf)

