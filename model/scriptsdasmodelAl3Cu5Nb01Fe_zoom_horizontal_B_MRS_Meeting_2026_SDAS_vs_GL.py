from numpy import trapezoid as trapz, abs, zeros, sqrt, exp, log, pi, cos, sin
import operator as op
import numpy as np
import matplotlib.pyplot as plt
from math  import exp, gamma, log
from scipy.optimize import brentq
import scipy.optimize as so
from midpoint import midpoint

# Alloy Al-0.8Si-0.6Mg-0.20Fe Horizontal

global Error

vc01 = np.zeros(300)
vk01 = np.zeros(300)
vbeta1 = np.zeros(300)
vgama1 = np.zeros(300)
vc02 = np.zeros(300)
vk02 = np.zeros(300)
vbeta2 = np.zeros(300)
vgama2 = np.zeros(300)
vc03 = np.zeros(300)
vk03 = np.zeros(300)
vbeta3 = np.zeros(300)
vgama3 = np.zeros(300)

vFo1 = np.zeros(300)
vFo2 = np.zeros(300)
vFo3 = np.zeros(300)
vtSL = np.zeros(300)


def Diffusion_coefficient(T,D0,Qd,R):
 
    Diff = D0 * exp(-Qd/(R*T))

    return Diff 


xsc = np.zeros(20)
ysc = np.zeros(20)
yscp = np.zeros(20)
yscm = np.zeros(20)
error_ysc = np.zeros(20)

# SDAS Experimental Data




xsc[0] = 20.19709
ysc[0] = 20.74
xsc[1] = 37.06885
ysc[1] = 27.542
xsc[2] = 52.87807
ysc[2] = 32.589
xsc[3] = 68.03454
ysc[3] = 36.011
xsc[4] = 82.72341
ysc[4] = 39.579
xsc[5] = 97.05009
ysc[5] = 42.509
xsc[6] = 111.08242
ysc[6] = 46.543
xsc[7] = 124.8676
ysc[7] = 49.284
xsc[8] = 138.44027
ysc[8] = 51.264
xsc[9] = 151.82691
ysc[9] = 49.523
xsc[10] = 178.12147
ysc[10] = 56.547
xsc[11] = 203.8758
ysc[11] = 59.444
xsc[12] = 229.17651
ysc[12] = 61.957
xsc[13] = 254.0872
ysc[13] = 64.785
xsc[14] = 278.65645
ysc[14] = 66.171
xsc[15] = 302.92253
ysc[15] = 68.273



yscp[0] = 10.244
yscp[1] = 8.128
yscp[2] = 21.214
yscp[3] = 32.701
yscp[4] = 24.239
yscp[5] = 11.782
yscp[6] = 16.846
yscp[7] = 19.652
yscp[8] = 18.737
yscp[9] = 17.184
yscp[10] = 33.294
yscp[11] = 15.457
yscp[12] = 22.357
yscp[13] = 19.529
yscp[14] = 19.455
yscp[15] = 26.935



yscm[0] = 5.798
yscm[1] = 9.594
yscm[2] = 12.934
yscm[3] = 14.001
yscm[4] = 13.577
yscm[5] = 16.987
yscm[6] = 11.701
yscm[7] = 19.552
yscm[8] = 14.585
yscm[9] = 15.962
yscm[10] = 23.224
yscm[11] = 18.638
yscm[12] = 18.22
yscm[13] = 20.387
yscm[14] = 14.48
yscm[15] = 17.128




"""

xsc[0] = 48.8 
ysc[0] = 31.718275 
xsc[1] = 85.4 
ysc[1] = 33.6492615 
xsc[2] = 128 
ysc[2] = 36.525645 
xsc[3] = 208
ysc[3] = 42.912604  
xsc[4] = 269 
ysc[4] = 50.38717 
xsc[5] = 328
ysc[5] = 54.3  




yscp[0] = 5.918675
yscp[1] = 2.94236
yscp[2] = 1.96885
yscp[3] = 2.9276
yscp[4] = 7.37857
yscp[5] = 7.04


yscm[0] = 4.65143
yscm[1] = 1.89694
yscm[2] = 1.77025
yscm[3] = 2.6122
yscm[4] = 8.88003
yscm[5] = 4.92

"""

for i in range(12):
    error_ysc[i] = 0.1 * ysc[i]

R = 8.31451

# Input data
TL = 936.58 #922.608 #925.32 #923.86 #866.5 #933.15   # Verify after
T_Al = TL
TF = 933.15
Ts = 817.75 #855.00 #525.0 + 273.15
DeltaH = 335300.0 #293900.0 #335300 #260000.0 
rhos = 2737.04 #2620.24 #2555.72 #2556.0
rhol = 2527.05 #2433.79 #2378.33 #2386.0
cps  = 1270.0
sigmasl = 1.3761 #1.209 #0.914
gamma_0 = 0.296  #0.168 #0.154
mu0 = 1.6
lambda0 = -5.0
rho0 = 2.56e-6
mu = 2.594e10
lambda_ = 5.034e10
sigma = sigmasl



Delta_Sv_hom = rhos * DeltaH/TL


# Data solute Cu
Dol_Cu = 1.06e-7
Qol_Cu = 24000
Dos_Cu = 6.5e-5
Qos_Cu = 136100
c01=3.0
cf1= 28.816 
DL1= Diffusion_coefficient(TL,Dol_Cu,Qol_Cu,R) 
DS1= Diffusion_coefficient(Ts,Dos_Cu,Qos_Cu,R) 
k01= 0.1308  
mL1= -3.64103  



# Data solute Nb
Dol_Nb = 1.1e-5 #9.9e-5 #1.06e-7
Qol_Nb = 28400.0 #71600.0 #24000
Dos_Nb = 2.0e-4 #6.23e-6 #6.5e-5
Qos_Nb = 134500 #152700 #115000.0 #136100
c02= 5.0 #0.5
cf2= 11.02 #34.2 #28.816 
DL2= Diffusion_coefficient(TL,Dol_Nb,Qol_Nb,R) 
DS2= Diffusion_coefficient(Ts,Dos_Nb,Qos_Nb,R) 
k02= 2.049 #0.36 #0.1308  
mL2= -2.96 #-5.1 #-3.64103  

"""
# Data solute Si
Dol_Si = 1.34E-07
Qol_Si = 30000.0 
Dos_Si = 2.48E-04
Qos_Si = 137000.0
c02= 0.8 #4.0
cf2= 12.7 #10.7275 
DL2= Diffusion_coefficient(TL,Dol_Si,Qol_Si,R) 
DS2= Diffusion_coefficient(Ts,Dos_Si,Qos_Si,R) 
k02= 0.11 #0.1075
mL2= -6.2 #-7.198322069
"""

# Data solute Fe
Dol_Fe = 2.34E-07
Qol_Fe = 35000.0
Dos_Fe = 5.30E-03
Qos_Fe = 183400.0
c03= 0.2  #0.11
cf3= 1.81 #1.655
DL3= Diffusion_coefficient(TL,Dol_Fe,Qol_Fe,R) 
DS3= Diffusion_coefficient(Ts,Dos_Fe,Qos_Fe,R) 
k03= 0.03 #0.0162 
mL3= -4.2 #-2.888217523
  


# Nucleation Calculation data


r_hom = 5.103e-06 #4.78571e-06 #4.85277e-06 #3.9821e-06  #4.84982e-06              #4.86322e-06              # [m]
r_het = 1.84819e-06 #1.92353e-06 #1.84446e-06 #1.93141e-06 #1.87165e-06              #1.86834e-06              # [m]
thet = 0.745625  #0.810377 #0.773305 #0.952871 # 0.783201                     #0.780197                 # [rad]
dTdr = 3294 #4501.95 #3708.2 #7902.38 #4267.76                   #4205                     # [K.m^-1]
DT_= 0.420615 #0.541478 #0.452259 #0.79083 #0.52019                     #0.51692                  # [K]
sigma = -1.96801 #-2.13273 #-1.84011 #-2.09008 #-1.66031                 #-1.65425                 # [N.m^-1]
surface_stress = -0.591909 #-0.631105 #-1.00008 #-0.570309       #-0.564251                # [N.m^-1]
gamsl = 1.04701 #1.26959 #0.915939 #1.45814 #0.879113                 #0.871318                 # [J.m^-2]
dgamsldr =  -720831 #-789468 #-756628 #-615670 #-777538 #-576854               #-574516                  # [J.m^-3]
Gibbs_Thomson_hom = 3.88702e-07 #5.20781e-07 #4.17094e-07 #7.53446e-07 #4.86812e-07  #4.80167e-07              # [m.K]
Gibbs_Thomson_het_1st = 3.887e-07 #5.20776e-07  #4.17091e-07 #7.63702e-07 #4.86809e-07      #4.80e-07         # [m.K]
thet_2nd = 3.04729 #3.04634 #3.04734 #3.04624 #3.04699               #3.04703                  # [rad]
Gibbs_Thomson_het = 1.06852e-06 #1.29567e-06 #1.09735e-06 #1.57452e-06 #1.26141e-06 #1.25023e-06 ### 2nd-Order ### 3.41374e-07 # [m.K]      # 2nd Homogeneous --> 2.208e-07 
Gibbs_Thomson_non_Eq = Gibbs_Thomson_het
Gibbs_Thomson = Gibbs_Thomson_het
GIBBSTHOMSON = Gibbs_Thomson_het 




# start guess on the solvent concentration in solid alpha phase

csv = 100 - k01 * c01 - k02 * c02 - k03 * c03

# Experimental solidification kinetics
P  = np.arange(1,251, 1)
#VL = 1.2*P**(-0.27)     # VL Horizontal
VL = 3.4*P**(-0.56) #2.4*P**(-0.305)     # VL Horizontal
#V  = VL / 1000.0
#tSL= 3.58 * P ** (1.09)       # tSL Horizontal 
#tSL= 1.60 * P ** (1.00)
tSL= 4.93116 * P ** (0.87606)


# Experimental solidification kinetics
P_non_Eq  = np.arange(1,101, 1.0)
VL_non_Eq = 5.0*P**(-0.43)
#V  = VL / 1000.0
tSL_non_Eq= 0.32*P**(1.36)



def diff_length_scale(P, VL, tSL):
   V  = VL / 1000.0
   a = min(V)
   b = max(V)
   # Model Equations
   # Liquid diffusion length scale
   delta1 = DL1/V
   delta2 = DL2/V
   delta3 = DL3/V
   # Mean liquid boundary layer affected diffusion
   a1 = trapz(V,delta1)
   a2 = trapz(V,delta2)
   a3 = trapz(V,delta3)

   minV = min(V)
   maxV = max(V)

   delta1m=abs((1.0/(maxV-minV))*a1)
   delta2m=abs((1.0/(maxV-minV))*a2)
   delta3m=abs((1.0/(maxV-minV))*a3)


   # Effective partition coefficient
   kef1= k01/(k01+(1.0 - k01)*(np.exp((-V*delta1m/DL1))))
   kef2= k02/(k02+(1.0 - k02)*(np.exp((-V*delta2m/DL2))))
   kef3= k03/(k03+(1.0 - k03)*(np.exp((-V*delta3m/DL3))))


   # Solute solid fractions
   w1=(c01/100)/((c01/100)+(c02/100)+(c03/100))
   w2=(c02/100)/((c01/100)+(c02/100)+(c03/100))
   w3=(c03/100)/((c01/100)+(c02/100)+(c03/100))


   return delta1m, delta2m, delta3m, kef1, kef2, kef3, w1, w2, w3


def func_f(f, c01, k01, beta1, c02, k02, beta2, c03, k03, beta3, csv):
    csv = -c01 * k01 * (1 - f * (1 - beta1 * k01)) ** ((k01-1)/(1-beta1 * k01)) - c02 * k02 * (1 - f * (1 - beta2 * k02)) ** ((k02-1)/(1-beta2 * k02)) - c03 * k03 * (1 - f * (1 - beta3 * k03)) ** ((k03-1)/(1-beta3 * k03)) + 100
#    print(f'function_f({f}) = {function_f}')
    return csv

def cs(f, c0, k0, beta):
    cs = c0 * k0 * (1 - f * (1 - beta * k0)) ** ((k0-1)/(1-beta * k0))
    return cs

def T(f, TF, TL, k0, beta):
    T = TF + (TL-TF)* (1 - f * (1 - beta * k0)) ** ((k0-1)/(1-beta * k0))
    return T


def dfdtau(f, Fo, gama, beta):
    df_dtau = Fo * (gama - beta)/(f * gama * beta)
    return df_dtau

#def func_DH(f, rhos, cps, TL, Tf, beta, k):
    funcDH = ( rhos * cps * (TF - TL) * (1 - k) * (1 - f * (1 - beta * k)) ** ((k - 1)/(1 - beta * k)) )/((1 - beta * k)*(1 - f*(1- beta * k)))
    return funcDH
 
def dcs_approx(f, c0, k0, beta):
    dcs_app = c0*k0*(1 - k0)*(-f*(-beta*k0 + 1) + 1)**((k0 - 1)/(-beta*k0 + 1))/((-f*(-beta*k0 + 1) + 1))
    return dcs_app

def dcs(f, c0, k0, beta, dbeta):
    dcs_nonap = c0*k0*(-f*(-k0*beta + 1) + 1)**((k0 - 1)/(-k0*beta + 1)) * (k0*(k0 - 1)*log(-f*(-k0*beta + 1) + 1)*dbeta/(-k0*beta + 1)**2 + (k0 - 1)*(f*k0*dbeta + k0*beta - 1)/((-f*(-k0*beta + 1) + 1)*(-k0*beta + 1)))
    return dcs_nonap

def dbetadf(f, dfdtau, Fo, gama, beta):
    dbeta = -Fo * gama**2 * dfdtau / ((Fo + gama * f * dfdtau)**2)
    return dbeta


SDAS2 = 1.0
beta_CU = 0.0
beta_SI = 0.0
Error=1.0
def calc(G_T, delta1m, delta2m, delta3m, kef1, kef2, kef3, w1, w2, w3):
  GIBBSTHOMSON = G_T 
  SDAS2 = 1.0
  Error=1.0
  while np.all(op.gt(Error, 1.e-10)):
    SDAS = SDAS2
    # Fourier number / SDAS initial guess
    Fo1=DS1*tSL/((SDAS/2)**2.0)
    Fo2=DS2*tSL/((SDAS/2)**2.0)
    Fo3=DS3*tSL/((SDAS/2)**2.0)
    # Macrosegregation  parameters Gama
    gama1=(4.0*Fo1*k01)/(4.0*Fo1*k01+1.0)
    gama2=(4.0*Fo2*k02)/(4.0*Fo2*k02+1.0)
    gama3=(4.0*Fo3*k03)/(4.0*Fo3*k03+1.0)
    # Segregation parameter betha
    beta1=(2.0*gama1*Fo1)/(2.0*Fo1+gama1)
    beta2=(2.0*gama2*Fo2)/(2.0*Fo2+gama2) 
    beta3=(2.0*gama3*Fo2)/(2.0*Fo3+gama3)  
    # Weighted betha
    bethaweighted=(beta1*w1)+(beta2*w2)+(beta3*w3)  
    # Rappaz and Boettinger Modified Model
    p1=(mL1*(1.0 - kef1)*(cf1-c01)/DL1)+(mL2*(1.0-kef2)*(cf2-c02)/DL2)+(mL3*(1.0-kef3)*(cf3-c03)/DL3)
    s2=(mL1*(1.0 - kef1)*( cf1)/ DL1)+(mL2*(1.0-kef2)*(cf2)/DL2)+(mL3*(1.0-kef3)*(cf3)/DL3)
    t3=(mL1*(1.0 - kef1)*( c01)/ DL1)+(mL2*(1.0-kef2)*(c02)/DL2)+(mL3*(1.0-kef3)*(c03)/DL3)
    # M e M^1/3
    mrbm=bethaweighted*((-GIBBSTHOMSON/p1)*np.log(s2/t3))
    mrbm2=5.5*(mrbm**(1/3))
    # Corrected SDAS value
    SDAS2= mrbm2*(tSL)**(1/3)
    Error = np.abs(SDAS-SDAS2)/SDAS
    
  return SDAS, beta1, beta2, beta3

  


def calc_RB(G_T, k01, k02, k03, c01, cf1, c02, cf2,c03,cf3, mL1, mL2, mL3, DL1, DL2, DL3, tSL):
  GIBBSTHOMSON = G_T
  # Rappaz and Boettinger
  p1=(mL1*(1.0 - k01)*(cf1-c01)/DL1)+(mL2*(1.0-k02)*(cf2-c02)/DL2)+(mL3*(1.0-k03)*(cf3-c03)/DL3)
  s2=(mL1*(1.0 - k01)*( cf1)/ DL1)+(mL2*(1.0-k02)*(cf2)/DL2)+(mL3*(1.0-k03)*(cf3)/DL3)
  t3=(mL1*(1.0 - k01)*( c01)/ DL1)+(mL2*(1.0-k02)*(c02)/DL2)+(mL3*(1.0-k03)*(c03)/DL3)
  # M e M^1/3
  mrbm=((-GIBBSTHOMSON/p1)*np.log(s2/t3))
  mrbm2=5.5*(mrbm**(1/3))
  # Calculate SDAS value
  SDAS= mrbm2*(tSL)**(1/3)
  return SDAS


SDAS_RB = calc_RB(Gibbs_Thomson_het, k01, k02, k03, c01, cf1, c02, cf2, c03, cf3, mL1, mL2, mL3, DL1, DL2,DL3, tSL)
print(f'Gibbs-Thomson = {Gibbs_Thomson}\n')
print(f'k01 = {k01}')
print(f'k02 = {k02}')
print(f'k03 = {k03}')
print(f'c01 = {c01}')
print(f'cf1 = {cf1}')
print(f'c02 = {c02}')
print(f'cf2 = {cf2}')
print(f'c03 = {c03}')
print(f'cf3 = {cf3}')
print(f'mL1 = {mL1}')
print(f'mL2 = {mL2}')
print(f'mL3 = {mL3}')
print(f'DL1 = {DL1}')
print(f'DL2 = {DL2}')
print(f'DL3 = {DL3}')
print(f'tSL = {tSL}')

print(f'SDAS_RB = {SDAS_RB} [K.m]')
 

delta1m_Eq, delta2m_Eq, delta3m_Eq, kef1_Eq, kef2_Eq, kef3_Eq, w1_Eq, w2_Eq, w3_Eq = diff_length_scale(P, VL, tSL)
print(f'delta1m_Eq = {delta1m_Eq}')
SDAS_Eq, beta_CU, beta_MG, beta_FE = calc(Gibbs_Thomson_het, delta1m_Eq, delta2m_Eq, delta3m_Eq, kef1_Eq, kef2_Eq, kef3_Eq, w1_Eq, w2_Eq, w3_Eq)
print(f'SDAS_Eq = {SDAS_Eq} [K.m]')

SDAS_het_1st_Eq, beta_het_1st_CU, beta_het_1st_SI, beta_het_1st_FE = calc(Gibbs_Thomson_het_1st, delta1m_Eq, delta2m_Eq, delta3m_Eq, kef1_Eq, kef2_Eq, kef3_Eq, w1_Eq, w2_Eq, w3_Eq)
print(f'SDAS_Eq = {SDAS_Eq} [K.m]')



delta1m_non_Eq, delta2m_non_Eq, delta3m_non_Eq, kef1_non_Eq, kef2_non_Eq, kef3_non_Eq, w1_non_Eq, w2_non_Eq, w3_non_Eq = diff_length_scale(P_non_Eq, VL_non_Eq, tSL_non_Eq)
print(f'delta1m_non_Eq = {delta1m_non_Eq}')
SDAS_non_Eq, beta_CU_non_Eq, beta_Mg_non_Eq, beta_FE_non_Eq = calc(Gibbs_Thomson_het, delta1m_non_Eq, delta2m_non_Eq, delta3m_non_Eq, kef1_non_Eq, kef2_non_Eq, kef3_non_Eq, w1_non_Eq, w2_non_Eq, w3_non_Eq)
print(f'SDAS_non_Eq = {SDAS_non_Eq} [K.m]')

print('### Este arquivo ####')

# plot data
lx_RB = np.log(tSL)
ly_RB = np.log(SDAS_RB*1e06)
z3_RB= np.polyfit(lx_RB,ly_RB,1)
z3n_RB= np.poly1d(z3_RB)


lx_Eq = np.log(tSL)
ly_Eq = np.log(SDAS_Eq*1e06)
z3_Eq= np.polyfit(lx_Eq,ly_Eq,1)
z3n_Eq= np.poly1d(z3_Eq)


lx_1st_Eq = np.log(tSL)
ly_1st_Eq = np.log(SDAS_het_1st_Eq*1e06)
z3_1st_Eq= np.polyfit(lx_1st_Eq,ly_1st_Eq,1)
z3n_1st_Eq= np.poly1d(z3_1st_Eq)


lx_non_Eq = np.log(tSL_non_Eq)
ly_non_Eq = np.log(SDAS_non_Eq*1e06)
z3_non_Eq= np.polyfit(lx_non_Eq,ly_non_Eq,1)
z3n_non_Eq= np.poly1d(z3_non_Eq)

fig1,ax1 = plt.subplots(1)
ax1.set_yscale('log')
ax1.set_xscale('log')
#ax1.grid(color='k', which='both', linestyle='-', linewidth=0.5)
ax1.set_title(r'SDAS Models Predictions against Experimental Scatter',{'color': 'black', 'fontsize': 14})
###ax1.plot(tSL,(SDAS_RB*1e06),'r>',linewidth = 1.5)
ax1.plot(tSL,(np.exp(z3n_RB[0]))*tSL**z3n_RB[1],'-k',linewidth = 1.5)
###ax1.plot(tSL,(SDAS_Eq*1e06),'g^',linewidth = 1.5)
ax1.plot(tSL,(np.exp(z3n_1st_Eq[0]))*tSL**z3n_1st_Eq[1],'--k',linewidth = 1.5)   # plotted before 2nd-order so the legend order matches (was swapped)
###ax1.plot(tSL,(SDAS_het_1st_Eq*1e06),'b^',linewidth = 1.5)
ax1.plot(tSL,(np.exp(z3n_Eq[0]))*tSL**z3n_Eq[1],'-.k',linewidth = 1.5)


ax1.plot(xsc, ysc, 'sk',linewidth = 1.5)
#ax1.errorbar(xsc, ysc, uplims=yscp,linestyle=ls)



ax1.errorbar(xsc, ysc, yerr=yscp, fmt='.k')
ax1.errorbar(xsc, ysc, yerr=yscm, fmt='.k')

ax1.set_ylabel(u'Secondary Dendrite Arm Spacing, $SDAS$ [$\mathrm{{\mu}m}$]', {'color': 'black','fontweight':'normal', 'fontsize': 17})
ax1.set_xlabel(r'Local solidification time, $\mathrm{t_{SL}}$ [s]', {'color': 'black','fontweight':'normal', 'fontsize': 17})
plt.text(15,380,'Alloy Al-3wt%Cu-5wt%Nb-0.1wt%Fe Horizontal', {'color': 'black', 'fontsize': 16})
###plt.text(1.5,210,r'$\mathit{SDAS_{RB}\,=\,%.2f\,t_{SL}^{%.3f} \,\, \mathrm{[\mu m]}}$' % (np.exp(z3n_RB[0]), z3n_RB[1]), {'color': 'red', 'fontsize': 12})
###plt.text(1.5,120,r'$\mathit{SDAS_{MBR}\,=\,%.2f\,t_{SL}^{%.3f} \,\, \mathrm{[\mu m]}}$' % (np.exp(z3n_1st_Eq[0]), z3n_1st_Eq[1]), {'color': 'green', 'fontsize': 12})
###plt.text(1.5,70,r'$\mathit{SDAS_{MRB}\,=\,%.2f\,t_{SL}^{%.3f} \,\, \mathrm{[\mu m]}}$' % (np.exp(z3n_Eq[0]), z3n_Eq[1]), {'color': 'blue', 'fontsize': 12})
#plt.text(1.5,2000,r'$\mathit{V_{L}\,=\,5.0\,{P}^{-0.41} \,\, \mathrm{[mm.s^{-1}]}}$', {'color': 'black', 'fontsize': 12})
#plt.text(1.5,1100,r'$\mathit{t_{SL}\,=\,0.0036\,{P}^{\, 2.32} \,\, \mathrm{[s]}}$', {'color': 'black', 'fontsize': 12})
#plt.text(1.5,500,r'$\mathit{\overline{\Gamma } _{het} ^{\,1st-Order} \,=\, %g \,\, \mathrm{[m.K]}}$' % (Gibbs_Thomson_het_1st), {'color': 'green', 'fontsize': 12})
#plt.text(1.5,250,r'$\mathit{\theta  \,=\, %g \,\, \mathrm{[rad]}}$' % (thet), {'color': 'black', 'fontsize': 12})
#plt.text(1.5,130,r'$\mathit{\overline{\Gamma } ^{2nd-Order}_{het} \,=\, %g \,\, \mathrm{[m.K]}}$' % (Gibbs_Thomson_het), {'color': 'blue', 'fontsize': 12, 'weight':'bold' })
#plt.text(1.5,60,r'$\mathit{\theta ^{2nd-Order}  \,=\, %g \,\, \mathrm{[rad]}}$' % (thet_2nd), {'color': 'blue', 'fontsize': 12})
#plt.text(20.0,12000,r'$\mathit{\overline{\nabla T}_{EXP} \,= %g \, \mathrm{[K.m^{-1}]}}$' % (dTdr), {'color': 'black', 'fontsize': 12})
#plt.text(20.0,7000,r'$\mathit{r_{C,hom}\,= %g \, \mathrm{[m]}}$' % (r_hom), {'color': 'black', 'fontsize': 12})
#plt.text(20.0,4000,r'$\mathit{r_{C,het}\,= %g \, [m]}$' % (r_het), {'color': 'black', 'fontsize': 12})
#plt.text(20.0,2000,r'$\mathit{\sigma _{SL}\,= %g \, \mathrm{[N.m^{-1}]}}$' % (sigma), {'color': 'black', 'fontsize': 12})
#plt.text(20.0,1100,r'$\mathit{s \, = %g \, [N.m^{-1}]}$' % (surface_stress), {'color': 'black', 'fontsize': 12})
#plt.text(20.0,500,r'$\mathit{\gamma _{SL}\,= %g \, \mathrm{[J.m^{-2}]}}$' % (gamsl), {'color': 'black', 'fontsize': 12})
#plt.text(20.0,250,r'$\mathit{\frac {\partial \gamma _{SL}} {\partial r} \,= %g \, \mathrm{[J.m^{-3}]}}$' % (dgamsldr), {'color': 'black', 'fontsize': 12})
#plt.text(150,12000,r'$\mathit{\sigma _{0}\,= %g \, \mathrm{[N.m^{-1}]}}$' % (-sigmasl), {'color': 'black', 'fontsize': 12})
#plt.text(150,7000,r'$\mathit{\gamma _{0}\,= %g \, \mathrm{[J.m^{-2}]}}$' % (gamma_0), {'color': 'black', 'fontsize': 12})
#plt.text(1.5,20000,r'$^*$Thermal gradient $\mathrm{\nabla T}$ for horizontal experiments were only measured at the symmetry line along the X axis.' , {'color': 'black','weight':'normal', 'fontsize': 10})



plt.setp(ax1.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax1.get_yticklabels(), fontsize=16, fontweight="normal")
plt.legend((r'RB SDAS Model (Rappaz and Boettinger, 1999), $\mathit{SDAS_{RB}\,=\,%.2f\,t_{SL}^{%.3f} \,\, \mathrm{[\mu m]}}$' % (np.exp(z3n_RB[0]), z3n_RB[1]),r'MRB Model (Ferreira, 2024) with $\Gamma ^{1nd-Order} _{Het}$ (Present work), $\mathit{SDAS_{MBR}\,=\,%.2f\,t_{SL}^{%.3f} \,\, \mathrm{[\mu m]}}$' % (np.exp(z3n_1st_Eq[0]), z3n_1st_Eq[1]),r'MRB Model (Ferreira, 2024) with $\Gamma ^{2nd-Order} _{Het}$ (Present work), $\mathit{SDAS_{MRB}\,=\,%.2f\,t_{SL}^{%.3f} \,\, \mathrm{[\mu m]}}$' % (np.exp(z3n_Eq[0]), z3n_Eq[1]), r'Experimental (Mendes et al, 2023)'),frameon=True,edgecolor='k',fontsize=10,shadow=True, loc= 'lower right')
###plt.legend((r'RB SDAS Model (Rappaz and Boettinger, 1999)','RB Model Curve Fit',r'MRB Model (Ferreira et al.,2018) with $\Gamma ^{1nd-Order} _{Het}$ (Present work)',r'MRB Model with $\Gamma ^{1nd-Order} _{Het}$ Curve Fit ',r'MRB Model (Ferreira et al.,2018) with $\Gamma ^{2nd-Order} _{Het}$ (Present work)',r'MRB Model with $\Gamma ^{2nd-Order} _{Het}$ Curve Fit ', r'Experimental (Costa et al., 2015)'),frameon=True,edgecolor='k',fontsize=10,shadow=True, loc= 'lower right')
#plt.legend(('SDAS Predictions Equilibrium, Rappaz and Boettinger (1999)','SDAS Fit, Rappaz and Boettinger (1999)', r'MRB, Ferreira et al. (2018) and Ferreira (2022)' ,r'SDAS Fit , Ferreira et al. (2018) and Ferreira (2022)','Experimental'),frameon=True,edgecolor='k',shadow=True, loc= 'lower right')
#plt.axis([0,0.35,-0.2e7,0])
plt.axis([1e1,1e3,1e0,1e3])


fig2,ax2 = plt.subplots(1)
ax2.set_xscale('log')
#ax2.grid(color='k', which='both', linestyle='-', linewidth=0.5)
ax2.set_title(r'Horizontal Solidification Back-Diffusion Parameter, $\beta _i$',{'color': 'black', 'fontsize': 14})
ax2.plot(tSL,beta_CU,'k-',linewidth = 1.5)
ax2.plot(tSL,beta_MG,'k-.',linewidth = 1.5)
ax2.plot(tSL,beta_FE,'k--',linewidth = 1.5)
plt.text(5,0.175,'Alloy Al-3wt%Cu-5wt%Nb-0.1wt%Fe Horizontal', {'color': 'black', 'fontsize': 16})
#plt.text(5,0.150,r'$\mathit{T_{L}\,= %g \, [K]}$' % (TL), {'color': 'black', 'fontsize': 14})
#plt.text(5,0.137,r'$\mathit{T_{S}\,= %g \, [K]}$' % (Ts), {'color': 'black', 'fontsize': 14})
#plt.text(5,.125,r'$\mathit{D _{L,Si}\, = %.3e  e^{\frac {- %g} {R T}} \, \mathrm{[m^{2}.s^{-1}]}}$' % (Dol_Si,Qol_Si), {'color': 'black', 'fontsize': 14})
#plt.text(5,.1125,r'$\mathit{D _{S,Si}\, = %.3e  e^{\frac {- %g} {R T}} \, \mathrm{[m^{2}.s^{-1}]}}$' % (Dos_Si,Qos_Si), {'color': 'black', 'fontsize': 14})
#plt.text(5,0.100,r'$\mathit{D _{L,Cu}\, = %.3e  e^{\frac {- %g} {R T}} \, \mathrm{[m^{2}.s^{-1}]}}$' % (Dol_Cu,Qol_Cu), {'color': 'black', 'fontsize': 14})
#plt.text(5,0.087 ,r'$\mathit{D _{S,Cu}\, = %.3e  e^{\frac {- %g} {R T}} \, \mathrm{[m^{2}.s^{-1}]}}$' % (Dos_Cu,Qos_Cu), {'color': 'black', 'fontsize': 14})
#plt.text(5,0.075,r'$\mathit{D _{L,Fe}\, = %.3e  e^{\frac {- %g} {R T}} \, \mathrm{[m^{2}.s^{-1}]}}$' % (Dol_Fe,Qol_Fe), {'color': 'black', 'fontsize': 14})
#plt.text(5,0.063 ,r'$\mathit{D _{S,Fe}\, = %.3e  e^{\frac {- %g} {R T}} \, \mathrm{[m^{2}.s^{-1}]}}$' % (Dos_Fe,Qos_Fe), {'color': 'black', 'fontsize': 14})


ax2.set_xlabel(r'Local solidification time, $t_{SL}$ [s]', {'color': 'black','fontweight':'normal', 'fontsize': 17})
ax2.set_ylabel(r'Back-Diffusion Parameter, $\beta _i$  ', {'color': 'black','fontweight':'normal', 'fontsize': 17})
plt.setp(ax2.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax2.get_yticklabels(), fontsize=16, fontweight="normal")
plt.legend((r'$\beta _{Cu}$',r'$\beta _{Nb}$', r'$\beta _{Fe}$'),frameon=True,edgecolor='k',shadow=True, loc= 'upper right',fontsize= 17)
plt.axis([1,1e3,0,0.2])



###
### SDAS written explicitly as a function of the thermal gradient G_L
###
### Casting laws (Thermal_data_Ivaldo_B_MRS_Meeting_2026_Al3Cu5Nb_GL_data_indev.py):
###     G_L(P) = A + B P + C P^2,  P [mm], monotonic for P < -B/(2C) = 95.6 mm
###     tSL(P) = 4.93116 P^0.87606
### Inverting the quadratic:
###     P(G) = [ -B - sqrt(B^2 - 4 C (A - G)) ] / (2 C),   G_min = A - B^2/(4C) <= G <= A
###     tSL(G) = 4.93116 P(G)^0.87606
### and the MRB model  SDAS = 5.5 [ M(beta) tSL ]^(1/3)  becomes
###     SDAS(G) = 5.5 [ M(beta) 4.93116 P(G)^0.87606 ]^(1/3)
### or, with the fitted MRB power law SDAS = a tSL^b,
###     SDAS(G) = a 4.93116^b P(G)^(0.87606 b)
###
A_G, B_G, C_G = 10405.54, -200.22, 1.04692
P_vertex = -B_G / (2 * C_G)
G_min = A_G - B_G ** 2 / (4 * C_G)

def P_from_G(G):
    return (-B_G - np.sqrt(B_G ** 2 - 4 * C_G * (A_G - G))) / (2 * C_G)

def tSL_from_G(G):
    return 4.93116 * P_from_G(G) ** 0.87606

def P_from_tSL(t):
    return (t / 4.93116) ** (1 / 0.87606)

def G_from_P(P_):
    return A_G + B_G * P_ + C_G * P_ ** 2

print('\n===== SDAS explicit in the thermal gradient =====')
print(f'G_L(P) invertible for P < {P_vertex:.1f} mm  ->  {G_min:.1f} <= G_L <= {A_G:.1f} K/m')
print(f'P(G)  = [{-B_G:g} - sqrt({B_G**2:.2f} - {4*C_G:.5f} ({A_G:g} - G))] / {2*C_G:.5f}   [mm]')
print(f'tSL(G) = 4.93116 P(G)^0.87606   [s]')

# experimental SDAS inside the range where G_L(P) is measured and invertible (P <= 90 mm)
P_sd = P_from_tSL(xsc[0:16])
in_G = P_sd <= 90.0
G_sd = G_from_P(P_sd[in_G])
y_sd = ysc[0:16][in_G]
p_exp, lnC_exp = np.polyfit(np.log(G_sd), np.log(y_sd), 1)
print(f'Experimental, {in_G.sum()} points with P <= 90 mm:  SDAS = {np.exp(lnC_exp):.1f} G^{p_exp:.4f} [um]')

###
### Calibration of the MRB prefactor (5.5 in calc, Rappaz-Boettinger with classical Gamma)
### for the 2nd-order Gibbs-Thomson: same Gamma, same exponent, beta iterated consistently.
###
def calc_prefactor(G_T, C_MRB, delta1m, delta2m, delta3m, kef1, kef2, kef3, w1, w2, w3):
  """calc() with the constant 5.5 replaced by C_MRB."""
  GIBBSTHOMSON = G_T
  SDAS2 = 1.0
  Error = 1.0
  while np.all(op.gt(Error, 1.e-10)):
    SDAS = SDAS2
    Fo1=DS1*tSL/((SDAS/2)**2.0)
    Fo2=DS2*tSL/((SDAS/2)**2.0)
    Fo3=DS3*tSL/((SDAS/2)**2.0)
    gama1=(4.0*Fo1*k01)/(4.0*Fo1*k01+1.0)
    gama2=(4.0*Fo2*k02)/(4.0*Fo2*k02+1.0)
    gama3=(4.0*Fo3*k03)/(4.0*Fo3*k03+1.0)
    beta1=(2.0*gama1*Fo1)/(2.0*Fo1+gama1)
    beta2=(2.0*gama2*Fo2)/(2.0*Fo2+gama2)
    beta3=(2.0*gama3*Fo2)/(2.0*Fo3+gama3)
    bethaweighted=(beta1*w1)+(beta2*w2)+(beta3*w3)
    p1=(mL1*(1.0 - kef1)*(cf1-c01)/DL1)+(mL2*(1.0-kef2)*(cf2-c02)/DL2)+(mL3*(1.0-kef3)*(cf3-c03)/DL3)
    s2=(mL1*(1.0 - kef1)*( cf1)/ DL1)+(mL2*(1.0-kef2)*(cf2)/DL2)+(mL3*(1.0-kef3)*(cf3)/DL3)
    t3=(mL1*(1.0 - kef1)*( c01)/ DL1)+(mL2*(1.0-kef2)*(c02)/DL2)+(mL3*(1.0-kef3)*(c03)/DL3)
    mrbm=bethaweighted*((-GIBBSTHOMSON/p1)*np.log(s2/t3))
    SDAS2= C_MRB*(mrbm**(1/3))*(tSL)**(1/3)
    Error = np.abs(SDAS-SDAS2)/SDAS
  return SDAS

args_Eq = (delta1m_Eq, delta2m_Eq, delta3m_Eq, kef1_Eq, kef2_Eq, kef3_Eq, w1_Eq, w2_Eq, w3_Eq)
def mean_log_residual(C_MRB):
    S_ = calc_prefactor(Gibbs_Thomson_het, C_MRB, *args_Eq) * 1e6
    return np.mean(np.log(np.interp(xsc[0:16][in_G], tSL, S_)) - np.log(y_sd))

assert abs(np.max(np.abs(calc_prefactor(Gibbs_Thomson_het, 5.5, *args_Eq) - SDAS_Eq))) < 1e-12   # reproduces calc()
C_MRB_2nd = brentq(mean_log_residual, 1.0, 10.0, xtol=1e-12)
SDAS_Eq_adj = calc_prefactor(Gibbs_Thomson_het, C_MRB_2nd, *args_Eq)
z3n_Eq_adj = np.poly1d(np.polyfit(np.log(tSL), np.log(SDAS_Eq_adj*1e06), 1))
print(f'\nMRB prefactor calibrated for Gamma_het 2nd-order: C = {C_MRB_2nd:.4f} (was 5.5, ratio {C_MRB_2nd/5.5:.4f})')

G_grid = np.linspace(G_from_P(90.0), G_from_P(5.0), 300)
t_grid = tSL_from_G(G_grid)
SDAS_vs_G = []
cases_G = [
    (r'MRB, $\Gamma ^{1st-Order} _{Het}$ = %.3e m.K' % (Gibbs_Thomson_het_1st), SDAS_het_1st_Eq, z3n_1st_Eq, '#2a78d6', '-'),
    (r'MRB, $\Gamma ^{2nd-Order} _{Het}$ = %.3e m.K' % (Gibbs_Thomson_het), SDAS_Eq, z3n_Eq, '#eb6834', '--'),
    (r'MRB, $\Gamma ^{2nd-Order} _{Het}$, prefactor %.2f (calibrated)' % (C_MRB_2nd), SDAS_Eq_adj, z3n_Eq_adj, '#1baf7a', '-.'),
]
for lab, SDAS_t, z_t, col, ls in cases_G:
    S_G = np.interp(t_grid, tSL, SDAS_t * 1e6)            # full MRB with beta(tSL)
    S_G_closed = np.exp(z_t[0]) * t_grid ** z_t[1]         # closed form from the MRB power law
    p_mod, lnC_mod = np.polyfit(np.log(G_grid), np.log(S_G), 1)
    S_at_exp = np.interp(G_sd, G_grid, S_G)
    mape = 100 * np.mean(np.abs(S_at_exp - y_sd) / y_sd)
    SDAS_vs_G.append((S_G, p_mod, lnC_mod))
    print(f'{lab.replace("$", "").replace(chr(92), ""):<45}  closed form: SDAS(G) = {np.exp(z_t[0]):.2f} 4.93116^{z_t[1]:.3f} P(G)^{0.87606*z_t[1]:.4f}'
          f'  | power fit: SDAS = {np.exp(lnC_mod):.1f} G^{p_mod:.4f}  | max|MRB - closed| = {np.max(np.abs(S_G - S_G_closed)):.2f} um  | MAPE = {mape:.1f}%')

fileG = open('SDAS_vs_GL_B_MRS_2026.csv', 'w')
fileG.write('GL;P_mm;tSL;SDAS_MRB_1st;SDAS_MRB_2nd;SDAS_MRB_2nd_calibrated\n')
for k in range(len(G_grid)):
    fileG.write(f'{G_grid[k]:.3f};{P_from_G(G_grid[k]):.4f};{t_grid[k]:.4f};{SDAS_vs_G[0][0][k]:.4f};{SDAS_vs_G[1][0][k]:.4f};{SDAS_vs_G[2][0][k]:.4f}\n')
fileG.close()

###
### Exponent -1/3:  tSL = DT_SL / (G_L V_L)  ->  SDAS = 5.5 [M DT_SL / (G_L V_L)]^(1/3)
###                 SDAS = K (G_L V_L)^(-1/3),  K = 5.5 (M DT_SL)^(1/3)
### V_L(G) = 3.4 P(G)^-0.56 [mm/s] (casting law), so SDAS is explicit in G_L.
### K is the only fitted constant (least squares in log with the slope fixed at -1/3).
###
def VL_from_G(G):
    return 3.4 * P_from_G(G) ** (-0.56) / 1000.0      # [m/s]

def K_one_third(G_, S_):
    """Prefactor of SDAS = K (G V_L)^(-1/3) by log least squares with fixed slope."""
    return np.exp(np.mean(np.log(S_) + np.log(G_ * VL_from_G(G_)) / 3.0))

Tdot_sd   = G_sd * VL_from_G(G_sd)
Tdot_grid = G_grid * VL_from_G(G_grid)
K_exp = K_one_third(G_sd, y_sd)
law_cases = [(r'Experimental fit', K_exp, 'k', ':')]
for (lab, _, _, col, ls), (S_G, _, _) in zip(cases_G, SDAS_vs_G):
    law_cases.append((lab, K_one_third(G_grid, S_G), col, ls))

def metrics(S_mod):
    rel = (S_mod - y_sd) / y_sd
    return 100 * np.mean(np.abs(rel)), np.sqrt(np.mean((S_mod - y_sd) ** 2))

print('\n===== Exponent -1/3 test:  SDAS = K (G_L V_L)^(-1/3) =====')
# (1) free slope against the cooling rate, with 95% confidence interval
X_ = np.log(Tdot_sd); Y_ = np.log(y_sd); n_ = len(X_)
slope_, icpt_ = np.polyfit(X_, Y_, 1)
res_ = Y_ - (slope_ * X_ + icpt_)
se_  = np.sqrt(np.sum(res_ ** 2) / (n_ - 2) / np.sum((X_ - X_.mean()) ** 2))
from scipy.stats import t as t_student
ci_ = t_student.ppf(0.975, n_ - 2) * se_
print(f'(1) free fit, experimental: SDAS = {np.exp(icpt_):.2f} Tdot^({slope_:.3f} +/- {ci_:.3f})  (95% CI, n = {n_})'
      f'  -> -1/3 {"inside" if abs(slope_ + 1/3) <= ci_ else "OUTSIDE"} the interval')
for (lab, _, _, _, _), (S_G, _, _) in zip(cases_G, SDAS_vs_G):
    s_mod = np.polyfit(np.log(Tdot_grid), np.log(S_G), 1)[0]
    print(f'    free fit, {lab.replace("$", "").replace(chr(92), ""):<42}: slope against Tdot = {s_mod:.3f}')
# (2) accuracy of each -1/3 law against the experimental points
print('(2) accuracy against the experimental points:')
for lab, K_, _, _ in law_cases:
    mape_, rmse_ = metrics(K_ * Tdot_sd ** (-1/3))
    print(f'    {lab.replace("$", "").replace(chr(92), ""):<44} K = {K_:7.3f} um.(K/s)^(1/3)   MAPE = {mape_:5.1f}%   RMSE = {rmse_:5.2f} um')
for (lab, _, _, _, _), (S_G, _, _) in zip(cases_G, SDAS_vs_G):
    mape_, rmse_ = metrics(np.interp(G_sd, G_grid, S_G))
    print(f'    {"full MRB (tSL law), " + lab.replace("$", "").replace(chr(92), ""):<44}                         MAPE = {mape_:5.1f}%   RMSE = {rmse_:5.2f} um')

fig5,ax5 = plt.subplots(1)
ax5.set_xscale('log')
ax5.set_yscale('log')
ax5.set_title(r'SDAS = K $(\nabla T_{L}\,V_{L})^{-1/3}$ as an explicit function of $\nabla T_{L}$',{'color': 'black', 'fontsize': 14})
for lab, K_, col, ls in law_cases:
    mape_, _ = metrics(K_ * Tdot_sd ** (-1/3))
    ax5.plot(G_grid, K_ * Tdot_grid ** (-1/3), color=col, linestyle=ls, linewidth=2.0,
             label=lab + r', K = %.2f, MAPE = %.1f%%' % (K_, mape_))
ax5.plot(G_sd, y_sd, 'sk', markersize=6, label=r'Experimental (Mendes et al, 2023)')
ax5.set_xlabel(r'Thermal Gradient, $\mathrm{\nabla T_{L} \,\, [K.m^{-1}] }$', {'color': 'black','fontweight':'normal', 'fontsize': 17})
ax5.set_ylabel(u'Secondary Dendrite Arm Spacing, $SDAS$ [$\mathrm{{\mu}m}$]', {'color': 'black','fontweight':'normal', 'fontsize': 17})
plt.setp(ax5.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax5.get_yticklabels(), fontsize=16, fontweight="normal")
ax5.grid(True, which='both', linestyle='--', color='darkgray', linewidth=0.5)
ax5.legend(frameon=True, edgecolor='k', fontsize=10, shadow=True, loc='upper center', bbox_to_anchor=(0.5, -0.2))
fig5.savefig('Figure_SDAS_vs_GL_B_MRS_2026.png', dpi=200, bbox_inches='tight')

# full MRB curves against G_L with their free (apparent) exponents
fig7,ax7 = plt.subplots(1)
ax7.set_xscale('log')
ax7.set_yscale('log')
ax7.set_title(r'SDAS as an explicit function of the thermal gradient, $\nabla T_{L}$',{'color': 'black', 'fontsize': 14})
for (lab, _, _, col, ls), (S_G, p_mod, lnC_mod) in zip(cases_G, SDAS_vs_G):
    mape_, _ = metrics(np.interp(G_sd, G_grid, S_G))
    ax7.plot(G_grid, S_G, color=col, linestyle=ls, linewidth=2.0,
             label=lab + r', $SDAS \approx %.0f\,\nabla T_{L}^{%.3f}$, MAPE = %.1f%%' % (np.exp(lnC_mod), p_mod, mape_))
ax7.plot(G_sd, y_sd, 'sk', markersize=6, label=r'Experimental (Mendes et al, 2023), $SDAS = %.0f\,\nabla T_{L}^{%.3f}$' % (np.exp(lnC_exp), p_exp))
ax7.plot(G_grid, np.exp(lnC_exp) * G_grid ** p_exp, ':k', linewidth=1.5)
ax7.set_xlabel(r'Thermal Gradient, $\mathrm{\nabla T_{L} \,\, [K.m^{-1}] }$', {'color': 'black','fontweight':'normal', 'fontsize': 17})
ax7.set_ylabel(u'Secondary Dendrite Arm Spacing, $SDAS$ [$\mathrm{{\mu}m}$]', {'color': 'black','fontweight':'normal', 'fontsize': 17})
plt.setp(ax7.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax7.get_yticklabels(), fontsize=16, fontweight="normal")
ax7.grid(True, which='both', linestyle='--', color='darkgray', linewidth=0.5)
ax7.legend(frameon=True, edgecolor='k', fontsize=9, shadow=True, loc='upper center', bbox_to_anchor=(0.5, -0.2))
fig7.savefig('Figure_SDAS_vs_GL_free_exponents_B_MRS_2026.png', dpi=200, bbox_inches='tight')

# test figure: against the cooling rate the -1/3 law is a straight line
fig6,ax6 = plt.subplots(1)
ax6.set_xscale('log')
ax6.set_yscale('log')
ax6.set_title(r'Test of the exponent $-1/3$ against the cooling rate $\dot T = \nabla T_{L}\,V_{L}$',{'color': 'black', 'fontsize': 14})
for lab, K_, col, ls in law_cases:
    ax6.plot(Tdot_grid, K_ * Tdot_grid ** (-1/3), color=col, linestyle=ls, linewidth=2.0, label=lab + r', slope $-1/3$')
ax6.plot(Tdot_grid, np.exp(icpt_) * Tdot_grid ** slope_, color='gray', linestyle='-', linewidth=1.0,
         label=r'Free fit, slope $%.3f \pm %.3f$ (95%%)' % (slope_, ci_))
ax6.plot(Tdot_sd, y_sd, 'sk', markersize=6, label=r'Experimental (Mendes et al, 2023)')
ax6.set_xlabel(r'Cooling rate, $\mathrm{\dot T = \nabla T_{L}\,V_{L} \,\, [K.s^{-1}] }$', {'color': 'black','fontweight':'normal', 'fontsize': 17})
ax6.set_ylabel(u'Secondary Dendrite Arm Spacing, $SDAS$ [$\mathrm{{\mu}m}$]', {'color': 'black','fontweight':'normal', 'fontsize': 17})
plt.setp(ax6.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax6.get_yticklabels(), fontsize=16, fontweight="normal")
ax6.grid(True, which='both', linestyle='--', color='darkgray', linewidth=0.5)
ax6.legend(frameon=True, edgecolor='k', fontsize=10, shadow=True, loc='upper center', bbox_to_anchor=(0.5, -0.2))
fig6.savefig('Figure_SDAS_vs_cooling_rate_one_third_test_B_MRS_2026.png', dpi=200, bbox_inches='tight')
fig1.savefig('Figure_SDAS_vs_tSL_legend_fixed_B_MRS_2026.png', dpi=200, bbox_inches='tight')


plt.show()

