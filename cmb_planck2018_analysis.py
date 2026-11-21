# !/ usr / bin / env python3
" " "
cmb_planck2018_analysis . py
Planck 2018 TT constraint on the consciousness - point power - spectrum
oscillation .
Data : https :// pla . esac . esa . int /# maps ( COM_PowerSpect_CMB - TT - full_R3 .01. txt )
" " "
import numpy as np
import emcee
import matplotlib . pyplot as plt
from scipy . interpolate import CubicSpline

def load_planck_tt ( path ) :
ell , Dl , err = np . loadtxt ( path , usecols =(0 ,1 ,2) , unpack = True )
mask = ( ell >= 30) & ( ell <= 2500)
return ell [ mask ] , Dl [ mask ] , err [ mask ]

def load_lcdm_template ( path = " lcdm_tt . txt " ) :
ell0 , Dl0 = np . loadtxt ( path , usecols =(0 ,1) , unpack = True )
return CubicSpline ( ell0 , Dl0 )

def model_Dl ( ell , A_T , phi , k_star , lcdm_spline ) :
Dl0 = lcdm_spline ( ell )
return Dl0 * (1.0 + A_T * np . sin (2* np . pi * ell / k_star + phi ) )

def log_likelihood ( theta , ell , Dl , err , lcdm_spline ) :
A_T , phi , k_star = theta
if not (0 < A_T < 0.5 and 0 <= phi < 2* np . pi and 50 < k_star < 2000) :
return - np . inf
mu = model_Dl ( ell , A_T , phi , k_star , lcdm_spline )
chisq = np . sum ((( Dl - mu ) / err ) **2)
return -0.5 * chisq

def run_mcmc ( data_file , lcdm_file , nwalkers =32 , nsteps =3000) :
ell , Dl , err = load_planck_tt ( data_file )
lcdm_spline = load_lcdm_template ( lcdm_file )
ndim = 3
p0 = np . array ([1 e -3 , 0.0 , 500.0]) + 1e -4* np . random . randn ( nwalkers , ndim )
sampler = emcee . EnsembleSampler ( nwalkers , ndim , log_likelihood ,
8
args =( ell , Dl , err , lcdm_spline ) )
sampler . run_mcmc ( p0 , nsteps , progress = True )
return sampler

def plot_constraints ( sampler , burn_in =1000) :
samples = sampler . get_chain ( discard = burn_in , flat = True )
A_T = samples [: ,0]
k_star = samples [: ,2]
fig , ax = plt . subplots ( figsize =(6 ,4) )
A_95 = np . percentile ( A_T , 95)
ax . hist2d ( k_star , A_T , bins =50 , cmap = ’ Blues ’)
ax . axhline ( A_95 , color = ’r ’ , ls = ’ -- ’ , label = f ’ 95% UL A_T ={ A_95 :.2 e } ’)
ax . set_xlabel ( r ’ Characteristic scale $ \ ell_ * $ ’)
ax . set_ylabel ( r ’ Amplitude $A_T$ ’)
 ax . set_title ( ’ Planck 2018 constraint on consciousness - point oscillation ’
)
ax . legend ()
plt . tight_layout ()
plt . savefig ( ’ cmb_constraint . png ’ , dpi =300)
print ( f " Planck 2018 95% UL on A_T : { A_95 :.2 e } " )

if __name__ == " __main__ " :
sampler = run_mcmc ( " COM_PowerSpect_CMB - TT - full_R3 .01. txt " , " lcdm_tt . txt "
)
plot_constraints ( sampler )
