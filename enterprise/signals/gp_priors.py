# gp_priors.py
"""Utilities module containing various useful
functions for use in other modules. Adapted for my Master's project.
"""

import numpy as np
import scipy.stats
import scipy.integrate as it

from enterprise.signals import parameter
from enterprise.signals.parameter import function
import enterprise.constants as const
import enterprise.signals.DOF_Adjustments as DOF

@function
def f_rh(log10_T_rh=9):
    return 3 * 10**(log10_T_rh) / (10**17)

@function
def Tf(f,log10_T_rh=9):
    return (
        const.T_0*f**2*const.eta_0**2/4*np.heaviside(f - 2/const.eta_0,1)*np.heaviside(const.f_eq - f,1) + 1/const.f_eq*const.T_eq*f*np.heaviside(f - const.f_eq,1)*np.heaviside(f_rh(log10_T_rh) - f,1) + (f/f_rh(log10_T_rh))**2 * 10**(log10_T_rh) * np.heaviside(f - f_rh(log10_T_rh),1) 
    )

@function
def fT(T,log10_T_rh=9):
    return (
        2*np.sqrt(T/const.T_0)/const.eta_0*np.heaviside(T-const.T_0,1)*np.heaviside(const.T_eq - T,1) + T/(1/const.f_eq * const.T_eq)*np.heaviside(T - const.T_eq,1)*np.heaviside(10**(log10_T_rh) - T,1) + np.sqrt(T/10**(log10_T_rh))*f_rh(log10_T_rh)*np.heaviside(T - 10**(log10_T_rh),1)
    )


@function
def Transfer_function(f, log10_T_rh=9, log10_f_inf=10):
    return(
       (DOF.return_DOFge(Tf(f,log10_T_rh))/DOF.return_DOFge(0))*(DOF.return_DOFgs(0)/DOF.return_DOFgs(Tf(f,log10_T_rh)))**(4/3) * const.Om_Mat**2 * 9 / (2 * np.pi * f * const.eta_0)**4 * (1 + 1.57*(f/const.f_eq) + 3.42*(f/const.f_eq)**2) * 1 / (1 - 0.22*(f/f_rh(log10_T_rh=log10_T_rh))**1.5 + 0.65*(f/f_rh(log10_T_rh=log10_T_rh))**2) * np.heaviside(log10_f_inf - np.log10(f),1) 
    )


@function
def Power_Spectrum(f, log10_r=-1.6, n_t=6, components=2):
    return(
        10**(float(log10_r))*const.A_s*(f/const.f_ref)**n_t
    )


@function
def custom_powerlaw(f, log10_r=-1.6, n_t=6, log10_T_rh=9, log10_f_inf=10, components=2):
    df = np.diff(np.concatenate((np.array([0]), f[::components])))
    #BBN prior
    if const.hc**2 * it.quad(lambda x: Power_Spectrum(x, log10_r=log10_r, n_t= n_t) * Transfer_function(x, log10_T_rh=log10_T_rh, log10_f_inf=log10_f_inf) * 1 / 12 *(2*np.pi*x)**2 / (const.H_0*1000)**2 /x, const.f_BBN, 10**(log10_f_inf))[0] > 5.6*10**(-6)*const.DelN:
        p = np.inf
    #LVK_prior
    elif Power_Spectrum(const.f_LVK, log10_r=log10_r, n_t= n_t) * Transfer_function(const.f_LVK, log10_T_rh=log10_T_rh, log10_f_inf=log10_f_inf) * 1 / 12 *(2*np.pi*const.f_LVK)**2 / (const.H_0*1000)**2 > const.Om_LVK:
        p = np.inf
    #f_inf_prior
    elif 10**(log10_f_inf) < fT(5*10**(6), log10_T_rh=log10_T_rh):
        p = np.inf
    else:
        p = 1 / 24 / np.pi**2 / f**3 *Power_Spectrum(f, log10_r=log10_r, n_t=n_t)*Transfer_function(f, log10_T_rh=log10_T_rh, log10_f_inf=log10_f_inf) 
    return(
        p * np.repeat(df, components)
    )

#BBN prior function

@function
def BBN_prior(f,log10_r=-1.6, n_t=6, log10_T_rh=9, log10_f_inf=10, components=2):
    df = np.diff(np.concatenate((np.array([0]), f[::components])))
    if const.hc**2 * it.quad(lambda x: Power_Spectrum(x, log10_r=log10_r, n_t= n_t) * Transfer_function(x, log10_T_rh=log10_T_rh, log10_f_inf=log10_f_inf) * 1 / 12 *(2*np.pi*x)**2 / (const.H_0*1000)**2 /x, const.f_BBN, 10**(log10_f_inf))[0] > 5.6*10**(-6)*const.DelN:
        #1 is much larger than the standard output
        p = np.ones(len(df))
    else:
        #0 does not change the standard output
        p = 0*df
    return(
        np.repeat(p,components)
    )

#LVK prior function

@function
def LVK_prior(f, log10_r=-1.6, n_t=6, log10_T_rh=9, log10_f_inf=10, components=2):
    df = np.diff(np.concatenate((np.array([0]),f[::components])))
    if Power_Spectrum(const.f_LVK, log10_r=log10_r, n_t= n_t) * Transfer_function(const.f_LVK, log10_T_rh=log10_T_rh, log10_f_inf=log10_f_inf) * 1 / 12 *(2*np.pi*const.f_LVK)**2 / (const.H_0*1000)**2 > const.Om_LVK:
        #1 is much larger than the standard output
        p = np.ones(len(df))
    else:
        #0 does not change the standard output
        p = 0*df
    return(
        np.repeat(p,components)
    )
#Prior on f_inf (lower bound at T = 5 MeV)

@function
def f_inf_prior(f, log10_T_rh=9, log10_f_inf=10, components=2):
    df = np.diff(np.concatenate((np.array([0]),f[::components])))
    if 10**(log10_f_inf) < fT(5*10**(6), log10_T_rh=log10_T_rh):
        p = np.ones(len(df))
    else:
        p = 0*df
    return(
        np.repeat(p,components)
    )
    
@function
def powerlaw(f, log10_A=-16, gamma=5, components=2):
    df = np.diff(np.concatenate((np.array([0]), f[::components])))
    return (
        (float(10**(float(log10_A))) ** 2 / 12.0 / np.pi**2 * const.fyr ** (gamma - 3) * f ** (-gamma) * np.repeat(df, components))
    )
    

@function
def turnover(f, log10_A=-15, gamma=4.33, lf0=-8.5, kappa=10 / 3, beta=0.5):
    df = np.diff(np.concatenate((np.array([0]), f[::2])))
    hcf = 10**(float(log10_A)) * (f / const.fyr) ** ((3 - gamma) / 2) / (1 + (10**lf0 / f) ** kappa) ** beta
    return hcf**2 / 12 / np.pi**2 / f**3 * np.repeat(df, 2)


@function
def free_spectrum(f, log10_rho=None):
    """
    Free spectral model. PSD  amplitude at each frequency
    is a free parameter. Model is parameterized by
    S(f_i) = \rho_i^2 * T,
    where \rho_i is the free parameter and T is the observation
    length.
    """
    return np.repeat(10 ** (2 * np.array(log10_rho)), 2)


@function
def t_process(f, log10_A=-15, gamma=4.33, alphas=None):
    """
    t-process model. PSD  amplitude at each frequency
    is a fuzzy power-law.
    """
    alphas = np.ones_like(f) if alphas is None else np.repeat(alphas, 2)
    return powerlaw(f, log10_A=log10_A, gamma=gamma) * alphas


@function
def t_process_adapt(f, log10_A=-15, gamma=4.33, alphas_adapt=None, nfreq=None):
    """
    t-process model. PSD  amplitude at each frequency
    is a fuzzy power-law.
    """
    if alphas_adapt is None:
        alpha_model = np.ones_like(f)
    else:
        if nfreq is None:
            alpha_model = np.repeat(alphas_adapt, 2)
        else:
            alpha_model = np.ones_like(f)
            alpha_model[2 * int(np.rint(nfreq))] = alphas_adapt
            alpha_model[2 * int(np.rint(nfreq)) + 1] = alphas_adapt

    return powerlaw(f, log10_A=log10_A, gamma=gamma) * alpha_model


def InvGammaPrior(value, alpha=1, gamma=1):
    """Prior function for InvGamma parameters."""
    return scipy.stats.invgamma.pdf(value, alpha, scale=gamma)


def InvGammaSampler(alpha=1, gamma=1, size=None):
    """Sampling function for Uniform parameters."""
    return scipy.stats.invgamma.rvs(alpha, scale=gamma, size=size)


def InvGamma(alpha=1, gamma=1, size=None):
    """Class factory for Inverse Gamma parameters."""

    class InvGamma(parameter.Parameter):
        _size = size
        _prior = parameter.Function(InvGammaPrior, alpha=alpha, gamma=gamma)
        _sampler = staticmethod(InvGammaSampler)
        _alpha = alpha
        _gamma = gamma

        def __repr__(self):
            return '"{}": InvGamma({},{})'.format(self.name, alpha, gamma) + (
                "" if self._size is None else "[{}]".format(self._size)
            )

    return InvGamma


@function
def turnover_knee(f, log10_A, gamma, lfb, lfk, kappa, delta):
    """
    Generic turnover spectrum with a high-frequency knee.
    :param f: sampling frequencies of GWB
    :param A: characteristic strain amplitude at f=1/yr
    :param gamma: negative slope of PSD around f=1/yr (usually 13/3)
    :param lfb: log10 transition frequency at which environment dominates GWs
    :param lfk: log10 knee frequency due to population finiteness
    :param kappa: smoothness of turnover (10/3 for 3-body stellar scattering)
    :param delta: slope at higher frequencies
    """
    df = np.diff(np.concatenate((np.array([0]), f[::2])))
    hcf = (
        10**log10_A
        * (f / const.fyr) ** ((3 - gamma) / 2)
        * (1.0 + (f / 10**lfk)) ** delta
        / np.sqrt(1 + (10**lfb / f) ** kappa)
    )
    return hcf**2 / 12 / np.pi**2 / f**3 * np.repeat(df, 2)


@function
def broken_powerlaw(f, log10_A, gamma, delta, log10_fb, kappa=0.1):
    """
    Generic broken powerlaw spectrum.
    :param f: sampling frequencies
    :param A: characteristic strain amplitude [set for gamma at f=1/yr]
    :param gamma: negative slope of PSD for f > f_break [set for comparison
        at f=1/yr (default 13/3)]
    :param delta: slope for frequencies < f_break
    :param log10_fb: log10 transition frequency at which slope switches from
        gamma to delta
    :param kappa: smoothness of transition (Default = 0.1)
    """
    df = np.diff(np.concatenate((np.array([0]), f[::2])))
    hcf = (
        10**log10_A
        * (f / const.fyr) ** ((3 - gamma) / 2)
        * (1 + (f / 10**log10_fb) ** (1 / kappa)) ** (kappa * (gamma - delta) / 2)
    )
    return hcf**2 / 12 / np.pi**2 / f**3 * np.repeat(df, 2)


@function
def powerlaw_genmodes(f, log10_A=-16, gamma=5, components=2, wgts=None):
    if wgts is not None:
        df = wgts**2
    else:
        df = np.diff(np.concatenate((np.array([0]), f[::components])))
    return (
        (10**log10_A) ** 2 / 12.0 / np.pi**2 * const.fyr ** (gamma - 3) * f ** (-gamma) * np.repeat(df, components)
    )


@function
def infinitepower(f):
    return np.full_like(f, 1e40, dtype="d")
