#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_gp_coefficients
----------------------------------

Tests for GP signals used with deterministic coefficients.
"""

import logging
import unittest
import pytest

import numpy as np
import scipy.sparse as sps
from sksparse.cholmod import cholesky

from enterprise.pulsar import Pulsar
from enterprise.signals import (
    deterministic_signals,
    gp_signals,
    parameter,
    selections,
    signal_base,
    utils,
    white_signals,
)
from enterprise.signals.selections import Selection
from tests.enterprise_test_data import datadir
from tests.enterprise_test_data import LIBSTEMPO_INSTALLED, PINT_INSTALLED

logging.basicConfig(format="%(levelname)s: %(name)s: %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)


@signal_base.function
def create_quant_matrix(toas, dt=1):
    U, _ = utils.create_quantization_matrix(toas, dt=dt, nmin=1)
    avetoas = np.array([toas[idx.astype(bool)].mean() for idx in U.T])
    # return value slightly different than 1 to get around ECORR columns
    return U * 1.0000001, avetoas


@signal_base.function
def se_kernel(etoas, log10_sigma=-7, log10_lam=np.log10(30 * 86400)):
    tm = np.abs(etoas[None, :] - etoas[:, None])
    d = np.eye(tm.shape[0]) * 10 ** (2 * (log10_sigma - 1.5))
    return 10 ** (2 * log10_sigma) * np.exp(-(tm**2) / 2 / 10 ** (2 * log10_lam)) + d


class TestGPCoefficients(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Setup the Pulsar object."""

        # initialize Pulsar class
        cls.psr = Pulsar(datadir + "/B1855+09_NANOGrav_9yv1.t2.feather")
        cls.psr2 = Pulsar(datadir + "/B1937+21_NANOGrav_9yv1.t2.feather")

    def test_ephemeris(self):
        """Test physical-ephemeris delay, made three ways: from
        marginalized GP, from coefficient-based GP, from
        deterministic model."""

        ef = white_signals.MeasurementNoise(efac=parameter.Uniform(0.1, 5.0))

        eph = gp_signals.FourierBasisCommonGP_physicalephem(sat_orb_elements=None)

        ephc = gp_signals.FourierBasisCommonGP_physicalephem(sat_orb_elements=None, coefficients=True)

        ephd = deterministic_signals.PhysicalEphemerisSignal(sat_orb_elements=False)

        model = ef + eph
        modelc = ef + ephc
        modeld = ef + ephd

        pta = signal_base.PTA([model(self.psr), model(self.psr2)])
        ptac = signal_base.PTA([modelc(self.psr), modelc(self.psr2)])
        ptad = signal_base.PTA([modeld(self.psr), modeld(self.psr2)])

        cf = 1e-3 * np.random.randn(11)
        cf[0] = 1e-5  # this is more sensitive to linearity

        bs = pta.get_basis()
        da = [np.dot(bs[0], cf), np.dot(bs[1], cf)]

        params = {
            "B1855+09_efac": 1,
            "B1937+21_efac": 1,
            "B1855+09_phys_ephem_gp_coefficients": cf,
            "B1937+21_phys_ephem_gp_coefficients": cf,
        }
        db = ptac.get_delay(params=params)

        dparams = {
            "B1855+09_efac": 1,
            "B1937+21_efac": 1,
            "frame_drift_rate": cf[0],
            "d_jupiter_mass": cf[1],
            "d_saturn_mass": cf[2],
            "d_uranus_mass": cf[3],
            "d_neptune_mass": cf[4],
            "jup_orb_elements": cf[5:],
        }
        dc = ptad.get_delay(params=dparams)

        msg = "Reconstructed ephemeris signals differ!"

        assert np.allclose(da[0], db[0]), msg
        assert np.allclose(da[1], db[1]), msg

        # we don't expect an exact match since we are linearizing
        assert np.allclose(da[0], dc[0], atol=1e-3), msg
        assert np.allclose(da[1], dc[1], atol=1e-3), msg

    def test_common_red_noise(self):
        """Test of a coefficient-based common GP."""
        pl = utils.powerlaw(log10_A=parameter.Uniform(-18, -12), gamma=parameter.Uniform(1, 7))

        ef = white_signals.MeasurementNoise(efac=parameter.Uniform(0.1, 5.0))

        Tspan = max(self.psr.toas.max(), self.psr2.toas.max()) - min(self.psr.toas.max(), self.psr2.toas.max())

        pl = utils.powerlaw(log10_A=parameter.Uniform(-18, -12), gamma=parameter.Uniform(1, 7))

        rn = gp_signals.FourierBasisCommonGP(spectrum=pl, orf=utils.hd_orf(), components=20, Tspan=Tspan)

        model = ef + rn

        rnc = gp_signals.FourierBasisCommonGP(
            spectrum=pl, orf=utils.hd_orf(), components=20, Tspan=Tspan, coefficients=True
        )

        modelc = ef + rnc

        pta = signal_base.PTA([model(self.psr), model(self.psr2)])
        ptac = signal_base.PTA([modelc(self.psr), modelc(self.psr2)])

        params = {
            "B1855+09_efac": 1.0,
            "B1937+21_efac": 1.0,
            "common_fourier_gamma": 5,
            "common_fourier_log10_A": -15,
        }

        # get GP delays in two different ways

        cf, cf2 = np.random.randn(40), np.random.randn(40)

        bs = pta.get_basis(params)
        da = [np.dot(bs[0], cf), np.dot(bs[1], cf2)]

        params.update({"B1855+09_common_fourier_coefficients": cf, "B1937+21_common_fourier_coefficients": cf2})

        db = ptac.get_delay(params)

        msg = "Implicit and explicit GP delays are different."
        assert np.allclose(da[0], db[0]), msg
        assert np.allclose(da[1], db[1]), msg

        cpar = [p for p in ptac.params if "coefficients" in p.name]

        def shouldfail():
            return cpar[0].get_logpdf(params)

        self.assertRaises(NotImplementedError, shouldfail)

    def test_fourier_red_noise(self):
        """Test that implicit and explicit GP delays are the same."""
        # set up signal parameter
        pl = utils.powerlaw(log10_A=parameter.Uniform(-18, -12), gamma=parameter.Uniform(1, 7))
        rn = gp_signals.FourierBasisGP(spectrum=pl, components=20)
        rnm = rn(self.psr)

        rnc = gp_signals.FourierBasisGP(spectrum=pl, components=20, coefficients=True)
        rnmc = rnc(self.psr)

        # parameters
        log10_A, gamma = -14.5, 4.33
        params = {"B1855+09_red_noise_log10_A": log10_A, "B1855+09_red_noise_gamma": gamma}

        # get the GP delays in two different ways
        cf = np.random.randn(40)
        d1 = np.dot(rnm.get_basis(params), cf)

        params.update({"B1855+09_red_noise_coefficients": cf})
        d2 = rnmc.get_delay(params)

        msg = "Implicit and explicit GP delays are different."
        assert np.allclose(d1, d2), msg

        # np.array cast is needed because we get a KernelArray
        phimat = np.array(rnm.get_phi(params))
        pr1 = -0.5 * np.sum(cf * cf / phimat) - 0.5 * np.sum(np.log(phimat)) - 0.5 * len(phimat) * np.log(2 * np.pi)

        cpar = [p for p in rnmc.params if "coefficients" in p.name][0]
        pr2 = cpar.get_logpdf(params=params)

        msg = "Implicit and explicit GP priors are different."
        assert np.allclose(pr1, pr2), msg

    def test_ecorr_backend(self):
        """Test that ecorr-backend signal returns correct values."""
        # set up signal parameter
        ecorr = parameter.Uniform(-10, -5)
        selection = Selection(selections.by_backend)
        ec = gp_signals.EcorrBasisModel(log10_ecorr=ecorr, selection=selection)
        ecm = ec(self.psr)

        ecc = gp_signals.EcorrBasisModel(log10_ecorr=ecorr, selection=selection, coefficients=True)
        eccm = ecc(self.psr)

        # parameters
        ecorrs = [-6.1, -6.2, -6.3, -6.4]
        params = {
            "B1855+09_basis_ecorr_430_ASP_log10_ecorr": ecorrs[0],
            "B1855+09_basis_ecorr_430_PUPPI_log10_ecorr": ecorrs[1],
            "B1855+09_basis_ecorr_L-wide_ASP_log10_ecorr": ecorrs[2],
            "B1855+09_basis_ecorr_L-wide_PUPPI_log10_ecorr": ecorrs[3],
        }

        fmat = ecm.get_basis(params)
        cf = 1e-6 * np.random.randn(fmat.shape[1])
        d1 = np.dot(fmat, cf)

        for key in ecm._keys:
            parname = "B1855+09_basis_ecorr_" + key + "_coefficients"
            params[parname] = cf[ecm._slices[key]]
        d2 = eccm.get_delay(params)

        msg = "Implicit and explicit ecorr-basis delays are different."
        assert np.allclose(d1, d2), msg

    def test_formalism(self):
        # create marginalized model
        ef = white_signals.MeasurementNoise(efac=parameter.Uniform(0.1, 5.0))
        tm = gp_signals.TimingModel()
        ec = gp_signals.EcorrBasisModel(log10_ecorr=parameter.Uniform(-10, -5))
        pl = utils.powerlaw(log10_A=parameter.Uniform(-18, -12), gamma=parameter.Uniform(1, 7))
        rn = gp_signals.FourierBasisGP(spectrum=pl, components=10)
        model = ef + tm + ec + rn
        pta = signal_base.PTA([model(self.psr)])

        # create hierarchical model
        tmc = gp_signals.TimingModel(coefficients=True)
        ecc = gp_signals.EcorrBasisModel(log10_ecorr=parameter.Uniform(-10, -5), coefficients=True)
        rnc = gp_signals.FourierBasisGP(spectrum=pl, components=10, coefficients=True)
        modelc = ef + tmc + ecc + rnc
        ptac = signal_base.PTA([modelc(self.psr)])

        ps = {
            "B1855+09_efac": 1,
            "B1855+09_basis_ecorr_log10_ecorr": -6,
            "B1855+09_red_noise_log10_A": -14,
            "B1855+09_red_noise_gamma": 3,
        }
        psc = utils.get_coefficients(pta, ps)

        d1 = ptac.get_delay(psc)[0]
        d2 = (
            np.dot(pta.pulsarmodels[0].signals[1].get_basis(ps), psc["B1855+09_linear_timing_model_coefficients"])
            + np.dot(pta.pulsarmodels[0].signals[2].get_basis(ps), psc["B1855+09_basis_ecorr_coefficients"])
            + np.dot(pta.pulsarmodels[0].signals[3].get_basis(ps), psc["B1855+09_red_noise_coefficients"])
        )

        msg = "Implicit and explicit PTA delays are different."
        assert np.allclose(d1, d2), msg

        l1 = pta.get_lnlikelihood(ps)
        l2 = ptac.get_lnlikelihood(psc)

        # I don't know how to integrate l2 to match l1...
        msg = "Marginal and hierarchical likelihoods should be different."
        assert l1 != l2, msg

    def test_conditional_gp(self):
        ef = white_signals.MeasurementNoise(efac=parameter.Uniform(0.1, 5.0))
        tm = gp_signals.TimingModel()
        ec = gp_signals.EcorrBasisModel(log10_ecorr=parameter.Uniform(-10, -5))
        pl = utils.powerlaw(log10_A=parameter.Uniform(-18, -12), gamma=parameter.Uniform(1, 7))
        rn = gp_signals.FourierBasisGP(spectrum=pl, components=10, combine=False)

        model = ef + tm + ec + rn
        pta = signal_base.PTA([model(self.psr), model(self.psr2)])

        p0 = {
            "B1855+09_basis_ecorr_log10_ecorr": -6.051740765663904,
            "B1855+09_efac": 2.9027266737466095,
            "B1855+09_red_noise_gamma": 6.9720332277819725,
            "B1855+09_red_noise_log10_A": -16.749192700991543,
            "B1937+21_basis_ecorr_log10_ecorr": -9.726747733721872,
            "B1937+21_efac": 3.959178240268702,
            "B1937+21_red_noise_gamma": 2.9030772884814797,
            "B1937+21_red_noise_log10_A": -17.978562921948992,
        }

        c = utils.ConditionalGP(pta)
        cmean = c.get_mean_coefficients(p0)

        # build index for the global coefficient vector
        idx, ntot = {}, 0
        for c_name, v in cmean.items():
            idx[c_name] = slice(ntot, ntot + len(v))
            ntot = ntot + len(v)

        # repeat the computation using the common-signal formalism
        TNrs = pta.get_TNr(p0)
        TNTs = pta.get_TNT(p0)
        phiinvs = pta.get_phiinv(p0, logdet=False, method="cliques")

        TNr = np.concatenate(TNrs)
        Sigma = sps.block_diag(TNTs, "csc") + sps.block_diag([np.diag(phiinvs[0]), np.diag(phiinvs[1])])

        ch = cholesky(Sigma)
        mn = ch(TNr)
        iSigma = sps.linalg.inv(Sigma)

        # check mean values
        msg = "Conditional GP coefficient value does not match"
        for c_name, v in cmean.items():
            assert np.allclose(mn[idx[c_name]], v, atol=1e-4, rtol=1e-4), msg

        # check variances
        par = "B1937+21_linear_timing_model_coefficients"
        c1 = np.cov(np.array([cs[par] for cs in c.sample_coefficients(p0, n=10000)]).T)
        c2 = iSigma[idx[par], idx[par]].toarray().T
        msg = "Conditional GP coefficient variance does not match"
        assert np.allclose(c1, c2, atol=1e-4, rtol=1e-4), msg

        # check mean processes
        proc = "B1937+21_linear_timing_model"
        p1 = c.get_mean_processes(p0)[proc]
        p2 = np.dot(pta["B1937+21"]["linear_timing_model"].get_basis(), mn[idx[par]])
        msg = "Conditional GP time series does not match"
        assert np.allclose(p1, p2, atol=1e-4, rtol=1e-4), msg

        # check mean of sampled processes
        p2 = np.mean(np.array([pc[proc] for pc in c.sample_processes(p0, n=1000)]), axis=0)
        msg = "Mean of sampled conditional GP processes does not match"
        assert np.allclose(p1, p2, atol=1e-4, rtol=1e-4)

        # now try with a common process

        crn = gp_signals.FourierBasisCommonGP(spectrum=pl, orf=utils.hd_orf(), components=10, combine=False)

        model = ef + tm + ec + crn
        pta = signal_base.PTA([model(self.psr), model(self.psr2)])

        p0 = {
            "B1855+09_basis_ecorr_log10_ecorr": -5.861847220080768,
            "B1855+09_efac": 4.588342210948306,
            "B1937+21_basis_ecorr_log10_ecorr": -9.151872649912377,
            "B1937+21_efac": 0.8947815819783302,
            "common_fourier_gamma": 6.638289750637263,
            "common_fourier_log10_A": -15.68180643904114,
        }

        c = utils.ConditionalGP(pta)
        cmean = c.get_mean_coefficients(p0)

        idx, ntot = {}, 0
        for c_name, v in cmean.items():
            idx[c_name] = slice(ntot, ntot + len(v))
            ntot = ntot + len(v)

        TNrs = pta.get_TNr(p0)
        TNTs = pta.get_TNT(p0)
        phiinvs = pta.get_phiinv(p0, logdet=False, method="cliques")

        TNr = np.concatenate(TNrs)
        Sigma = sps.block_diag(TNTs, "csc") + sps.csc_matrix(phiinvs)

        ch = cholesky(Sigma)
        mn = ch(TNr)

        msg = "Conditional GP coefficient value does not match for common GP"
        for c_name, v in cmean.items():
            assert np.allclose(mn[idx[c_name]], v)


@pytest.mark.skipif(not PINT_INSTALLED, reason="Skipping tests that require PINT because it isn't installed")
class TestGPCoefficientsPint(TestGPCoefficients):
    @classmethod
    def setUpClass(cls):
        """Setup the Pulsar object."""

        # initialize Pulsar class
        cls.psr = Pulsar(
            datadir + "/B1855+09_NANOGrav_9yv1.gls.par",
            datadir + "/B1855+09_NANOGrav_9yv1.tim",
            ephem="DE430",
            timing_package="pint",
        )

    def test_ephemeris(self):
        # skipping ephemeris with PINT
        pass


@pytest.mark.skipif(not LIBSTEMPO_INSTALLED, reason="Skipping tests that require libstempo because it isn't installed")
class TestGPCoefficientsTempo2(TestGPCoefficients):
    @classmethod
    def setUpClass(cls):
        """Setup the Pulsar object."""

        # initialize Pulsar class
        cls.psr = Pulsar(datadir + "/B1855+09_NANOGrav_9yv1.gls.par", datadir + "/B1855+09_NANOGrav_9yv1.tim")

        cls.psr2 = Pulsar(datadir + "/B1937+21_NANOGrav_9yv1.gls.par", datadir + "/B1937+21_NANOGrav_9yv1.tim")
