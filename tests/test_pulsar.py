#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_pulsar
----------------------------------

Tests for `pulsar` module. Will eventually want to add tests
for time slicing, PINT integration and pickling.
"""

import sys
import os
import shutil
import unittest
import pickle
import pytest

import numpy as np

from enterprise.pulsar import Pulsar
from tests.enterprise_test_data import datadir
from tests.enterprise_test_data import LIBSTEMPO_INSTALLED, PINT_INSTALLED

if PINT_INSTALLED:
    import pint.models.timing_model
    from pint.models import get_model_and_toas


@pytest.mark.skipif(not LIBSTEMPO_INSTALLED, reason="Skipping tests that require libstempo because it isn't installed")
class TestTimingPackageExceptions(unittest.TestCase):
    def test_unkown_timing_package(self):
        # initialize Pulsar class
        with self.assertRaises(ValueError):
            self.psr = Pulsar(
                datadir + "/B1855+09_NANOGrav_9yv1.gls.par",
                datadir + "/B1855+09_NANOGrav_9yv1.tim",
                timing_package="foobar",
            )

    def test_clk_but_no_bipm(self):
        self.psr = Pulsar(
            datadir + "/B1855+09_NANOGrav_9yv1.gls.par",
            datadir + "/B1855+09_NANOGrav_9yv1.tim",
            clk="TT(BIPM2020)",
            timing_package="pint",
        )


@pytest.mark.skipif(not LIBSTEMPO_INSTALLED, reason="Skipping tests that require libstempo because it isn't installed")
class TestPulsar(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Setup the Pulsar object."""

        # initialize Pulsar class
        cls.psr = Pulsar(
            datadir + "/B1855+09_NANOGrav_9yv1.gls.par", datadir + "/B1855+09_NANOGrav_9yv1.tim", drop_t2pulsar=True
        )

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree("pickle_dir", ignore_errors=True)

    def test_droppsr(self):
        self.psr_nodrop = Pulsar(
            datadir + "/B1855+09_NANOGrav_9yv1.gls.par", datadir + "/B1855+09_NANOGrav_9yv1.tim", drop_t2pulsar=False
        )

        self.psr_nodrop.drop_tempopsr()

        with self.assertRaises(AttributeError):
            _ = self.psr.t2pulsar

    def test_residuals(self):
        """Check Residual shape."""

        msg = "Residuals shape incorrect"
        assert self.psr.residuals.shape == (4005,), msg

    def test_toaerrs(self):
        """Check TOA errors shape."""

        msg = "TOA errors shape incorrect"
        assert self.psr.toaerrs.shape == (4005,), msg

    def test_toas(self):
        """Check TOA shape."""

        msg = "TOA shape incorrect"
        assert self.psr.toas.shape == (4005,), msg

    def test_stoas(self):
        """Check STOA shape."""

        msg = "stoa shape incorrect"
        assert self.psr.stoas.shape == (4005,), msg

    def test_dm(self):
        """Check DM/DMX access."""

        msg = "dm value incorrect"
        assert self.psr.dm == 13.299393, msg

        msg = "dmx struct incorrect (spotcheck)"
        assert len(self.psr.dmx) == 72, msg
        assert self.psr.dmx["DMX_0001"]["DMX"] == 0.015161863, msg
        assert self.psr.dmx["DMX_0001"]["fit"], msg

    def test_freqs(self):
        """Check frequencies shape."""

        msg = "Frequencies shape incorrect"
        assert self.psr.freqs.shape == (4005,), msg

    def test_flags(self):
        """Check flags shape, content, and setting"""

        msg = "Flags shape incorrect"
        assert self.psr.flags["f"].shape == (4005,), msg

        msg = "Flag content or sorting incorrect"
        assert np.all(self.psr._flags["fe"][self.psr._isort] == self.psr.flags["fe"]), msg

        # only possible if flags are stored as dict
        if isinstance(self.psr._flags, dict):
            self.psr.set_flags("name2", self.psr.flags["name"])

            msg = "Setting flags returns incorrect match"
            assert np.all(self.psr.flags["name"] == self.psr.flags["name2"])

    def test_backend_flags(self):
        """Check backend_flags shape and content"""

        msg = "Backend flags shape incorrect"
        assert self.psr.backend_flags.shape == (4005,), msg

        # for the test pulsar, backend should be the same as 'f'
        msg = "Flag content or sorting incorrect"
        assert np.all(self.psr._flags["f"][self.psr._isort] == self.psr.backend_flags), msg

    def test_sky(self):
        """Check Sky location."""

        sky = (1.4023093811712661, 4.9533700839400492)

        msg = "Incorrect sky location"
        assert np.allclose(self.psr.theta, sky[0]), msg
        assert np.allclose(self.psr.phi, sky[1]), msg

    def test_design_matrix(self):
        """Check design matrix shape."""

        msg = "Design matrix shape incorrect."
        assert self.psr.Mmat.shape == (4005, 91), msg

    def test_filter_data(self):
        """Place holder for filter_data tests."""
        assert self.psr.filter_data() is None

    def test_planetssb(self):
        """Place holder for filter_data tests."""
        assert hasattr(self.psr, "planetssb")

    def test_sunssb(self):
        """Place holder for filter_data tests."""
        assert hasattr(self.psr, "sunssb")

    @pytest.mark.skipif(sys.version_info < (3, 8), reason="Requires Python >= 3.8")
    def test_deflate_inflate(self):
        psr = Pulsar(datadir + "/B1855+09_NANOGrav_9yv1.gls.par", datadir + "/B1855+09_NANOGrav_9yv1.tim")

        dm = psr._designmatrix.copy()

        psr.deflate()
        psr.to_pickle()

        with open("B1855+09.pkl", "rb") as f:
            pkl_psr = pickle.load(f)
        pkl_psr.inflate()

        assert np.allclose(dm, pkl_psr._designmatrix)

        del pkl_psr

        psr.destroy()

        with open("B1855+09.pkl", "rb") as f:
            pkl_psr = pickle.load(f)

        with self.assertRaises(FileNotFoundError):
            pkl_psr.inflate()

        os.remove("B1855+09.pkl")

    def test_wrong_input(self):
        """Test exception when incorrect par(tim) file given."""

        with self.assertRaises(IOError) as context:
            Pulsar("wrong.par", "wrong.tim")

            msg = "Cannot find parfile wrong.par or timfile wrong.tim!"
            self.assertTrue(msg in context.exception)

    def test_value_error(self):
        """Test exception when unknown argument is given"""

        with self.assertRaises(ValueError):
            Pulsar(datadir + "/B1855+09_NANOGrav_9yv1.gls.par", datadir + "/B1855+09_NANOGrav_9yv1.time")

    def test_to_feather(self):
        """Test creating feather file from Pulsar method"""

        self.psr.to_feather("test.feather")
        assert os.path.exists("test.feather")

        loaded_psr = Pulsar("test.feather")
        assert np.allclose(self.psr.residuals, loaded_psr.residuals, rtol=1e-10)

        os.remove("test.feather")


@pytest.mark.skipif(not PINT_INSTALLED, reason="Skipping tests that require PINT because it isn't installed")
class TestPulsarPint(TestPulsar):
    @classmethod
    def setUpClass(cls):
        """Setup the Pulsar object."""

        # initialize Pulsar class
        cls.psr = Pulsar(
            datadir + "/B1855+09_NANOGrav_9yv1.gls.par",
            datadir + "/B1855+09_NANOGrav_9yv1.tim",
            ephem="DE430",
            drop_pintpsr=True,
            timing_package="pint",
        )

    def test_droppsr(self):
        self.psr_nodrop = Pulsar(
            datadir + "/B1855+09_NANOGrav_9yv1.gls.par",
            datadir + "/B1855+09_NANOGrav_9yv1.tim",
            ephem="DE430",
            drop_pintpsr=False,
            timing_package="pint",
        )

        self.psr_nodrop.drop_pintpsr()

        with self.assertRaises(AttributeError):
            _ = self.psr_nodrop.model

        with self.assertRaises(AttributeError):
            _ = self.psr_nodrop.parfile

        with self.assertRaises(AttributeError):
            _ = self.psr_nodrop.pint_toas

        with self.assertRaises(AttributeError):
            _ = self.psr_nodrop.timfile

    def test_drop_not_picklable(self):
        self.psr_nodrop = Pulsar(
            datadir + "/B1855+09_NANOGrav_9yv1.gls.par",
            datadir + "/B1855+09_NANOGrav_9yv1.tim",
            ephem="DE430",
            drop_pintpsr=False,
            timing_package="pint",
        )

        self.psr_nodrop.drop_not_picklable()

        with self.assertRaises(AttributeError):
            _ = self.psr_nodrop.model

        with self.assertRaises(AttributeError):
            _ = self.psr_nodrop.pint_toas

    def test_deflate_inflate(self):
        pass

    def test_load_radec_psr(self):
        """Setup the Pulsar object."""

        with self.assertRaises(pint.models.timing_model.TimingModelError):
            # initialize Pulsar class with RAJ DECJ and PMLAMBDA, PMBETA
            Pulsar(
                datadir + "/J0030+0451_RADEC_wrong.par",
                datadir + "/J0030+0451_NANOGrav_9yv1.tim",
                ephem="DE430",
                drop_pintpsr=False,
                timing_package="pint",
            )

    def test_load_radec_psr_mdc(self):
        """Setup the Pulsar object."""

        # initialize Pulsar class with RAJ DECJ so _get_radec can be covered
        psr = Pulsar(
            datadir + "/mdc1/J0030+0451.par",
            datadir + "/mdc1/J0030+0451.tim",
            ephem="DE430",
            drop_pintpsr=False,
            timing_package="pint",
        )

        msg = f"Pulsar not loaded properly {self.psr.Mmat.shape}"
        assert psr.Mmat.shape == (130, 8), msg

    def test_no_planet(self):
        """Test exception when incorrect par(tim) file given."""

        with self.assertRaises(ValueError) as context:
            model, toas = get_model_and_toas(
                datadir + "/J0030+0451_NANOGrav_9yv1.gls.par", datadir + "/J0030+0451_NANOGrav_9yv1.tim", planets=False
            )
            Pulsar(model, toas, planets=True, drop_pintpsr=False)
            msg = "obs_earth_pos is not in toas.table.colnames. Either "
            msg += "`planet` flag is not True in `toas` or further Pint "
            msg += "development to add additional planets is needed."
            self.assertTrue(msg in context.exception)
