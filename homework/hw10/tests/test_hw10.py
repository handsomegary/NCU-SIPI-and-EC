"""Behavioral checks for mixed modes, causal pulses, convolution and DFE."""
import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

spec = importlib.util.spec_from_file_location('hw10', Path(__file__).resolve().parents[1] / 'src' / 'hw10.py')
hw10 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hw10)


class DifferentialChannelTests(unittest.TestCase):
    def test_mode_conversion_matches_excited_waves(self):
        rng = np.random.default_rng(12)
        s = rng.normal(size=(7,4,4)) + 1j*rng.normal(size=(7,4,4))
        mm = hw10.mixed_mode(s)
        ad = np.array([1,0,-1,0])/np.sqrt(2)
        output = s @ ad
        assert_allclose(mm[:,3,1],(output[:,1]-output[:,3])/np.sqrt(2),atol=1e-15)

    def test_ideal_matched_channel_has_unit_received_differential_level(self):
        s = np.zeros((2,4,4),complex)
        s[:,1,0]=s[:,0,1]=s[:,3,2]=s[:,2,3]=1
        assert_allclose(hw10.mixed_mode(s)[:,3,1],1)
        # +/-2 V source differential signal, matched source divisor 2.
        assert_allclose(hw10.mixed_mode(s)[:,3,1]/2*hw10.VPP,1)

    def test_trapezoid_is_causal_and_has_correct_area_and_endpoints(self):
        p = hw10.trapezoid(50)
        self.assertEqual(p[0],0)
        self.assertEqual(p[-1],0)
        assert_allclose(p[10:51],2)
        self.assertAlmostEqual(float(p.sum()*.5e-12),2*25e-12,places=20)
        assert_allclose(p[1:11],np.arange(1,11)*.2,atol=1e-15)

    def test_polyphase_matches_direct_linear_convolution(self):
        rng = np.random.default_rng(4)
        symbols=rng.choice([-1,1],43)
        h=rng.normal(size=37)
        p=8
        impulses=np.zeros((len(symbols)-1)*p+1)
        impulses[::p]=symbols
        expected=np.convolve(impulses,h)
        assert_allclose(hw10.symbol_waveform(symbols,h,p),expected,atol=4e-15)

    def test_ffe_is_applied_once_as_symbol_delays(self):
        symbols=np.array([1,-1,1,1,-1])
        sbr=np.array([0,.2,.5,.1,0])
        taps=np.array([-.2,.8]); p=4
        left=hw10.symbol_waveform(symbols,hw10.apply_ffe(sbr,taps,p),p)
        right=hw10.symbol_waveform(np.convolve(symbols,taps),sbr,p)
        assert_allclose(left,right,atol=1e-15)

    def test_dfe_cancels_postcursor_isi_recursively(self):
        truth=np.array([1,-1,-1,1,1,-1,1,-1],np.int8)
        taps=np.array([.30,-.12,.05])
        received=np.convolve(truth,np.r_[1,taps])[:len(truth)]
        decisions,feedback=hw10.dfe_decisions(received,taps)
        assert_array_equal(decisions,truth)
        assert_allclose(received-feedback,truth,atol=1e-15)
        changed=received.copy();changed[0]=-1
        _,changed_feedback=hw10.dfe_decisions(changed,taps)
        self.assertAlmostEqual(changed_feedback[1],-.30)

    def test_reader_preserves_four_port_row_order(self):
        s=np.arange(1,17,dtype=float).reshape(4,4)
        row=' '.join(f'{v} 0' for v in s.ravel())
        fixture='# Hz S RI R 50\n0 '+row+'\n10000000 '+row+'\n'
        with patch.object(Path,'read_text',return_value=fixture):
            f,loaded,z=hw10.read_touchstone(Path('fixture.s4p'))
        self.assertEqual(z,50)
        assert_array_equal(f,[0,1e7]);assert_array_equal(loaded[0],s)
        self.assertEqual(loaded[0,1,0],5)

    def test_width_requires_zero_threshold_opening(self):
        p=np.array([-2.,-1.,0.,1.,2.])
        one=np.array([-.2,.2,.5,.2,-.2]);zero=-one
        width,left,right=hw10._opening_width(p,one,zero,2)
        assert_allclose([width,left,right],[3,-1.5,1.5])
        # A vertical gap entirely above zero is not a valid fixed-zero eye.
        self.assertEqual(hw10._opening_width(p,one+2,zero+2,2)[0],0)

    def test_equalizer_constraints(self):
        peak,_=hw10.ctle_peak()
        self.assertLess(peak,6)
        self.assertGreater(peak,5.98)
        self.assertAlmostEqual(np.abs(hw10.DEFAULT_FFE).sum(),1)


if __name__=='__main__':
    unittest.main()
