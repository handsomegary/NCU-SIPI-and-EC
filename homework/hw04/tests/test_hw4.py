"""Independent correctness checks for the HW4 numerical pipeline.

Run from homework/hw04 with: python -m unittest discover -s tests -v
The SciPy reference uses its separate spectral implementation with detrend disabled;
the waveform oracle is a convolution, not a copy of the ramp construction.
"""
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from scipy.signal import periodogram as scipy_periodogram

HOME = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HOME / "src"))
import hw4


def reference_nrz(bits, previous_bit, samples_per_bit, rise_samples, vpp):
    """Convolve the ideal levels with an exclusive causal moving average."""
    ideal = np.repeat(np.asarray(bits, dtype=float) * vpp - vpp / 2, samples_per_bit)
    if rise_samples == 0:
        return ideal
    prehistory = np.full(rise_samples, previous_bit * vpp - vpp / 2)
    return np.convolve(np.r_[prehistory, ideal], np.ones(rise_samples) / rise_samples,
                       mode="valid")[:-1]


class WaveformTests(unittest.TestCase):
    def test_known_four_transition_cases(self):
        # 0->0, 0->1, 1->1, 1->0 with an explicit expected four-sample edge.
        got = hw4.build_nrz_waveform(np.array([0, 1, 1, 0]), previous_bit=0,
                                    samples_per_bit=6, rise_samples=4, vpp=2)
        expected = np.array([-1]*6 + [-1, -.5, 0, .5, 1, 1] + [1]*6
                            + [1, .5, 0, -.5, -1, -1], dtype=float)
        assert_array_equal(got, expected)

    def test_convolution_oracle_all_edge_lengths(self):
        bits = np.random.default_rng(476).integers(0, 2, 61, dtype=np.uint8)
        for previous in (0, 1):
            for samples in (1, 7, 50):
                for rise in sorted({0, 1, samples // 2, samples}):
                    with self.subTest(previous=previous, samples=samples, rise=rise):
                        got = hw4.build_nrz_waveform(bits, previous_bit=previous,
                                                    samples_per_bit=samples,
                                                    rise_samples=rise, vpp=1.7)
                        expected = reference_nrz(bits, previous, samples, rise, 1.7)
                        assert_allclose(got, expected, rtol=0, atol=8e-16)

    def test_sections_equal_continuous_construction(self):
        bits = np.random.default_rng(477).integers(0, 2, 37, dtype=np.uint8)
        bits[12:14] = (1, 0)  # Force a falling edge on the split.
        bits[22:24] = (0, 1)  # Force a rising edge on the split.
        whole = hw4.build_nrz_waveform(bits, previous_bit=1, samples_per_bit=9,
                                      rise_samples=4, vpp=1.2)
        split = [0, 13, 23, len(bits)]
        parts = [hw4.build_nrz_waveform(bits[a:b], previous_bit=1 if a == 0 else bits[a-1],
                                      samples_per_bit=9, rise_samples=4, vpp=1.2)
                 for a, b in zip(split[:-1], split[1:])]
        assert_array_equal(np.concatenate(parts), whole)

    def test_analytic_edge_energy(self):
        bits = np.random.default_rng(478).integers(0, 2, 211, dtype=np.uint8)
        n, r, vpp, z0 = 50, 10, 1.3, 75
        x = hw4.build_nrz_waveform(bits, previous_bit=0, samples_per_bit=n,
                                  rise_samples=r, vpp=vpp)
        transitions = np.count_nonzero(np.diff(np.r_[0, bits]))
        expected = (vpp/2)**2 / z0 * (1 - transitions * (2/3)*(r-1/r)/(len(bits)*n))
        self.assertAlmostEqual(np.mean(x*x)/z0, expected, delta=2e-17)

    def test_invalid_waveform_parameters(self):
        for kwargs in ({"samples_per_bit":0}, {"samples_per_bit":3.5},
                       {"rise_samples":-1}, {"rise_samples":51},
                       {"rise_samples":2.5}, {"previous_bit":2},
                       {"vpp":0}, {"vpp":-1}, {"vpp":np.nan}):
            with self.subTest(kwargs=kwargs), self.assertRaises((ValueError, TypeError)):
                hw4.build_nrz_waveform(np.array([0,1]), **kwargs)
        for bits in ([0,2], [0,.5], [[0,1],[1,0]], []):
            with self.subTest(bits=bits), self.assertRaises((ValueError, TypeError)):
                hw4.build_nrz_waveform(np.asarray(bits))


class PeriodogramTests(unittest.TestCase):
    def test_independent_scipy_reference_even_and_odd(self):
        rng = np.random.default_rng(479)
        for count, nfft in ((127,127), (127,128), (128,255), (128,256), (257,1024)):
            x = rng.normal(.37, .8, count)
            for win in (None, np.hanning(count), np.linspace(.2,1.2,count)):
                with self.subTest(count=count,nfft=nfft,window="boxcar" if win is None else "array"):
                    w = np.ones(count) if win is None else win
                    f, s = hw4.one_sided_periodogram(x, dt=.002, z0=75, nfft=nfft, window=win)
                    sf, ss = scipy_periodogram(x, fs=500, window=w, nfft=nfft,
                                              detrend=False, return_onesided=True, scaling="density")
                    assert_allclose(f, sf, rtol=2e-15, atol=0)
                    assert_allclose(s, ss/75, rtol=2e-12, atol=2e-20)

    def test_parseval_for_arbitrary_windows_and_padding(self):
        rng = np.random.default_rng(480)
        x = rng.normal(.25,1,193)
        for nfft in (193,194,512,1024):
            for w in (np.ones(len(x)), np.hanning(len(x)), np.linspace(.2,1,len(x))):
                with self.subTest(nfft=nfft):
                    f,s = hw4.one_sided_periodogram(x,dt=2e-12,z0=50,nfft=nfft,window=w)
                    integrated = np.sum(s)*(f[1]-f[0])
                    expected = np.sum((x*w)**2)/np.sum(w*w)/50
                    self.assertAlmostEqual(integrated,expected,delta=2e-16)

    def test_dc_is_not_doubled(self):
        f,s=hw4.one_sided_periodogram(np.full(64,2.),dt=.01,z0=50,nfft=64)
        self.assertAlmostEqual(s[0]*(f[1]-f[0]),4/50,delta=2e-16)
        assert_allclose(s[1:],0,atol=1e-30)

    def test_nyquist_is_not_doubled(self):
        x=3*(-1.)**np.arange(64)
        f,s=hw4.one_sided_periodogram(x,dt=.01,z0=75,nfft=64)
        self.assertAlmostEqual(f[-1],50)
        self.assertAlmostEqual(s[-1]*(f[1]-f[0]),9/75,delta=2e-16)
        assert_allclose(s[:-1],0,atol=1e-30)

    def test_on_bin_sine_has_expected_total_power(self):
        x=1.8*np.sin(2*np.pi*13*np.arange(256)/256)
        f,s=hw4.one_sided_periodogram(x,dt=1e-6,z0=60,nfft=256)
        self.assertEqual(np.argmax(s),13)
        self.assertAlmostEqual(np.sum(s)*(f[1]-f[0]),1.8**2/(2*60),delta=2e-16)

    def test_zero_padding_preserves_sampled_psd_and_integrated_power(self):
        x=np.random.default_rng(481).normal(size=64)
        f1,s1=hw4.one_sided_periodogram(x,dt=.002,z0=50,nfft=64)
        f4,s4=hw4.one_sided_periodogram(x,dt=.002,z0=50,nfft=256)
        assert_allclose(f4[::4],f1,rtol=0,atol=0)
        assert_allclose(s4[::4],s1,rtol=2e-14,atol=1e-20)
        self.assertAlmostEqual(s1.sum()*(f1[1]-f1[0]),s4.sum()*(f4[1]-f4[0]),delta=2e-16)

    def test_window_rescaling_does_not_change_density(self):
        x=np.random.default_rng(482).normal(size=127);w=np.hanning(len(x))
        f1,s1=hw4.one_sided_periodogram(x,.002,window=w,nfft=512)
        f2,s2=hw4.one_sided_periodogram(x,.002,window=17*w,nfft=512)
        assert_array_equal(f1,f2)
        assert_allclose(s1,s2,rtol=4e-14,atol=1e-20)

    def test_invalid_periodogram_parameters(self):
        x=np.arange(8,dtype=float)
        for args,kwargs in (((x,0),{}),((x,-1),{}),((x,np.nan),{}),
                            ((x,.1),{"z0":0}),((x,.1),{"z0":-1}),
                            ((x,.1),{"nfft":7}),((x,.1),{"nfft":8.5}),
                            ((x,.1),{"window":np.zeros(8)}),
                            ((x,.1),{"window":np.ones(7)}),
                            ((np.array([1.,np.nan]),.1),{}),
                            ((np.array([1.,np.inf]),.1),{}),
                            ((np.array([]),.1),{}),((np.ones((2,4)),.1),{})):
            with self.subTest(args=args,kwargs=kwargs), self.assertRaises((ValueError,TypeError)):
                hw4.one_sided_periodogram(*args,**kwargs)


class TheoreticalPSDTests(unittest.TestCase):
    def test_random_nrz_theory_and_dc_convention(self):
        f=np.array([0.,1e9,2.5e9,10e9,50e9])
        tb,tr,vpp,z0=100e-12,20e-12,1.,50.
        expected=(vpp/2)**2*tb/z0*np.sinc(f*tb)**2*np.sinc(f*tr)**2
        expected[1:]*=2
        got=hw4.theoretical_random_nrz_psd(f,bit_duration=tb,rise_time=tr,vpp=vpp,z0=z0)
        assert_allclose(got,expected,rtol=2e-14,atol=1e-100)

    def test_theoretical_scaling_and_zero_rise(self):
        f=np.array([0.,1e9,3e9]);tb=1e-10
        baseline=hw4.theoretical_random_nrz_psd(f,rise_time=0,vpp=1,z0=50)
        expected=tb/(4*50)*np.sinc(f*tb)**2;expected[1:]*=2
        assert_allclose(baseline,expected,rtol=2e-14)
        scaled=hw4.theoretical_random_nrz_psd(f,rise_time=0,vpp=2,z0=100)
        assert_allclose(scaled,2*baseline,rtol=2e-14)


class InputDataTests(unittest.TestCase):
    def test_load_bits_whitespace_and_invalid_characters(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as d:
            p=Path(d)/"bits.txt";p.write_text(" 01\n10\t0\r\n")
            assert_array_equal(hw4.load_bits(p),[0,1,1,0,0])
            self.assertEqual(hw4.load_bits(p).dtype,np.uint8)
            for text in ("", "\n \t", "0102", "0,1", "0x01"):
                p.write_text(text)
                with self.subTest(text=text),self.assertRaises(ValueError):hw4.load_bits(p)

    def test_supplied_bit_statistics_and_hybrid_structure(self):
        expected={1:(99981,99843,17),2:(100294,121900,6),3:(100007,110130,9)}
        for number,(ones,transitions,longest) in expected.items():
            b=hw4.load_bits(HOME/"instructions"/f"HW4_bit_pattern_{number}.txt")
            self.assertEqual(len(b),200000);self.assertEqual(int(b.sum()),ones)
            cuts=np.flatnonzero(np.diff(b))+1
            self.assertEqual(len(cuts),transitions)
            self.assertEqual(int(np.diff(np.r_[0,cuts,len(b)]).max()),longest)
        b=hw4.load_bits(HOME/"instructions"/"HW4_bit_pattern_3.txt").reshape(-1,10)
        assert_array_equal(b[:,8],0);assert_array_equal(b[:,9],1)
        self.assertEqual(len(np.unique(b[:,:8],axis=0)),256)

    def test_supplied_pattern2_disparity_anomaly(self):
        # Finite source data are preserved; a plot is insufficient to certify a line code.
        b=hw4.load_bits(HOME/"instructions"/"HW4_bit_pattern_2.txt").astype(np.int16)
        d=2*b.reshape(-1,10).sum(axis=1)-10
        self.assertEqual(dict(zip(*np.unique(d,return_counts=True))),{-2:4699,0:10427,2:4755,4:119})
        self.assertEqual(int((2*b-1).sum()),588)
        self.assertEqual(int(np.cumsum(2*b-1).max()),590)
        self.assertEqual("".join(map(str,b[7780:7790])),"1001111011")


class CompleteAnalysisTests(unittest.TestCase):
    def test_checkpoints_against_independent_continuous_waveform_and_scipy(self):
        bits=np.random.default_rng(483).integers(0,2,25,dtype=np.uint8)
        bits[4:6]=(1,0);bits[9:11]=(0,1)
        config=hw4.Config(bit_duration=1.,rise_time=.25,dt=.125,vpp=1.3,z0=75,
                          block_bits=5,nfft=64,kaiser_beta=1.1,averages=(1,3,5))
        # Construct all 25 bits at once by a convolution and then split the samples.
        whole=reference_nrz(bits,0,8,2,1.3).reshape(5,40)
        window=np.kaiser(40,1.1)
        reference=[]
        for row in whole:
            sf,ss=scipy_periodogram(row,fs=8,window=window,nfft=64,
                                   detrend=False,scaling="density")
            reference.append(ss/75)
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as d:
            p=Path(d)/"bits.txt";p.write_text("".join(map(str,bits)))
            result=hw4.analyze_file(p,config)
        assert_allclose(result["frequency_hz"],sf,rtol=0,atol=0)
        self.assertEqual(set(result["estimates"]),{1,3,5})
        for k in config.averages:
            with self.subTest(k=k):
                expected=np.mean(reference[:k],axis=0)
                assert_allclose(result["estimates"][k],expected,rtol=2e-13,atol=1e-20)
                cp=result["checkpoints"][str(k)]
                self.assertEqual(cp["bits_used"],k*5)
                self.assertAlmostEqual(cp["mean_time_domain_power_w"],
                                       np.mean(whole[:k]**2)/75,delta=2e-17)
                self.assertAlmostEqual(cp["integrated_full_band_psd_w"],
                                       expected.sum()*(sf[1]-sf[0]),delta=2e-17)

    def test_config_rejects_invalid_physics_and_discrete_counts(self):
        for kwargs in ({"dt":0},{"dt":3e-12},{"rise_time":101e-12},
                       {"z0":np.inf},{"block_bits":0},{"block_bits":2.5},
                       {"nfft":524288.5},{"nfft":1000},{"kaiser_beta":-1},
                       {"averages":()},{"averages":(0,)},{"averages":(2.5,)},
                       {"averages":(1,1)}):
            with self.subTest(kwargs=kwargs),self.assertRaises((ValueError,TypeError)):
                hw4.Config(**kwargs).validate()


if __name__=="__main__":unittest.main()
