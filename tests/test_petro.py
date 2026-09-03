import numpy as np
import pytest

from lasanalysis.petro import (
    MATRIX_DENSITY,
    archie_sw,
    archie_sw_lines,
    density_porosity,
    frac_to_pu,
    gamma_ray_index,
    neutron_density_crossover,
    pu_to_frac,
    vshale_larionov,
    vshale_linear,
)


def test_gamma_ray_index_clips_and_propagates_nan():
    gr = np.array([0.0, 20.0, 65.0, 110.0, 200.0, np.nan])
    igr = gamma_ray_index(gr, 20, 110)
    np.testing.assert_allclose(igr[:5], [0.0, 0.0, 0.5, 1.0, 1.0])
    assert np.isnan(igr[5])
    with pytest.raises(ValueError):
        gamma_ray_index(gr, 110, 20)


def test_vshale_linear_equals_igr():
    gr = np.linspace(0, 150, 7)
    np.testing.assert_allclose(vshale_linear(gr, 20, 110), gamma_ray_index(gr, 20, 110))


def test_vshale_larionov_bounds_and_ordering():
    gr = np.linspace(20, 110, 11)
    tert = vshale_larionov(gr, 20, 110)
    old = vshale_larionov(gr, 20, 110, older=True)
    assert tert[0] == 0.0 and old[0] == 0.0
    # The published coefficients do not reach exactly 1 at IGR=1:
    # tertiary 0.083*(2**3.7 - 1) = 0.9957, older 0.33*(2**2 - 1) = 0.99.
    assert tert[-1] == pytest.approx(0.083 * (2**3.7 - 1))
    assert old[-1] == pytest.approx(0.99)
    # Larionov is always <= linear (it suppresses low IGR), and tertiary <= older in the middle.
    assert np.all(tert <= vshale_linear(gr, 20, 110) + 1e-12)
    assert np.all(tert[1:-1] < old[1:-1])
    assert np.all((tert >= 0) & (tert <= 1))


def test_density_porosity():
    assert density_porosity(2.65) == 0.0
    assert density_porosity(1.0) == pytest.approx(1.0)
    assert density_porosity(2.71, "limestone") == 0.0
    assert density_porosity(2.45, rho_matrix=2.65, rho_fluid=1.0) == pytest.approx(0.2 / 1.65)
    assert density_porosity(2.80) < 0  # not clipped, on purpose
    with pytest.raises(ValueError):
        density_porosity(2.0, rho_matrix=1.0, rho_fluid=1.0)
    with pytest.raises(KeyError):
        density_porosity(2.0, rho_matrix="granite")
    assert set(MATRIX_DENSITY) >= {"sandstone", "limestone", "dolomite", "salt"}


def test_porosity_unit_helpers():
    np.testing.assert_allclose(pu_to_frac([25.0, 0.0]), [0.25, 0.0])
    np.testing.assert_allclose(frac_to_pu([0.25]), [25.0])


def test_neutron_density_crossover():
    nphi = np.array([0.20, 0.20, 0.20, np.nan])
    dphi = np.array([0.30, 0.20, 0.10, 0.30])
    sep, flag = neutron_density_crossover(nphi, dphi)
    np.testing.assert_allclose(sep[:3], [0.10, 0.0, -0.10])
    assert flag.tolist() == [True, False, False, False]
    _, flag2 = neutron_density_crossover(nphi, dphi, threshold=0.15)
    assert not flag2.any()


def test_archie_sw():
    # phi=0.2, Rt=20, Rw=0.05, a=1, m=n=2  ->  Sw = sqrt(0.05 / (0.04*20)) = 0.25
    assert archie_sw(20.0, 0.2, rw=0.05) == pytest.approx(0.25)
    # Sw > 1 clipped by default, raw with clip=False
    assert archie_sw(0.5, 0.2, rw=0.05) == 1.0
    assert archie_sw(0.5, 0.2, rw=0.05, clip=False) == pytest.approx(np.sqrt(0.05 / (0.04 * 0.5)))
    # non-physical inputs -> NaN, not inf / complex
    out = archie_sw([20.0, 20.0, -5.0, np.nan], [0.2, 0.0, 0.2, 0.2], rw=0.05)
    assert out[0] == pytest.approx(0.25)
    assert np.isnan(out[1:]).all()
    with pytest.raises(ValueError):
        archie_sw(20.0, 0.2, rw=0.0)


def test_fit_water_line_recovers_m_and_rw():
    from lasanalysis.petro import fit_water_line

    rng = np.random.default_rng(1)
    n = 4000
    phi = 10 ** rng.uniform(np.log10(0.06), np.log10(0.35), n)
    # true water line m=2, Rw=0.03; scatter points *above* it (Sw < 1) with a wet floor
    sw = np.where(rng.uniform(size=n) < 0.3, 1.0, rng.uniform(0.2, 1.0, n))
    rt = 0.03 / (phi**2 * sw**2) * 10 ** rng.normal(0, 0.02, n)
    fit = fit_water_line(rt, phi, q=5)
    assert fit["m"] == pytest.approx(2.0, abs=0.1)
    assert fit["rw"] == pytest.approx(0.03, rel=0.15)
    assert fit["n_points"] == n
    assert fit["envelope"].shape[1] == 3 and len(fit["envelope"]) >= 3
    # NaN / non-positive / out-of-range inputs are ignored, not fatal
    rt2 = np.concatenate([rt, [np.nan, -1, 5, 5]])
    phi2 = np.concatenate([phi, [0.1, 0.1, 0.01, 0.9]])
    assert fit_water_line(rt2, phi2)["n_points"] == n
    with pytest.raises(ValueError, match="bins"):
        fit_water_line(rt[:20], phi[:20])


def test_archie_sw_lines_are_consistent_with_archie_sw():
    phi = np.array([0.05, 0.1, 0.3])
    lines = archie_sw_lines(phi, rw=0.05, sw_values=(1.0, 0.5))
    for sw, rt in lines.items():
        np.testing.assert_allclose(archie_sw(rt, phi, rw=0.05, clip=False), sw)
