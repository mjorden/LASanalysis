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
    with pytest.raises(ValueError, match="granite.*sandstone"):
        density_porosity(2.0, rho_matrix="granite")
    assert density_porosity(2.71, " Limestone ") == 0.0  # case/whitespace tolerant
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


def test_neutron_lithology_correction():
    from lasanalysis.petro import neutron_lithology_correction as nlc

    np.testing.assert_allclose(nlc([0.20, 0.10], "limestone", "sandstone"), [0.24, 0.14])
    np.testing.assert_allclose(nlc([0.20], "limestone", "dolomite"), [0.14])
    np.testing.assert_allclose(nlc([0.20], "sandstone", "dolomite"), [0.10])
    np.testing.assert_allclose(nlc([0.20], "Limestone", " limestone "), [0.20])
    assert np.isnan(nlc([np.nan], "limestone", "sandstone")[0])
    with pytest.raises(ValueError, match="salt"):
        nlc([0.2], "limestone", "salt")


def test_rwa_and_rw_pick():
    from lasanalysis.petro import pick_rw_from_rwa, rwa

    np.testing.assert_allclose(rwa([20.0, 5.0], [0.2, 0.1]), [0.8, 0.05])
    assert np.isnan(rwa([20.0, -1.0], [0.0, 0.1])).all()
    rng = np.random.default_rng(3)
    n = 500
    phi = rng.uniform(0.06, 0.3, n)
    sw = np.where(rng.uniform(size=n) < 0.3, 1.0, rng.uniform(0.3, 1.0, n))
    rt = 0.05 / (phi**2 * sw**2) * 10 ** rng.normal(0, 0.02, n)
    vsh = rng.uniform(0, 0.3, n)
    depth = np.linspace(3000, 3500, n)
    r = pick_rw_from_rwa(rt, phi, vsh, depth=depth, vsh_cut=0.15)
    assert r["rw"] == pytest.approx(0.05, rel=0.1)
    assert r["n_points"] == int((vsh < 0.15).sum())
    assert r["interval"] is not None and 3000 <= r["interval"][0] < r["interval"][1] <= 3500
    assert pick_rw_from_rwa(rt, phi)["interval"] is None
    with pytest.raises(ValueError, match="need 20"):
        pick_rw_from_rwa(rt[:5], phi[:5])


def test_pick_rsh():
    from lasanalysis.petro import pick_rsh

    rt = np.array([2.0, 3.0, 4.0, 50.0, 60.0] * 4)   # 12 shale samples
    vsh = np.array([0.9, 0.95, 1.0, 0.1, 0.2] * 4)
    assert pick_rsh(rt, vsh) == 3.0
    with pytest.raises(ValueError):
        pick_rsh(rt[:4], vsh[:4])


def test_shaly_sand_models_reduce_to_archie_and_lower_sw_in_shale():
    from lasanalysis.petro import sw_indonesia, sw_simandoux, water_saturation

    rt, phi, rw, rsh = np.array([20.0, 5.0, 2.0]), np.array([0.2, 0.15, 0.1]), 0.05, 3.0
    arch = archie_sw(rt, phi, rw, clip=False)
    np.testing.assert_allclose(sw_simandoux(rt, phi, 0.0, rw, rsh, clip=False), arch, rtol=1e-9)
    np.testing.assert_allclose(sw_indonesia(rt, phi, 0.0, rw, rsh, clip=False), arch, rtol=1e-9)
    # with shale present, both give lower Sw than Archie (shale conductivity explained away)
    sim = sw_simandoux(rt, phi, 0.4, rw, rsh)
    ind = sw_indonesia(rt, phi, 0.4, rw, rsh)
    sim_raw = sw_simandoux(rt, phi, 0.4, rw, rsh, clip=False)
    ind_raw = sw_indonesia(rt, phi, 0.4, rw, rsh, clip=False)
    assert np.all(sim_raw < arch) and np.all(ind_raw < arch)
    assert np.all(sim <= 1) and np.all(ind <= 1)
    # n != 2 Newton path agrees with the closed form at n = 2 and is close to it nearby
    np.testing.assert_allclose(sw_simandoux(rt, phi, 0.4, rw, rsh, n=2.0000001), sim, rtol=1e-5)
    # the Simandoux solution satisfies the defining equation
    s = sw_simandoux(rt, phi, 0.4, rw, rsh, n=2.3, clip=False)
    np.testing.assert_allclose(phi**2 * s**2.3 / rw + 0.4 * s / rsh, 1 / rt, rtol=1e-6)
    # invalid inputs -> NaN; dispatcher
    assert np.isnan(sw_simandoux(-1.0, 0.2, 0.1, rw, rsh)) and np.isnan(sw_indonesia(5.0, 0.0, 0.1, rw, rsh))
    np.testing.assert_allclose(water_saturation("archie", rt, phi, 0.4, rw), archie_sw(rt, phi, rw))
    np.testing.assert_allclose(water_saturation("Indonesia", rt, phi, 0.4, rw, rsh=rsh), ind)
    with pytest.raises(ValueError, match="needs rsh"):
        water_saturation("simandoux", rt, phi, 0.4, rw)
    with pytest.raises(ValueError, match="unknown Sw model"):
        water_saturation("waxman", rt, phi, 0.4, rw, rsh=rsh)


def test_archie_sw_lines_are_consistent_with_archie_sw():
    phi = np.array([0.05, 0.1, 0.3])
    lines = archie_sw_lines(phi, rw=0.05, sw_values=(1.0, 0.5))
    for sw, rt in lines.items():
        np.testing.assert_allclose(archie_sw(rt, phi, rw=0.05, clip=False), sw)
