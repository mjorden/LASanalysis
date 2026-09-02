import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from lasanalysis import crossplot_neutron_density, curves, pickett_plot, plot_tracks


@pytest.fixture
def df():
    depth = np.arange(1000.0, 1100.0, 0.5)
    rng = np.random.default_rng(0)
    return pd.DataFrame(
        {
            "GR": rng.uniform(20, 110, depth.size),
            "RT": rng.uniform(1, 500, depth.size),
            "RM": rng.uniform(1, 500, depth.size),
            "RHOB": rng.uniform(2.0, 2.8, depth.size),
            "NPHI": rng.uniform(5, 35, depth.size),
        },
        index=pd.Index(depth, name="DEPT"),
    )


TRACKS = [
    {"curves": ["GR"], "xlim": (0, 150), "fill": "left"},
    {"curves": ["RT", "RM"], "xlim": (1, 1000), "log": True},
    {"curves": ["RHOB", "NPHI"], "xlim": [(1.95, 2.95), (45, -15)], "twin": True},
]


def test_plot_tracks_depth_axis_is_inverted_and_shared(df):
    fig = plot_tracks(df, TRACKS, depth_range=(1020, 1080))
    try:
        main_axes = [ax for ax in fig.axes if ax.get_shared_y_axes().joined(ax, fig.axes[0])]
        # 3 tracks share y; the twin adds a 4th axes that also shares y.
        assert len(fig.axes) == 4
        for ax in fig.axes:
            lo, hi = ax.get_ylim()
            assert (lo, hi) == (1080, 1020), f"depth must increase downward, got {(lo, hi)}"
        assert len(main_axes) >= 3
        # twin track: second curve got its own x axis with its own (reversed) limits
        rhob_ax, nphi_ax = fig.axes[2], fig.axes[3]
        assert rhob_ax.get_xlim() == (1.95, 2.95)
        assert nphi_ax.get_xlim() == (45, -15)
        assert fig.axes[1].get_xscale() == "log"
    finally:
        plt.close(fig)


def test_plot_tracks_defaults_to_full_depth(df):
    fig = plot_tracks(df, TRACKS[:1])
    try:
        assert fig.axes[0].get_ylim() == (df.index.max(), df.index.min())
    finally:
        plt.close(fig)


def test_plot_tracks_rejects_bad_input(df):
    with pytest.raises(KeyError):
        plot_tracks(df, [{"curves": ["NOPE"]}])
    with pytest.raises(ValueError):
        plot_tracks(df, TRACKS, depth_range=(1080, 1020))
    with pytest.raises(ValueError):
        plot_tracks(df, [])
    with pytest.raises(ValueError):
        plot_tracks(df, [{"curves": ["RT", "RM"], "xlim": [(1, 10)]}])


def test_plot_tracks_accepts_lasfile(pearson):
    fig = plot_tracks(pearson, [{"curves": ["GR"]}, {"curves": ["RILD", "RILM"], "log": True}], depth_range=(3400, 4200))
    try:
        assert fig.axes[0].get_ylim() == (4200, 3400)
        assert len(curves(pearson)) == 17444
    finally:
        plt.close(fig)


def test_crossplot_and_pickett_smoke(df):
    ax = crossplot_neutron_density(df["NPHI"], df["RHOB"], color_by=df["GR"])
    try:
        assert ax.get_ylim()[0] > ax.get_ylim()[1]  # density axis inverted
        assert len(ax.lines) == 3  # three matrix lines
    finally:
        plt.close(ax.figure)
    ax = pickett_plot(df["RT"], df["NPHI"] / 100, rw=0.05, sw_lines=(1.0, 0.5))
    try:
        assert ax.get_xscale() == "log" and ax.get_yscale() == "log"
        assert len(ax.lines) == 2
    finally:
        plt.close(ax.figure)
