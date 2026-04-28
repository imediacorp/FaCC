"""Shared pytest fixtures for FaCC test suite."""

import numpy as np
import pytest


@pytest.fixture
def k_array():
    return np.logspace(-2, 0, 50)


@pytest.fixture
def pk_array(k_array):
    return 1e4 * k_array ** (-2)


@pytest.fixture
def sigma_pk(pk_array):
    return pk_array * 0.05


@pytest.fixture
def hz_csv(tmp_path):
    """Write a minimal H(z) CSV to a temp file."""
    path = tmp_path / "hz.csv"
    path.write_text(
        "z,H,sigma_H\n"
        "0.1,69.0,2.5\n"
        "0.3,78.0,3.0\n"
        "0.5,90.0,4.0\n"
        "0.8,113.0,5.0\n"
        "1.0,131.0,6.0\n"
        "1.5,170.0,8.0\n"
        "2.0,224.0,10.0\n"
    )
    return str(path)


@pytest.fixture
def cmb_csv(tmp_path):
    """Write a minimal CMB low-ℓ CSV to a temp file."""
    ells = np.arange(2, 30)
    cl = 5000 * np.exp(-0.02 * ells) + 50 * np.sin(ells)
    sigma = np.abs(cl) * 0.1 + 10.0
    lines = ["ell,C_ell,sigma\n"]
    for l, c, s in zip(ells, cl, sigma):
        lines.append(f"{l},{c:.4f},{s:.4f}\n")
    path = tmp_path / "cmb.csv"
    path.write_text("".join(lines))
    return str(path)


@pytest.fixture
def pk_csv(tmp_path):
    """Write a minimal P(k) CSV to a temp file."""
    k = np.logspace(-2, -0.3, 40)
    pk = 3e4 * k ** (-2)
    sigma = pk * 0.05
    lines = ["k,Pk,sigma_Pk\n"]
    for ki, pi, si in zip(k, pk, sigma):
        lines.append(f"{ki:.6f},{pi:.4f},{si:.4f}\n")
    path = tmp_path / "pk.csv"
    path.write_text("".join(lines))
    return str(path)
