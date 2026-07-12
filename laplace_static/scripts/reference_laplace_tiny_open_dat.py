#!/usr/bin/env python3
"""Independent DAT/openQCD reference for the tiny Laplace correlator."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np


NT = 8
NS = 4
AVAILABLE_NV = 16
CFG_ID = 1


def parse_args() -> argparse.Namespace:
    default_root = Path(
        os.environ.get(
            "TINY_DIST_ROOT",
            "/home/m2130292/Masterarbeit/laplace_static/validation/"
            "tiny_open8x4x4x4/data",
        )
    )

    parser = argparse.ArgumentParser(
        description=(
            "Evaluate the tiny open-boundary Laplace correlator independently "
            "from the raw DAT eigenvectors and openQCD gauge file."
        )
    )
    parser.add_argument("--data-root", type=Path, default=default_root)
    parser.add_argument("--nvecs", type=int, default=10)
    parser.add_argument("--t-src", type=int, default=2)
    parser.add_argument("--t-sink", type=int, default=6)
    parser.add_argument("--r-sep", type=int, default=2)
    parser.add_argument("--print-site-values", action="store_true")
    return parser.parse_args()


def load_openqcd_gauge(path: Path) -> tuple[np.ndarray, float]:
    """Load qcd_GF_OPENQCD into U[t,x,y,z,mu,a,b]."""

    with path.open("rb") as stream:
        dims = np.fromfile(stream, dtype="<i4", count=4)
        plaquette = np.fromfile(stream, dtype="<f8", count=1)
        payload = np.fromfile(stream, dtype="<f8")

    expected_dims = np.array([NT, NS, NS, NS], dtype=np.int32)
    if dims.shape != (4,) or not np.array_equal(dims, expected_dims):
        raise ValueError(f"unexpected gauge dimensions {dims.tolist()}")
    if plaquette.shape != (1,):
        raise ValueError("missing openQCD plaquette header")

    odd_sites = NT * NS**3 // 2
    expected_doubles = odd_sites * 8 * 3 * 3 * 2
    if payload.size != expected_doubles:
        raise ValueError(
            f"gauge payload has {payload.size} doubles, expected {expected_doubles}"
        )

    raw = payload.reshape(odd_sites, 8, 3, 3, 2)
    blocks = raw[..., 0] + 1j * raw[..., 1]

    gauge = np.empty((NT, NS, NS, NS, 4, 3, 3), dtype=np.complex128)
    assigned = np.zeros((NT, NS, NS, NS, 4), dtype=bool)

    block_index = 0
    extents = (NT, NS, NS, NS)
    for t in range(NT):
        for x in range(NS):
            for y in range(NS):
                for z in range(NS):
                    coord = [t, x, y, z]
                    if sum(coord) % 2 != 1:
                        continue

                    block = blocks[block_index]
                    block_index += 1

                    for mu in range(4):
                        here = tuple(coord)
                        gauge[here + (mu,)] = block[2 * mu]
                        assigned[here + (mu,)] = True

                        previous = coord.copy()
                        previous[mu] = (previous[mu] - 1) % extents[mu]
                        previous_tuple = tuple(previous)
                        gauge[previous_tuple + (mu,)] = block[2 * mu + 1]
                        assigned[previous_tuple + (mu,)] = True

    if block_index != odd_sites or not np.all(assigned):
        raise RuntimeError("openQCD gauge reconstruction did not assign every link")

    return gauge, float(plaquette[0])


def load_dat_eigenvectors(path: Path, nvecs: int) -> np.ndarray:
    """Load mode-major complex doubles as v[mode,x,y,z,color]."""

    raw = np.fromfile(path, dtype="<f8")
    expected_doubles = AVAILABLE_NV * NS**3 * 3 * 2
    if raw.size != expected_doubles:
        raise ValueError(
            f"eigenvector file has {raw.size} doubles, expected {expected_doubles}"
        )

    complex_values = raw.reshape(-1, 2)
    eigenvectors = (
        complex_values[:, 0] + 1j * complex_values[:, 1]
    ).reshape(AVAILABLE_NV, NS, NS, NS, 3)

    selected = eigenvectors[:nvecs]
    flat = selected.reshape(nvecs, -1)
    gram = flat.conj() @ flat.T
    residual = float(np.max(np.abs(gram - np.eye(nvecs))))
    if residual > 1.0e-8:
        raise ValueError(
            "DAT eigenvectors fail the orthonormality check; "
            f"max residual={residual:.6e}"
        )

    return selected


def compute_tau(
    gauge: np.ndarray,
    v_src: np.ndarray,
    v_sink: np.ndarray,
    t_src: int,
    t_sink: int,
) -> np.ndarray:
    nvecs = v_src.shape[0]
    tau = np.empty((NS, NS, NS, nvecs, nvecs), dtype=np.complex128)

    for x in range(NS):
        for y in range(NS):
            for z in range(NS):
                transporter = np.eye(3, dtype=np.complex128)
                for time in range(t_src, t_sink):
                    transporter = transporter @ gauge[time, x, y, z, 0]

                source_modes = v_src[:, x, y, z, :].T
                sink_modes = v_sink[:, x, y, z, :].T
                tau[x, y, z] = (
                    source_modes.conj().T @ transporter @ sink_modes
                )

    return tau


def pair_value(
    tau: np.ndarray,
    source: tuple[int, int, int],
    sink: tuple[int, int, int],
) -> complex:
    return complex(np.vdot(tau[source], tau[sink]))


def main() -> None:
    args = parse_args()

    if not 1 <= args.nvecs <= AVAILABLE_NV:
        raise ValueError(f"--nvecs must be in [1,{AVAILABLE_NV}]")
    if not 0 < args.t_src < args.t_sink < NT - 1:
        raise ValueError("require 0 < --t-src < --t-sink < 7")
    if not 0 <= args.r_sep <= NS // 2:
        raise ValueError(f"--r-sep must be in [0,{NS // 2}]")

    gauge_path = args.data_root / "cnfg" / "open8x4x4x4n1"
    src_path = (
        args.data_root
        / "disteigvecs"
        / f"disteigvecs_open8x4x4x4n1_t{args.t_src}.dat"
    )
    sink_path = (
        args.data_root
        / "disteigvecs"
        / f"disteigvecs_open8x4x4x4n1_t{args.t_sink}.dat"
    )

    gauge, plaquette = load_openqcd_gauge(gauge_path)
    v_src = load_dat_eigenvectors(src_path, args.nvecs)
    v_sink = load_dat_eigenvectors(sink_path, args.nvecs)
    tau = compute_tau(gauge, v_src, v_sink, args.t_src, args.t_sink)

    print(
        "REF_META "
        f"cfg={CFG_ID} Nt={NT} Ns={NS} Nv={args.nvecs} r={args.r_sep} "
        f"t_src={args.t_src} t_sink={args.t_sink} "
        f"tau={args.t_sink - args.t_src} "
        "temporal_bc=open spatial_bc=periodic rho=constant "
        f"header_plaquette={plaquette:+.16e}"
    )
    print(f"REF_INPUT gauge={gauge_path}")
    print(f"REF_INPUT evec_src={src_path}")
    print(f"REF_INPUT evec_sink={sink_path}")

    axes = ("x", "y", "z")
    combined: list[complex] = []
    max_pair_conjugacy = 0.0
    max_axis_imaginary = 0.0

    for axis, name in enumerate(axes):
        values: list[complex] = []
        origin_value = 0.0j

        for x in range(NS):
            for y in range(NS):
                for z in range(NS):
                    source = (x, y, z)
                    sink_list = [x, y, z]
                    sink_list[axis] = (sink_list[axis] + args.r_sep) % NS
                    sink = tuple(sink_list)

                    value = pair_value(tau, source, sink)
                    reverse = pair_value(tau, sink, source)
                    max_pair_conjugacy = max(
                        max_pair_conjugacy,
                        abs(value - reverse.conjugate()),
                    )
                    values.append(value)

                    if source == (0, 0, 0):
                        origin_value = value

                    if args.print_site_values:
                        print(
                            f"REF_SITE axis={name} "
                            f"source={x},{y},{z} sink={sink[0]},{sink[1]},{sink[2]} "
                            f"Re={value.real:+.16e} Im={value.imag:+.16e}"
                        )

        axis_average = complex(np.mean(values))
        max_axis_imaginary = max(max_axis_imaginary, abs(axis_average.imag))
        combined.extend(values)

        print(
            f"REF_PAIR axis={name} source=0,0,0 "
            f"Re={origin_value.real:+.16e} Im={origin_value.imag:+.16e}"
        )
        print(
            f"REF_AXIS axis={name} Nsrc={len(values)} "
            f"Re={axis_average.real:+.16e} Im={axis_average.imag:+.16e}"
        )

    combined_average = complex(np.mean(combined))
    relative_imaginary = abs(combined_average.imag) / max(
        1.0, abs(combined_average.real)
    )

    print(
        f"REF_COMBINED Nsrc={len(combined)} "
        f"Re={combined_average.real:+.16e} "
        f"Im={combined_average.imag:+.16e}"
    )
    print(
        "REF_CHECK "
        f"max_pair_conjugacy_abs={max_pair_conjugacy:+.16e} "
        f"max_axis_avg_abs_im={max_axis_imaginary:+.16e} "
        f"combined_rel_im={relative_imaginary:+.16e}"
    )


if __name__ == "__main__":
    main()
