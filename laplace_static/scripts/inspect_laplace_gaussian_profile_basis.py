#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import h5py
import numpy as np


SIGMAS = np.array(
    [0.05, 0.0894, 0.1289, 0.1683, 0.2078, 0.2472, 0.2867],
    dtype=float,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Validate production eigenvalue files and inspect the "
            "seven Gaussian rho_i^2 weight profiles."
        )
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("mental/runs_Em1p4_Nv10_qcdnew_full"),
    )
    parser.add_argument(
        "--num-configurations",
        type=int,
        default=100,
    )
    parser.add_argument(
        "--nv",
        type=int,
        default=10,
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if args.num_configurations < 1:
        raise ValueError("--num-configurations must be positive")
    if args.nv < 1:
        raise ValueError("--nv must be positive")

    by_cfg = {}

    for path in args.root.glob("n*/eigenvalues/*.h5"):
        match = re.fullmatch(r"n(\d+)", path.parents[1].name)
        if match:
            by_cfg.setdefault(int(match.group(1)), []).append(path)

    expected = set(range(1, args.num_configurations + 1))
    found = set(by_cfg)
    missing = sorted(expected - found)
    unexpected = sorted(found - expected)
    non_unique = {
        cfg: paths
        for cfg, paths in by_cfg.items()
        if len(paths) != 1
    }

    if missing or unexpected or non_unique:
        print("missing configurations:", missing)
        print("unexpected configurations:", unexpected)
        print("non-unique files:", non_unique)
        raise SystemExit("Eigenvalue-file inventory failed")

    all_lambda = []

    for cfg in range(1, args.num_configurations + 1):
        path = by_cfg[cfg][0]

        with h5py.File(path, "r") as h5:
            dataset = h5["disteigvals"]

            if dataset.ndim != 2 or dataset.shape[1] < args.nv:
                raise ValueError(
                    f"{path}: unexpected eigenvalue shape "
                    f"{dataset.shape}"
                )

            eigenvalues = np.asarray(
                dataset[:, : args.nv],
                dtype=float,
            )
            times = np.asarray(h5["times"][:], dtype=int)
            cfg_id = int(h5["cnfg_id"][()])
            run_name = h5["run_name"][()]

            if isinstance(run_name, bytes):
                run_name = run_name.decode()

        if not np.array_equal(
            times,
            np.arange(eigenvalues.shape[0]),
        ):
            raise ValueError(f"{path}: unexpected times")
        if cfg_id != cfg:
            raise ValueError(
                f"{path}: cnfg_id={cfg_id}, expected {cfg}"
            )
        if run_name != "Em1p4":
            raise ValueError(
                f"{path}: run_name={run_name!r}"
            )
        if not np.all(np.isfinite(eigenvalues)):
            raise ValueError(f"{path}: non-finite eigenvalue")
        if np.any(np.diff(eigenvalues, axis=1) < -1.0e-13):
            raise ValueError(
                f"{path}: eigenvalues are not ordered"
            )

        all_lambda.append(eigenvalues)

    all_lambda = np.stack(all_lambda)

    # These are the rho_i^2 weights that enter Eq. (8) directly.
    weights = np.exp(
        -(all_lambda[None, :, :, :] ** 2)
        / (2.0 * SIGMAS[:, None, None, None] ** 2)
    )

    # Do not normalize this reshape in place: it is a view of weights.
    raw_profiles = weights.reshape(len(SIGMAS), -1)
    profile_norms = np.linalg.norm(
        raw_profiles,
        axis=1,
        keepdims=True,
    )
    normalized_profiles = raw_profiles / profile_norms

    singular_values = np.linalg.svd(
        normalized_profiles,
        compute_uv=False,
    )
    relative = singular_values / singular_values[0]

    print(
        f"Validated files: {args.num_configurations}"
    )
    print("Combined eigenvalue shape:", all_lambda.shape)
    print("Global lambda minimum:", all_lambda.min())
    print("Global lambda maximum:", all_lambda.max())

    print("\nDirect correlator weights rho_i^2:")
    for sigma, block in zip(SIGMAS, weights):
        print(
            f"sigma={sigma:.4f}: "
            f"min={block.min():.8e} "
            f"mean={block.mean():.8e} "
            f"max={block.max():.8e}"
        )

    print(
        "\nNormalized-profile singular values:",
        singular_values,
    )
    print("Relative singular values:", relative)
    print(
        "Condition estimate:",
        singular_values[0] / singular_values[-1],
    )

    for cutoff in (
        1.0e-2,
        1.0e-3,
        1.0e-4,
        1.0e-5,
        1.0e-6,
        1.0e-8,
    ):
        rank = int(np.count_nonzero(relative >= cutoff))
        print(
            "Effective rank at relative cutoff "
            f"{cutoff:.0e}: {rank}"
        )


if __name__ == "__main__":
    main()
