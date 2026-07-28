#!/usr/bin/env python3
"""Inspect the 10x10 mode-delta Laplace-static correlator kernel."""

from pathlib import Path
import sys

import numpy as np

if len(sys.argv) != 3:
    raise SystemExit(
        "usage: inspect_laplace_delta_kernel.py "
        "DELTA_PILOT_DIR SLURM_ARRAY_JOB_ID"
    )

root, jobid = Path(sys.argv[1]), sys.argv[2]
n_t, n_r, n_v, n_src = 8, 7, 10, 24**3

prefix = f"cfg{jobid}_"
by_cfg = {}
for path in root.glob(f"{prefix}*.out"):
    suffix = path.name[len(prefix) : -len(".out")]
    if not suffix.isdigit():
        continue
    cfg = int(suffix)
    if cfg in by_cfg:
        raise ValueError(f"duplicate output for configuration {cfg}")
    by_cfg[cfg] = path

if not by_cfg:
    raise ValueError(f"{root}: no outputs found for job {jobid}")

expected_cfgs = list(range(1, max(by_cfg) + 1))
if sorted(by_cfg) != expected_cfgs:
    missing = sorted(set(expected_cfgs) - set(by_cfg))
    raise ValueError(
        f"{root}: non-contiguous configuration inventory; missing={missing}"
    )

n_cfg = len(expected_cfgs)
if n_cfg < 2:
    raise ValueError("at least two configurations are required")

data = np.full((n_cfg, n_t, n_r), np.nan + 1j * np.nan)
mat = np.full(
    (n_cfg, n_t, n_r, n_v, n_v),
    np.nan + 1j * np.nan,
)

for cfg in range(1, n_cfg + 1):
    path = by_cfg[cfg]
    n_data = n_weighted = 0
    with path.open() as stream:
        for lineno, line in enumerate(stream, 1):
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "DATA":
                if len(fields) != 7:
                    raise ValueError(f"{path}:{lineno}: malformed DATA")
                _, c, time, radius, sources, real, imag = fields
                c, time, radius, sources = map(
                    int, (c, time, radius, sources)
                )
                if (
                    c != cfg
                    or sources != n_src
                    or not (1 <= time <= n_t)
                    or not (0 <= radius < n_r)
                    or np.isfinite(data[cfg - 1, time - 1, radius])
                ):
                    raise ValueError(f"{path}:{lineno}: invalid DATA")
                data[cfg - 1, time - 1, radius] = complex(
                    float(real), float(imag)
                )
                n_data += 1
            elif fields[0] == "WDATA":
                if len(fields) != 9:
                    raise ValueError(f"{path}:{lineno}: malformed WDATA")
                _, c, time, radius, k, ell, sources, real, imag = fields
                c, time, radius, k, ell, sources = map(
                    int, (c, time, radius, k, ell, sources)
                )
                if (
                    c != cfg
                    or sources != n_src
                    or not (1 <= time <= n_t)
                    or not (0 <= radius < n_r)
                    or not (0 <= k < n_v and 0 <= ell < n_v)
                    or np.isfinite(
                        mat[cfg - 1, time - 1, radius, k, ell]
                    )
                ):
                    raise ValueError(f"{path}:{lineno}: invalid WDATA")
                mat[cfg - 1, time - 1, radius, k, ell] = complex(
                    float(real), float(imag)
                )
                n_weighted += 1
    if n_data != n_t * n_r or n_weighted != n_t * n_r * n_v**2:
        raise ValueError(
            f"{path}: DATA={n_data}, WDATA={n_weighted}"
        )

if not np.all(np.isfinite(data)) or not np.all(np.isfinite(mat)):
    raise ValueError("non-finite or missing parsed value")

reconstructed = mat.sum(axis=(-2, -1))
absolute = np.abs(reconstructed - data)
relative = absolute / np.maximum(np.abs(data), 1.0e-300)
if not np.allclose(
    reconstructed, data, rtol=5.0e-12, atol=5.0e-14
):
    raise ValueError("constant-channel reconstruction failed")

def hermitian(matrix):
    return 0.5 * (matrix + matrix.conj().T)

def eigenpairs(matrix):
    values, vectors = np.linalg.eigh(hermitian(matrix))
    return values[::-1], vectors[:, ::-1]

def angle(left, right):
    cosine = np.linalg.svd(
        left.conj().T @ right, compute_uv=False
    ).min()
    return np.degrees(np.arccos(np.clip(cosine, 0.0, 1.0)))

print(f"Validated configurations: {n_cfg}")
print(f"Validated DATA/WDATA rows: {data.size}/{mat.size}")
if n_cfg <= n_v:
    print(
        "Statistical warning: Ncfg <= Nv; subleading eigendirections "
        "are exploratory."
    )
print(
    "Maximum delta reconstruction residual: "
    f"abs={absolute.max():.3e}, rel={relative.max():.3e}"
)
print(
    "\nT R herm     s2/s1    s3/s1    e2/e1    e3/e1    "
    "c2/c1    c3/c1    jk2       jk3       "
    "+2 +3 A2max A3max CV2/e1 zCV2 CV3/e1 zCV3"
)

mean = mat.mean(axis=0)
for time in range(1, 6):
    for radius in range(1, 7):
        block = mat[:, time - 1, radius]
        average = block.mean(axis=0)
        residual = (
            np.linalg.norm(average - average.conj().T)
            / np.linalg.norm(average)
        )
        singular = np.linalg.svd(average, compute_uv=False)
        singular /= singular[0]
        values, vectors = eigenpairs(average)
        if values[0] <= 0:
            raise ValueError(f"T={time}, R={radius}: e1 <= 0")
        ratios = values / values[0]
        diagonal = np.real(np.diag(hermitian(average)))
        if np.all(diagonal > 0):
            corr = hermitian(average) / np.sqrt(
                np.outer(diagonal, diagonal)
            )
            corr_values = np.linalg.eigvalsh(hermitian(corr))[::-1]
            corr_values /= corr_values[0]
        else:
            corr_values = np.full(n_v, np.nan)

        loo_ratios, cv, angles_2, angles_3 = [], [], [], []
        for omitted in range(n_cfg):
            training = (block.sum(axis=0) - block[omitted]) / (n_cfg - 1)
            loo_values, loo_vectors = eigenpairs(training)
            if loo_values[0] <= 0:
                raise ValueError(
                    f"T={time}, R={radius}, omit={omitted + 1}: e1 <= 0"
                )
            loo_ratios.append(loo_values / loo_values[0])
            heldout = hermitian(block[omitted])
            cv.append(
                np.real(
                    np.diag(
                        loo_vectors.conj().T @ heldout @ loo_vectors
                    )
                )
            )
            angles_2.append(angle(vectors[:, :2], loo_vectors[:, :2]))
            angles_3.append(angle(vectors[:, :3], loo_vectors[:, :3]))

        loo_ratios = np.asarray(loo_ratios)
        jk = np.sqrt(
            (n_cfg - 1)
            / n_cfg
            * np.sum(
                (loo_ratios - loo_ratios.mean(axis=0)) ** 2,
                axis=0,
            )
        )
        cv = np.asarray(cv)
        cv_mean = cv.mean(axis=0)
        cv_sem = cv.std(axis=0, ddof=1) / np.sqrt(n_cfg)
        cv_z = np.divide(
            cv_mean,
            cv_sem,
            out=np.full(n_v, np.inf),
            where=cv_sem > 0,
        )

        print(
            f"{time} {radius} {residual:8.1e} "
            f"{singular[1]:8.1e} {singular[2]:8.1e} "
            f"{ratios[1]:+8.1e} {ratios[2]:+8.1e} "
            f"{corr_values[1]:+8.1e} {corr_values[2]:+8.1e} "
            f"{jk[1]:8.1e} {jk[2]:8.1e} "
            f"{np.count_nonzero(loo_ratios[:, 1] > 0):2d} "
            f"{np.count_nonzero(loo_ratios[:, 2] > 0):2d} "
            f"{max(angles_2):5.1f} {max(angles_3):5.1f} "
            f"{cv_mean[1] / values[0]:+8.1e} {cv_z[1]:+5.1f} "
            f"{cv_mean[2] / values[0]:+8.1e} {cv_z[2]:+5.1f}"
        )

print("\nCross-time principal angles relative to T/a=4 (degrees)")
print("R  top1(T3,T5)  top2(T3,T5)  top3(T3,T5)")
for radius in range(1, 7):
    basis = {
        time: eigenpairs(mean[time - 1, radius])[1]
        for time in (3, 4, 5)
    }
    result = []
    for rank in (1, 2, 3):
        result.append(
            (
                angle(basis[4][:, :rank], basis[3][:, :rank]),
                angle(basis[4][:, :rank], basis[5][:, :rank]),
            )
        )
    print(
        f"{radius}  "
        + "  ".join(f"{left:5.1f},{right:5.1f}" for left, right in result)
    )
