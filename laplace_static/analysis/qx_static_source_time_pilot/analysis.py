"""Position-resolved descriptive diagnostics for the three technical cases."""
import math
from common import (CONFIGS, RADII, CUTOFFS, POSITIONS, BASES, ANCHORS, LABELS,
                    frobenius, normalized_matrix, delta_matrix, phase, ratio)


def summarize_projection(z, baseline, equal_time, baseline_equal_time):
    if equal_time.real <= 0 or baseline_equal_time.real <= 0:
        raise ValueError('Nonpositive equal-time fixed-source reference')
    return dict(Q_RE=z.real, Q_IM=z.imag, phase=phase(z),
                real_sign=1 if z.real > 0 else (-1 if z.real < 0 else 0),
                common_t90_normalized_re=z.real / baseline_equal_time.real,
                common_t90_normalized_im=z.imag / baseline_equal_time.real,
                own_equal_time_normalized_re=z.real / equal_time.real,
                own_equal_time_normalized_im=z.imag / equal_time.real,
                delta_Q_RE=z.real - baseline.real, delta_Q_IM=z.imag - baseline.imag,
                common_t90_normalized_delta_re=(z.real - baseline.real) / baseline_equal_time.real,
                common_t90_normalized_delta_im=(z.imag - baseline.imag) / baseline_equal_time.real,
                complex_relative_change=ratio(abs(z - baseline), abs(baseline)),
                relative_change_defined=baseline != 0)


def compute(data):
    lookup = {(r['configuration'], r['label'], r['R'], r['basis'], r['anchor'], r['tau']): complex(r['Q_RE'], r['Q_IM']) for r in data['projected']}
    projections, checkpoints = [], []
    for cfg in CONFIGS:
        for radius in RADII:
            for basis in BASES:
                for anchor in ANCHORS:
                    for t0 in POSITIONS:
                        zero = lookup[cfg, LABELS[t0], radius, basis, anchor, 0]
                        refzero = lookup[cfg, LABELS[90], radius, basis, anchor, 0]
                        for tau in range(33):
                            z = lookup[cfg, LABELS[t0], radius, basis, anchor, tau]
                            baseline = lookup[cfg, LABELS[90], radius, basis, anchor, tau]
                            row = dict(configuration=cfg, R=radius, basis=basis, anchor=anchor,
                                       t0=t0, t1=t0 + tau, tau=tau, piter=10)
                            row.update(summarize_projection(z, baseline, zero, refzero))
                            projections.append(row)
                            if basis == 'full' and anchor == 4 and tau in (0, 4, 8, 16, 24, 32):
                                checkpoints.append(dict(row))
    matrices = {}
    control = {}
    for row in data['values']:
        if row['axis'] != '-1':
            continue
        cfg, t0, piter, radius, tau = (int(row[k]) for k in ('cfg', 't0', 'piter', 'R', 'tau'))
        key = cfg, t0, radius, tau
        target = control if piter == 40 else matrices
        target.setdefault(key, [[0j] * 4 for _ in range(4)])[CUTOFFS.index(int(row['Nv_sink']))][CUTOFFS.index(int(row['Nv_src']))] = complex(float(row['L_RE']), float(row['L_IM']))
    for (cfg, src, sink, tau, radius), row in data['reference'].items():
        matrices.setdefault((cfg, 90, radius, tau), [[0j] * 4 for _ in range(4)])[CUTOFFS.index(sink)][CUTOFFS.index(src)] = complex(float(row['L_RE']), float(row['L_IM']))
    comparisons, controls = [], []
    for cfg in CONFIGS:
        for radius in RADII:
            diagonal = [matrices[cfg, 90, radius, 0][i][i].real for i in range(4)]
            for tau in range(33):
                baseline = normalized_matrix(matrices[cfg, 90, radius, tau], diagonal)
                for t0 in POSITIONS:
                    original = matrices[cfg, t0, radius, tau]
                    matrix = normalized_matrix(original, diagonal)
                    anti = [[(matrix[i][j] - matrix[j][i].conjugate()) / 2 for j in range(4)] for i in range(4)]
                    delta = frobenius(delta_matrix(matrix, baseline))
                    comparisons.append(dict(configuration=cfg, R=radius, t0=t0, t1=t0 + tau, tau=tau,
                        fixed_metric_matrix_norm=frobenius(matrix),
                        fixed_metric_delta_norm=delta,
                        relative_delta_to_baseline=ratio(delta, frobenius(baseline)),
                        relative_delta_defined=frobenius(baseline) != 0,
                        relative_antihermitian_norm=ratio(frobenius(anti), frobenius(matrix)),
                        relative_imaginary_norm=ratio(frobenius([[z.imag for z in row] for row in matrix]), frobenius(matrix))))
                    if cfg == 650 and t0 == 122:
                        high = control[cfg, 122, radius, tau]
                        difference = delta_matrix(original, high)
                        controls.append(dict(configuration=650, R=radius, t0=122, tau=tau,
                            maximum_absolute_component_difference=max(max(abs(z.real), abs(z.imag)) for row in difference for z in row),
                            relative_frobenius_difference=ratio(frobenius(difference), frobenius(original))))
    timings = []
    for launch in data['launches']:
        timings.append({k: launch[k] for k in ('configuration', 'kind', 'piter', 't0', 'stem', 'seconds', 'exit_code', 'job_id', 'array_task_id')})
    totals = [dict(configuration=cfg, mpi_launches=sum(r['configuration'] == cfg for r in timings),
                   mpi_launch_seconds=math.fsum(r['seconds'] for r in timings if r['configuration'] == cfg),
                   mpi_launch_minutes=math.fsum(r['seconds'] for r in timings if r['configuration'] == cfg) / 60)
              for cfg in CONFIGS]
    matrix_checkpoints = []
    for radius in RADII:
        for tau in (0, 4, 8, 16, 24, 32):
            selected = [r['relative_delta_to_baseline'] for r in comparisons
                        if r['R'] == radius and r['tau'] == tau and r['t0'] != 90]
            matrix_checkpoints.append(dict(R=radius, tau=tau, paired_changes=len(selected),
                minimum_relative_change=min(selected), maximum_relative_change=max(selected),
                interpretation='Range over three configurations and two shifts; not an uncertainty interval.'))
    counts = []
    for radius in RADII:
        for t0 in POSITIONS:
            for basis in BASES:
                for anchor in ANCHORS:
                    selected = [r for r in projections if (r['R'], r['t0'], r['basis'], r['anchor']) == (radius, t0, basis, anchor)]
                    late = [r for r in selected if r['tau'] >= 24]
                    counts.append(dict(R=radius, t0=t0, basis=basis, anchor=anchor,
                        configuration_time_entries=len(selected), late_entries=len(late),
                        late_negative_real_entries=sum(r['real_sign'] < 0 for r in late),
                        late_absolute_phase_above_005=sum(r['phase'] is not None and abs(r['phase']) > .05 for r in late),
                        interpretation='Descriptive counts of correlated configuration/time entries; not independent tests or probabilities.'))
    tables = dict(position_projections=projections, full_anchor4_checkpoints=checkpoints,
                  matrix_position_diagnostics=comparisons, projection_control=controls,
                  matrix_checkpoint_ranges=matrix_checkpoints,
                  launch_timings=timings, launch_totals=totals, descriptive_flags=counts)
    summary = dict(configurations=list(CONFIGS), source_positions=list(POSITIONS), radii=list(RADII),
        source_training_configurations=68, frozen_sources=36, statistical_pilot_configurations=3,
        maximum_saved_mode_relative_residual=max(c['audit']['maximum_relative_residual'] for c in data['cases']),
        audited_slices=sum(c['audit']['slices'] for c in data['cases']),
        audited_modes=sum(c['audit']['modes'] for c in data['cases']),
        fresh_reference_matches=sum(s['fresh_reference_matches'] for c in data['cases'] for s in c['scans']),
        maximum_projection_control_component_difference=max(r['maximum_absolute_component_difference'] for r in controls),
        maximum_projection_control_relative_matrix_difference=max(r['relative_frobenius_difference'] for r in controls),
        launch_totals=totals, table_rows={k: len(v) for k, v in tables.items()},
        inference_scope='Position-resolved technical/descriptive pilot, conditional on frozen training sources; no standard errors, significance tests, energy fits, GEVP or final potential.',
        next_measurement_campaign_submitted=False)
    return tables, summary
