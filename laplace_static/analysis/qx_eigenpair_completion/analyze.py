#!/usr/bin/env python3
"""Verify qX generation and saved-mode exports and plot measured coverage."""
import argparse
import csv
import hashlib
import importlib.util
import io
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def csv_bytes(rows, fields):
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode()


def placement(folder):
    groups = json.loads((folder/'groups.json').read_text())['groups']
    all_hosts = []
    for i, hosts in enumerate(groups):
        text = (folder/f'placement_{i}.out').read_text()
        rows = [dict(re.findall(r'(\w+)=([^\s]+)', line))
                for line in text.splitlines() if line.startswith('PLACEMENT ')]
        expected = {h.split('.')[0] for h in hosts}
        if len(expected) != 8 or len(rows) != 384 or {int(r['RANK']) for r in rows} != set(range(384)):
            raise ValueError('Incomplete placement')
        if {r['HOST'].split('.')[0] for r in rows} != expected:
            raise ValueError('Unexpected placement hosts')
        cores = {(r['HOST'].split('.')[0], int(r['SOCKET']), int(r['CORE'])) for r in rows}
        if len(cores) != 384 or any(int(r['BAD']) or int(r['AFFINITY_CPUS']) < 1 for r in rows):
            raise ValueError('Invalid physical-core placement')
        if any(sum(r['HOST'].split('.')[0] == h for r in rows) != 48 for h in expected):
            raise ValueError('Wrong ranks per host')
        markers = [dict(re.findall(r'(\w+)=(\d+)', line)) for line in text.splitlines()
                   if line.startswith('PASS_PLACEMENT ')]
        if len(markers) != 1 or any(int(markers[0].get(k, -1)) != v for k, v in
                [('WORLD_RANKS',384),('HOSTS',8),('PHYSICAL_CORES_PER_HOST',48),('NTHREAD',1)]):
            raise ValueError('Missing placement PASS')
        if int(markers[0]['PROVIDED_MIN']) < 1:
            raise ValueError('Insufficient MPI thread support')
        all_hosts.extend(expected)
    if len(set(all_hosts)) != 8*len(groups):
        raise ValueError('Worker node groups overlap')


def verify_and_reduce():
    inputs = HERE/'inputs'
    sources = json.loads((inputs/'SOURCES.json').read_text())
    for name, expected in sources['preserved_inputs'].items():
        if digest(inputs/name) != expected:
            raise ValueError('Input checksum changed: '+name)
    generated = {}
    for entry in sources['generation_jobs']:
        folder = inputs/entry['directory']
        summary = json.loads((folder/'SUMMARY.json').read_text())
        if summary['success'] is not True or int(summary['job_id']) != entry['job_id']:
            raise ValueError('Failed or mismatched generation job')
        placement(folder/'placement')
        if sorted(w['configuration'] for w in summary['workers']) != sorted(summary['configurations']):
            raise ValueError('Incomplete worker inventory')
        for worker in summary['workers']:
            cfg = worker['configuration']
            if cfg in generated or not worker['success'] or worker['mpi_exit_code'] != 0:
                raise ValueError('Duplicate or failed generation')
            run = folder/f'n{cfg}'
            text = (run/'log/mental.log').read_text()+next(run.glob('mental_*.out')).read_text()
            times = [int(x) for x in re.findall(r'Eigenvalues/vectors on x0=(\d+) computed', text)]
            if len(times) != 192 or set(times) != set(range(192)):
                raise ValueError('Incomplete convergence-return timeslice coverage')
            if 'MPI_INIT requested=FUNNELED' not in text or re.search(r'did not converge|MPI_Abort', text, re.I):
                raise ValueError('Failed MPI generation diagnostic')
            done = json.loads((run/'COMPLETED.json').read_text())
            if done != worker['completion'] or (done['configuration'],done['timeslices'],done['eigenvector_files']) != (cfg,192,192):
                raise ValueError('Completion inventory mismatch')
            generated[cfg] = dict(configuration=cfg, generation_job=entry['job_id'],
                                  seconds=worker['mpi_and_full_write_seconds'],
                                  eigenvector_files=192, eigenvector_bytes=done['eigenvector_bytes'])
    audits = []
    all_modes, all_slices = [], []
    audited = set()
    for entry in sources['audit_jobs']:
        folder = inputs/entry['directory']
        spec = importlib.util.spec_from_file_location('saved_mode_export_'+str(entry['job_id']),folder/'check_export.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        summary = json.loads((folder/'SUMMARY.json').read_text())
        launches = json.loads((folder/'LAUNCHES.json').read_text())
        prepared = json.loads((folder/'PREPARED.json').read_text())
        if summary['success'] is not True or [r['configuration'] for r in launches] != module.CONFIGS:
            raise ValueError('Failed audit or launch coverage')
        modes, slices = [], []
        for cfg in module.CONFIGS:
            if cfg in audited:
                raise ValueError('Repeated audited configuration')
            audited.add(cfg)
            launch = next(r for r in launches if r['configuration'] == cfg)
            if launch['mpi_exit_code'] != 0:
                raise ValueError('Saved-mode MPI failure')
            passed, current_modes, current_slices = module.parse(folder/f'n{cfg}.out', cfg)
            if not passed:
                raise ValueError('Saved-mode numerical audit failed')
            run = next(r for r in prepared['runs'] if r['configuration'] == cfg)
            if cfg in generated:
                genjob = generated[cfg]['generation_job']
                marker = inputs/f'job{genjob}'/f'n{cfg}/COMPLETED.json'
                if digest(marker) != run['completion_marker_sha256']:
                    raise ValueError('Audited run differs from generation output')
            residual = max(r['residual_relative'] for r in current_modes)
            claimed = next(r for r in summary['runs'] if r['configuration'] == cfg)
            if not claimed['success'] or claimed['mpi_exit_code'] != 0 or residual != claimed['maximum_relative_residual']:
                raise ValueError('Audit summary reduction mismatch')
            above = sum(r['residual_relative'] > 1e-9 for r in current_modes)
            if above != claimed['modes_above_requested_1e_minus9']:
                raise ValueError('Requested-residual counter mismatch')
            audits.append(dict(configuration=cfg, audit_job=entry['job_id'], selected_times=len(module.TIMES),
                modes=len(current_modes), max_relative_residual=residual,
                max_central_relative_residual=max(r['residual_relative'] for r in current_modes if 95<=r['t']<=101),
                max_norm_defect=max(r['max_norm_defect'] for r in current_slices),
                max_gram_diagonal_defect=max(r['gram_diag_defect'] for r in current_slices),
                max_gram_offdiagonal=max(r['gram_offdiag'] for r in current_slices),
                max_rayleigh_im=max(r['max_rayleigh_im'] for r in current_slices),
                modes_above_requested_1e_minus9=above))
            modes.extend(current_modes)
            slices.extend(current_slices)
        if csv_bytes(modes,['configuration','t','mode']+module.MODE_FIELDS) != (folder/'modes.csv').read_bytes():
            raise ValueError('Mode CSV differs from parsed MPI output')
        if csv_bytes(slices,['configuration','t','status']+module.SLICE_FIELDS) != (folder/'slices.csv').read_bytes():
            raise ValueError('Slice CSV differs from parsed MPI output')
        if (len(modes),len(slices)) != (summary['mode_records'],summary['slice_records']):
            raise ValueError('CSV inventory differs from summary')
        all_modes.extend(modes)
        all_slices.extend(slices)
    prior = sources['prior_generation']['configurations']
    completed = sorted(set(prior)|audited|set(generated))
    expected = list(range(100,801,10))
    if sorted(set(expected)-set(completed)) != [640] or set(generated) != set(range(680,801,10)):
        raise ValueError('Unexpected generation coverage')
    profile_cases = set()
    for source in sources['profile_sources']:
        path = REPO/source['path']
        if digest(path) != source['sha256']:
            raise ValueError('Archived profile source changed')
        with path.open(newline='') as stream:
            profile_cases.update((int(r['configuration']),int(r['tau']),int(r['piter'])) for r in csv.DictReader(stream))
    profiles = sorted({c for c,_,_ in profile_cases})
    if profile_cases != {(c,t,p) for c in profiles for t in (2,4,6) for p in (4,40)}:
        raise ValueError('Incomplete archived paired-profile cohort')
    summary = dict(date=sources['date'], expected_configurations=expected,
        generated_configurations=completed, generated_count=len(completed), unresolved_configuration=640,
        newly_generated_configurations=sorted(generated), newly_generated_count=len(generated),
        freshly_audited_configurations=sorted(audited), freshly_audited_count=len(audited),
        audited_mode_records=len(all_modes), audited_slice_records=len(all_slices),
        max_selected_relative_residual=max(r['max_relative_residual'] for r in audits),
        max_central_relative_residual=max(r['max_central_relative_residual'] for r in audits),
        max_gram_offdiagonal=max(r['max_gram_offdiagonal'] for r in audits),
        max_norm_defect=max(r['max_norm_defect'] for r in audits),
        modes_above_requested_1e_minus9=sum(r['modes_above_requested_1e_minus9'] for r in audits),
        new_eigenvector_bytes=sum(r['eigenvector_bytes'] for r in generated.values()),
        paired_profile_configurations=profiles, paired_profile_count=len(profiles),
        profile_piters=[4,40], target_profile_count=len(completed)-1,
        scope='Prior generation inventory is inherited from the sweep controller. New generation logs and placement are rechecked here. Saved-mode residual/normalization/all-pair Gram checks cover 100 modes on 13 selected times for 16 configurations, using the existing qcd operator; full-time eigenvalue metadata was checked by the archived audit. Gauge and eigenvector HDF5 contents are not included or independently reread here. n100 remains a regression control outside the profile ensemble.')
    return summary, sorted(generated.values(),key=lambda r:r['configuration']), sorted(audits,key=lambda r:r['configuration'])


def write_csv(path, rows):
    path.write_bytes(csv_bytes(rows,list(rows[0])))


def plot(summary, audits, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                         'svg.hashsalt':'qx-eigenpair-completion'})
    def save(fig,name):
        fig.savefig(output/(name+'.png'),dpi=180)
        fig.savefig(output/(name+'.svg'),metadata={'Date':None})
        path=output/(name+'.svg')
        path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
        plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(13,6.4),sharex=True,layout='constrained')
    expected=summary['expected_configurations']
    new=summary['newly_generated_configurations']
    previous=sorted(set(summary['generated_configurations'])-set(new))
    for ax in axes:
        ax.scatter(expected,[0]*len(expected),marker='s',s=35,color='#dce1e5')
        ax.set_ylim(-.65,.65);ax.set_yticks([]);ax.grid(axis='x',alpha=.15)
    axes[0].scatter(previous,[0]*len(previous),marker='s',s=40,color='#5485a4',label='Earlier completed inventory: 57')
    axes[0].scatter(new,[0]*len(new),marker='s',s=40,color='#278260',label='New generation checked: 13')
    axes[0].scatter([640],[0],marker='x',s=85,color='#bb413b',label='n640 unresolved')
    axes[0].set_title('Eigenpair outputs: 70 / 71 configurations; 192 time slices per completed configuration',loc='left',fontsize=12)
    audited=summary['freshly_audited_configurations']
    axes[1].scatter(audited,[0]*len(audited),marker='s',s=40,color='#278260',label='Numerical audit in these archives: 16')
    axes[1].set_title('Saved-mode checks in this review: 100 modes × 13 selected times × 16 configurations',loc='left',fontsize=12)
    profiles=summary['paired_profile_configurations']
    axes[2].scatter(profiles,[0]*len(profiles),marker='s',s=40,color='#9c5f9b',label='Paired piter=4/40 profiles: 13')
    axes[2].scatter([100],[0],marker='o',s=55,facecolors='white',edgecolors='#5485a4',label='n100 regression control')
    axes[2].set_title('Measured profile cohort: n110–230; 69 configurations eligible after excluding n100 and n640',loc='left',fontsize=12)
    for ax in axes:ax.legend(loc='upper left',ncol=3,fontsize=9)
    axes[2].set_xticks(range(100,801,50));axes[2].set_xlabel('Configuration n')
    fig.suptitle('qX generation, saved-mode validation and profile coverage — 13 October 2026',fontsize=14)
    fig.supxlabel('Grey: not covered by that row. Earlier numerical-audit archives are outside this review. Eigenpairs alone do not provide new profiles.',fontsize=10)
    save(fig,'configuration_coverage')
    fig,axes=plt.subplots(2,1,figsize=(12.5,7),sharex=True,layout='constrained')
    x=[r['configuration'] for r in audits]
    for field,label,color,marker in [('max_relative_residual','Maximum over all 13 selected times','#176b8e','o'),
        ('max_central_relative_residual','Maximum at central t/a=95…101','#278260','s')]:
        axes[0].plot(x,[r[field] for r in audits],marker+'-',ms=5,color=color,label=label)
    axes[0].axhline(1e-9,color='#bb413b',ls='--',label='Requested residual: 10⁻⁹')
    axes[0].set_ylabel(r'$\|Av-\lambda v\|_2/(|\lambda|\,\|v\|_2)$')
    axes[0].set_ylim(1e-12,3e-9)
    for field,label,color,marker in [('max_gram_offdiagonal','Maximum off-diagonal Gram entry','#9c5f9b','o'),
        ('max_gram_diagonal_defect','Maximum diagonal Gram defect','#cd5519','^'),
        ('max_norm_defect','Maximum norm defect','#278260','s')]:
        axes[1].plot(x,[r[field] for r in audits],marker+'-',ms=4,color=color,label=label)
    axes[1].axhline(1e-9,color='#bb413b',ls='--',label='Audit tolerance: 10⁻⁹')
    axes[1].set_ylabel('Dimensionless defect');axes[1].set_ylim(2e-16,3e-9)
    for ax in axes:ax.set_yscale('log');ax.grid(axis='y',alpha=.2);ax.legend(ncol=2,fontsize=9)
    axes[1].set_xticks(x);axes[1].tick_params(axis='x',rotation=45);axes[1].set_xlabel('Configuration n')
    fig.suptitle('Saved eigenpair quality: n650–800\n100 modes on 13 selected time slices per configuration',fontsize=14)
    fig.supxlabel('All-pair Gram checks and residuals use the existing qcd Laplacian. Sampled eigenvector validation; all-time eigenvalue metadata checked separately.',fontsize=10)
    save(fig,'saved_mode_quality_16configs')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--no-plots',action='store_true')
    args=parser.parse_args()
    summary,generation,audits=verify_and_reduce()
    output=HERE/'results';output.mkdir(exist_ok=True)
    (output/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    write_csv(output/'generation_13configs.csv',generation)
    write_csv(output/'saved_mode_quality_16configs.csv',audits)
    if not args.no_plots:plot(summary,audits,output)
    print('PASS:',summary['generated_count'],'completed inventories;',len(audits),'fresh audits;',summary['audited_mode_records'],'checked modes; profiles remain',summary['paired_profile_count'])


if __name__=='__main__':main()
