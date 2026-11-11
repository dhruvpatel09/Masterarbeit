"""Balanced held-out splits with complete source retraining in every deletion."""
import numpy as np
from common import HERE,PREVIOUS,RADII,CONFIGS,csv_rows,previous_projections,close


def make_assignments(plan):
    labels=[np.arange(71)%2]
    for seed in plan['split_seeds']:
        permutation=np.random.RandomState(seed).permutation(71)
        fold=np.ones(71,dtype=int);fold[permutation[:36]]=0;labels.append(fold)
    return np.array(labels,dtype=np.int8)


def previous_keys():
    rows=csv_rows(PREVIOUS/'results/projected_correlators.csv')
    if len(rows)!=72*33:raise ValueError('Published curve inventory')
    keys=[]
    for i in range(72):
        group=rows[i*33:(i+1)*33];r=group[0]
        key=(int(r['R']),r['basis'],int(r['anchor']),r['procedure'])
        if [int(a['tau']) for a in group]!=list(range(33)):raise ValueError('Published time inventory')
        keys.append(key)
    return keys


def compute_splits(C,plan,assignments,progress=False):
    native=previous_projections();cache=np.load(str(PREVIOUS/'results/projection_resamples.npz'),allow_pickle=False)
    oldkeys=previous_keys();keys=[];means=[];deletions=[];membership=[]
    for name,fold in zip(plan['split_names'],assignments):
        if sorted(np.bincount(fold).tolist())!=[35,36]:raise ValueError('Unbalanced split')
        for cfg,label in zip(CONFIGS,fold):membership.append(dict(split=name,configuration=cfg,heldout_fold=int(label)))
    ids=np.arange(71)
    for ri,r in enumerate(RADII):
        for basis,indices in native.BASES.items():
            selected=native.submatrix(C[:,ri],indices)
            for anchor in plan['training_times']:
                key=dict(R=r,basis=basis,anchor=anchor,variant='same_sample')
                old=oldkeys.index((r,basis,anchor,'same_sample'))
                central=native.estimate(selected,ids,anchor,'same_sample')
                close(central.real,cache['projected_means'][old].real,'same-sample central')
                close(central.imag,cache['projected_means'][old].imag,'same-sample central imaginary')
                keys.append(key);means.append(cache['projected_means'][old]);deletions.append(cache['projected_jackknife_means'][old])
                for si,(name,fold) in enumerate(zip(plan['split_names'],assignments)):
                    mean=native.estimate(selected,fold,anchor,'cross_fit')
                    reps=[]
                    for deleted in ids:
                        keep=ids!=deleted
                        reps.append(native.estimate(selected[keep],fold[keep],anchor,'cross_fit'))
                    reps=np.array(reps)
                    if si==0:
                        old=oldkeys.index((r,basis,anchor,'cross_fit'))
                        close(mean.real,cache['projected_means'][old].real,'published alternating central')
                        close(mean.imag,cache['projected_means'][old].imag,'published alternating imaginary')
                        close(reps.real,cache['projected_jackknife_means'][old].real,'published alternating deletions')
                        close(reps.imag,cache['projected_jackknife_means'][old].imag,'published alternating imaginary deletions')
                    keys.append(dict(R=r,basis=basis,anchor=anchor,variant=name));means.append(mean);deletions.append(reps)
        if progress:print('Reproduced alternating split and four new held-out splits: R={}'.format(r),flush=True)
    return keys,np.asarray(means),np.asarray(deletions),membership
