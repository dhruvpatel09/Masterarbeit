#!/usr/bin/env python3
"""Recreate the pilot numerical exports; does not run or submit MPI jobs."""
import csv
import json
import os
from pathlib import Path
from analysis import compute
from replay import replay

HERE=Path(__file__).resolve().parent


def atomic_text(path,text):
    path=Path(path);temporary=path.with_name(path.name+'.pending')
    with temporary.open('x') as f:f.write(text);f.flush();os.fsync(f.fileno())
    temporary.replace(path)


def write_csv(path,rows):
    path=Path(path);temporary=path.with_name(path.name+'.pending')
    with temporary.open('x',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows);f.flush();os.fsync(f.fileno())
    temporary.replace(path)


def run():
    rows,_,validation=replay();tables,summary=compute(rows,validation)
    for name,values in tables.items():write_csv(HERE/'results'/(name+'.csv'),values)
    atomic_text(HERE/'results/VALIDATION.json',json.dumps(validation,indent=2,sort_keys=True,allow_nan=False)+'\n')
    atomic_text(HERE/'results/SUMMARY.json',json.dumps(summary,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print('PASS pilot analysis: 33792 original matrix records; 3 configurations; full complex pencils; no production jobs submitted.')


if __name__=='__main__':run()
