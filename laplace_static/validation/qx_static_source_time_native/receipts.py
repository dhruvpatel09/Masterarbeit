"""Atomic receipts and deterministic table exports; no numerical dependencies."""
import csv,hashlib,json,os
from pathlib import Path


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):h.update(block)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())


def write(path,value,exclusive=False):
    path=Path(path)
    if exclusive and path.exists():raise ValueError('Receipt already exists: '+str(path))
    temporary=path.with_name(path.name+'.pending')
    with temporary.open('x') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n');stream.flush();os.fsync(stream.fileno())
    temporary.replace(path)


def csv_write(path,values):
    path=Path(path);temporary=path.with_name(path.name+'.pending')
    with temporary.open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(values[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(values);stream.flush();os.fsync(stream.fileno())
    temporary.replace(path)


def csv_equal(path,values):
    with Path(path).open(newline='') as stream:old=list(csv.DictReader(stream))
    if old!=[{key:str(value) for key,value in row.items()} for row in values]:
        raise ValueError('Validated CSV differs from native log replay: '+str(path))
