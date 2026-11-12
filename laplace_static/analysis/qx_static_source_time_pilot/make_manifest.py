"""Intentional receipt generation after producing and validating a release."""
from common import HERE, digest, write_json

files = {}
for path in sorted(HERE.rglob('*')):
    if (not path.is_file() or '__pycache__' in path.parts or path.name == 'FILE_MANIFEST.json'
            or path.suffix == '.pyc' or 'work' in path.relative_to(HERE).parts):
        continue
    files[path.relative_to(HERE).as_posix()] = dict(bytes=path.stat().st_size, sha256=digest(path))
write_json(HERE / 'FILE_MANIFEST.json', dict(files=files))
print('Wrote package integrity receipt for %d files' % len(files))
