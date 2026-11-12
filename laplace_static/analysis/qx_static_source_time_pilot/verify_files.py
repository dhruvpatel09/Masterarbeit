"""SHA-256 integrity check; standard-library only."""
from common import HERE, digest, read


def verify():
    manifest = read(HERE / 'FILE_MANIFEST.json')
    for name, record in manifest['files'].items():
        path = HERE / name
        if not path.is_file() or path.stat().st_size != record['bytes'] or digest(path) != record['sha256']:
            raise ValueError('Source-time analysis package changed: ' + name)
    print('PASS source-time pilot package: %d pinned files, original review and published primary cache.' % len(manifest['files']))


if __name__ == '__main__':
    verify()
