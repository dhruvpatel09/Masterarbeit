#!/usr/bin/env python3
"""Standard-library integrity check for full71 orchestration and native sources."""
from campaign import verify_sources


if __name__=='__main__':
    verify_sources()
    print('PASS full71 matrix package and unchanged native source/reference dependencies')
