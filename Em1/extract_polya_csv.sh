#!/usr/bin/env bash
set -euo pipefail

BDIO_FILE="${1:-dat/Em1p4.polya.bdio}"
OUT_FILE="${2:-dat/Em1p4.polya.csv}"

echo "nc,corr_id,ReP0,ImP0,AbsP0,ReP1,ImP1,AbsP1,ReP2,ImP2,AbsP2,ReP3,ImP3,AbsP3" > "$OUT_FILE"

# Binary records are 4, 6, 8, ..., 62 for nc = 1,...,30.
# More generally, after metadata records 1 and 2, each config has:
# ASCII record, then binary f64 record.
for rec in $(seq 4 2 62); do
    vals=$(lsbdio -c 0 -d "$rec" "$BDIO_FILE")

    echo "$vals" | awk '
    {
        a[NR]=$1
    }
    END {
        nc=a[1]
        corr_id=a[2]

        re0=a[3];  im0=a[4]
        re1=a[5];  im1=a[6]
        re2=a[7];  im2=a[8]
        re3=a[9];  im3=a[10]

        abs0=sqrt(re0*re0 + im0*im0)
        abs1=sqrt(re1*re1 + im1*im1)
        abs2=sqrt(re2*re2 + im2*im2)
        abs3=sqrt(re3*re3 + im3*im3)

        printf "%.0f,%.0f,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n", \
               nc,corr_id,re0,im0,abs0,re1,im1,abs1,re2,im2,abs2,re3,im3,abs3
    }' >> "$OUT_FILE"
done

echo "Wrote: $OUT_FILE"