"""Reference values of the Gating-ML transforms from FlowKit, for checking
the R implementation in AutoGating (scripts/data_layer_check_test.R).

Run in the Honeychrome virtual environment:
    python ~/dev/claude-sandboxes/AutoGating-sandbox/scripts/transform_reference.py
Writes transform_reference.csv next to this script.
"""

import csv
from pathlib import Path

import numpy as np
from flowkit import transforms as fk

X = np.array([-1e5, -1e4, -1000, -100, -10, -1, -0.01, 0, 0.01, 1, 10, 100,
              1000, 5000, 1e4, 5e4, 1e5, 262144, 1e6, 4194304], dtype=float)

LOGICLE = [
    (262144, 1.0, 4.418539922, 0.0),
    (4194304, 1.5, 5.62, 0.0),
    (262144, 0.5, 4.5, 0.0),
    (10000, 0.5, 4.5, 1.0),
    (262144, 0.0, 4.5, 0.0),
    (1000000, np.log10(250) / 2, 5.0, 0.0),
]
LINEAR = [(262144, 0.0), (4194304, 100.0)]
LOG = [(262144, 4.5)]

rows = []
for t, w, m, a in LOGICLE:
    y = fk.LogicleTransform(param_t=t, param_w=w, param_m=m, param_a=a).apply(X)
    rows += [dict(id=1, scale_t=t, linear_a='', logicle_w=w, logicle_m=m, logicle_a=a,
                  log_m='', x=x, y=v) for x, v in zip(X, np.asarray(y, dtype=float))]
for t, a in LINEAR:
    y = fk.LinearTransform(param_t=t, param_a=a).apply(X)
    rows += [dict(id=0, scale_t=t, linear_a=a, logicle_w='', logicle_m='', logicle_a='',
                  log_m='', x=x, y=v) for x, v in zip(X, np.asarray(y, dtype=float))]
pos = X[X > 0]
for t, m in LOG:
    y = fk.LogTransform(param_t=t, param_m=m).apply(pos)
    rows += [dict(id=2, scale_t=t, linear_a='', logicle_w='', logicle_m='', logicle_a='',
                  log_m=m, x=x, y=v) for x, v in zip(pos, np.asarray(y, dtype=float))]

out = Path(__file__).with_name('transform_reference.csv')
with open(out, 'w', newline='') as fh:
    writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
    writer.writeheader()
    for r in rows:
        writer.writerow({k: (repr(float(v)) if isinstance(v, (float, np.floating)) else v)
                         for k, v in r.items()})
print(f"{len(rows)} rows written to {out}")
