"""Honeychrome side of the differential statistics check
(scripts/stats_check_test.R).

Uses the snapshots of ag_results and honeychrome.controller_components
(differential_stats, moderated_stats, anova_stats) in
AutoGating-sandbox/dev/reference. Run in the Honeychrome virtual environment:

    python stats_reference.py

Writes to scripts/stats_reference/: a synthetic experiment (gates.json,
population_table.csv and medians.csv in AutoGating's run.autogating
format, samples.csv) and, for each design, freq_<design>.csv and
markers_<design>.csv from ag_results.run_statistics.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'dev' / 'reference'))
import ag_results  # noqa: E402

OUT = HERE / 'stats_reference'
OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(11)


def pop(label, region):
    return {'label': label, 'label_pos': 1, 'region': region}


gate_defs = [
    dict(gate_name='cells', gate_number=1, gate_type='free', gate_marker_x='FSC-A',
         gate_marker_y='SSC-A', parent_gate=None, parent_popul='root',
         populations={'cells': pop('Cells', 'inside')}, algorithm='cells', gate_param={}),
    dict(gate_name='cd3', gate_number=2, gate_type='1dsep', gate_marker_x='CD3', gate_marker_y=None,
         parent_gate='cells', parent_popul='cells',
         populations={'T': pop('T cells', 'above'), 'nonT': pop('Non-T', 'below')},
         algorithm='mixture', gate_param={}),
    dict(gate_name='cd4cd8', gate_number=3, gate_type='2dsep', gate_marker_x='CD4',
         gate_marker_y='CD8', parent_gate='cd3', parent_popul='T',
         populations={'CD4': pop('CD4', 'x+y-'), 'CD8': pop('CD8', 'x-y+'), 'DN': pop('DN', 'x-y-'),
                      'DP': pop('DP', 'x+y+')},
         algorithm='mixture', gate_param={}, stats_parent={'DP': 'cells/cells'}),
    dict(gate_name='cd69', gate_number=4, gate_type='1dsep', gate_marker_x='CD69', gate_marker_y=None,
         parent_gate='cd4cd8', parent_popul='CD4',
         populations={'act': pop('CD69+', 'above'), 'rest': pop('CD69-', 'below')},
         algorithm='tail', gate_param={}),
    dict(gate_name='rare', gate_number=5, gate_type='1dsep', gate_marker_x='CD25', gate_marker_y=None,
         parent_gate='cd4cd8', parent_popul='DP',
         populations={'hi': pop('CD25 hi', 'above'), 'lo': pop('CD25 lo', 'below')},
         algorithm='tail', gate_param={}),
]
groups = ['A'] * 5 + ['B'] * 5 + ['C'] * 5
samples = [f's{i + 1:02d}' for i in range(15)]
donor = [f'd{(i % 5) + 1}' for i in range(15)]                    # crossed with group
mouse = [f'{g}{k}' for g, k in zip(groups, [1, 1, 2, 2, 3] * 3)]  # nested in group
age = [float(a) for a in rng.uniform(20, 70, 15).round(1)]
effect = {'A': 0.0, 'B': 0.6, 'C': -0.4}
channels = ['CD3', 'CD4', 'CD8', 'CD69', 'CD25', 'Ki67']

summaries, rows, med_rows = {}, [], []
for s, g, d in zip(samples, groups, donor):
    shift = effect[g] + 0.15 * int(d[1:])
    n_total = int(rng.integers(30000, 60000))
    cells = int(n_total * rng.uniform(0.7, 0.8))
    t = int(rng.binomial(cells, 1 / (1 + np.exp(-(0.2 + 0.3 * shift)))))
    cd4 = int(rng.binomial(t, 0.6 + 0.05 * effect[g]))
    cd8 = int(rng.binomial(t - cd4, 0.8))
    dp = int(rng.binomial(t - cd4 - cd8, 0.1))
    dn = t - cd4 - cd8 - dp
    act = int(rng.binomial(cd4, 1 / (1 + np.exp(-(-2 + shift)))))
    hi = int(rng.binomial(dp, 0.02 if g != 'C' else 0.0))
    counts = {'cells/cells': (cells, n_total, 'root'), 'cd3/T': (t, cells, 'cells/cells'),
              'cd3/nonT': (cells - t, cells, 'cells/cells'), 'cd4cd8/CD4': (cd4, t, 'cd3/T'),
              'cd4cd8/CD8': (cd8, t, 'cd3/T'), 'cd4cd8/DN': (dn, t, 'cd3/T'),
              'cd4cd8/DP': (dp, t, 'cd3/T'), 'cd69/act': (act, cd4, 'cd4cd8/CD4'),
              'cd69/rest': (cd4 - act, cd4, 'cd4cd8/CD4'), 'rare/hi': (hi, dp, 'cd4cd8/DP'),
              'rare/lo': (dp - hi, dp, 'cd4cd8/DP')}
    pops, markers = [], {}
    for key, (n, parent, pkey) in counts.items():
        sp = cells if key == 'cd4cd8/DP' else None
        pops.append({'key': key, 'n_events': n, 'n_parent': parent, 'parent_key': pkey,
                     'n_stats_parent': sp, 'stats_parent_missing': False})
        gate, popn = key.split('/')
        rows.append([s, key, gate, popn, key, n, pkey, parent, n / parent if parent else '',
                     'cells/cells' if sp else '', sp if sp else '',
                     (n / sp) if sp else '', 'FALSE'])
        if n >= 1:
            markers[key] = {}
            for j, ch in enumerate(channels):
                m = 0.3 + 0.05 * j + 0.1 * effect[g] * (j % 3) + rng.normal(0, 0.01)
                markers[key][ch] = [m, m - 0.05, m + 0.05]
                med_rows.append([s, key, ch, repr(m), repr(m)])
    summaries[s] = {'populations': pops, 'markers': markers, 'n_total': n_total}

with open(OUT / 'gates.json', 'w') as fh:
    json.dump(gate_defs, fh, indent=1)


def write(name, header, data):
    with open(OUT / name, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(data)


write('population_table.csv', ['sample', 'key', 'gate', 'population', 'label', 'n.events', 'parent.key',
                               'n.parent', 'fraction.of.parent', 'stats.parent', 'n.stats.parent',
                               'fraction.of.stats.parent', 'stats.parent.missing'], rows)
write('medians.csv', ['sample', 'key', 'channel', 'median', 'median.raw'], med_rows)
write('samples.csv', ['sample', 'group', 'donor', 'mouse', 'age'],
      [[s, g, d, m, a] for s, g, d, m, a in zip(samples, groups, donor, mouse, age)])

designs = {
    'reference': dict(contrast_mode='reference', reference='A'),
    'pairwise': dict(contrast_mode='pairwise'),
    'fixed': dict(contrast_mode='reference', reference='A', pairing=dict(zip(samples, donor)),
                  pairing_mode='auto'),
    'random': dict(contrast_mode='reference', reference='A', pairing=dict(zip(samples, mouse)),
                   pairing_mode='auto'),
    'covariate': dict(contrast_mode='reference', reference='A',
                      covariates=pd.DataFrame({'age': age}, index=samples),
                      numeric_covariates=['age']),
}
for name, kw in designs.items():
    res = ag_results.run_statistics(summaries, samples, dict(zip(samples, groups)), gate_defs,
                                    group_order=['A', 'B', 'C'], **kw)
    res['freq'].to_csv(OUT / f'freq_{name}.csv', index=False)
    res['markers'].to_csv(OUT / f'markers_{name}.csv', index=False)
    print(name, res['pairing_mode'], len(res['freq']), len(res['markers']),
          int(res['freq']['significant'].sum()), int(res['markers']['significant'].sum()))
