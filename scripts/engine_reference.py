"""Reference results of the Honeychrome gating engine (ag_core), for checking
the R port in AutoGating (scripts/engine_check_test.R).

Uses the ag_core.py snapshot in AutoGating-sandbox/dev/reference. Run in the
Honeychrome virtual environment:
    python ~/dev/claude-sandboxes/AutoGating-sandbox/scripts/engine_reference.py
Writes CSV and JSON files to scripts/engine_reference/.
"""

import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'dev' / 'reference'))
import ag_core  # noqa: E402

OUT = HERE / 'engine_reference'
OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(20261111)


def mix(parts):
    return np.concatenate([rng.normal(m, s, n) for n, m, s in parts])


def num(v):
    return '' if v is None else repr(float(v))


def write_csv(name, header, rows):
    with open(OUT / name, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


# ---- 1. Threshold algorithms on vectors --------------------------------------

vectors = {
    'bimodal_small': mix([(3000, .25, .05), (1000, .6, .06)]),
    'bimodal_large': mix([(45000, .25, .05), (15000, .6, .06)]),
    'unimodal_tail': mix([(10000, .3, .04), (300, .55, .05)]),
    'trimodal': mix([(4000, .2, .04), (3000, .45, .04), (2000, .7, .05)]),
    'dbic_unimodal': mix([(3000, .3, .05)]),
    'debris_fsc': mix([(2000, .08, .03), (8000, .45, .06)]),
}
write_csv('vectors.csv', ['vector', 'value'],
          [[k, repr(float(v))] for k, vs in vectors.items() for v in vs])

rows = []
for name, v in vectors.items():
    res = {
        'tail_1e-4': ag_core.compute_threshold(v, 'tail', {'tail_fraction': 1e-4}),
        'tail_1e-2': ag_core.compute_threshold(v, 'tail', {'tail_fraction': 1e-2}),
        'kde_min': ag_core.compute_threshold(v, 'kde_min'),
        'bimodal': ag_core.compute_threshold(v, 'bimodal'),
        'otsu': ag_core.compute_threshold(v, 'otsu'),
        'mixture2': ag_core.compute_threshold(v, 'mixture'),
        'mixture3_outer': ag_core.compute_threshold(v, 'mixture', {'mixture_components': 3}),
        'mixture3_lower': ag_core.compute_threshold(
            v, 'mixture', {'mixture_components': 3, 'mixture_cut': 'lower'}),
        'mixture3_upper': ag_core.compute_threshold(
            v, 'mixture', {'mixture_components': 3, 'mixture_cut': 'upper'}),
        'kde_mode': ag_core.kde_mode(v),
        'delta_bic': ag_core.bimodality_delta_bic(v),
        'debris_cut': ag_core.debris_cut(v, {}),
        'extent_lo': ag_core.density_edges(v)[0],
        'extent_hi': ag_core.density_edges(v)[1],
    }
    rows += [[name, k, num(val)] for k, val in res.items()]
write_csv('thresholds.csv', ['vector', 'quantity', 'value'], rows)


# ---- 2. A gate hierarchy, trained on one sample and applied to another -------

CHANNELS = ['FSC-A', 'FSC-H', 'SSC-A', 'Viab', 'CD3', 'CD4', 'CD8']


def make_events(n, seed, cd3_pos_mean):
    r = np.random.default_rng(seed)
    n_debris = int(0.1 * n)
    n_cells = n - n_debris
    fsc = np.concatenate([r.normal(0.08, 0.03, n_debris), r.normal(0.45, 0.06, n_cells)])
    doublet = r.random(n) < 0.05
    fsch = 0.8 * fsc + r.normal(0, 0.01, n)
    fsch[doublet] = 0.55 * fsc[doublet] + r.normal(0, 0.01, doublet.sum())
    ssc = np.concatenate([r.normal(0.1, 0.04, n_debris), r.normal(0.3, 0.07, n_cells)])
    dead = r.random(n) < 0.15
    viab = np.where(dead, r.normal(0.7, 0.05, n), r.normal(0.2, 0.04, n))
    t = r.random(n) < 0.35
    cd3 = np.where(t, r.normal(cd3_pos_mean, 0.05, n), r.normal(0.25, 0.04, n))
    u = r.random(n)
    cd4p = (t & (u < 0.6)) | (t & (u >= 0.9) & (u < 0.95))
    cd8p = t & (u >= 0.6) & (u < 0.95)
    cd4 = np.where(cd4p, r.normal(0.6, 0.05, n), r.normal(0.22, 0.04, n))
    cd8 = np.where(cd8p, r.normal(0.62, 0.05, n), r.normal(0.22, 0.04, n))
    return np.column_stack([fsc, fsch, ssc, viab, cd3, cd4, cd8])


data_a = make_events(20000, 1, 0.65)
data_b = make_events(20000, 2, 0.75)
for name, d in (('data_a.csv', data_a), ('data_b.csv', data_b)):
    write_csv(name, CHANNELS, [[repr(float(v)) for v in row] for row in d])
index = {c: i for i, c in enumerate(CHANNELS)}


def pop(label, region, **extra):
    return dict(label=label, label_pos=1, region=region, **extra)


gate_defs = [
    dict(gate_name='singlets', gate_number=1, gate_type='singlets', gate_marker_x='FSC-A',
         gate_marker_y='FSC-H', parent_gate=None, parent_popul=None,
         populations={'singlets': pop('Singlets', 'inside')}, algorithm='singlets',
         gate_param={}),
    dict(gate_name='cells', gate_number=2, gate_type='free', gate_marker_x='FSC-A',
         gate_marker_y='SSC-A', parent_gate='singlets', parent_popul='singlets',
         populations={'cells': pop('Cells', 'inside')}, algorithm='cells', gate_param={}),
    dict(gate_name='live', gate_number=3, gate_type='1dsep', gate_marker_x='Viab',
         gate_marker_y=None, parent_gate='cells', parent_popul='cells',
         populations={'live': pop('Live', 'below'), 'dead': pop('Dead', 'above')},
         algorithm='mixture', gate_param={}),
    dict(gate_name='cd3', gate_number=4, gate_type='1dsep', gate_marker_x='CD3',
         gate_marker_y='SSC-A', sep_axis='x', parent_gate='live', parent_popul='live',
         populations={'T': pop('T cells', 'above'), 'nonT': pop('Non-T', 'below')},
         algorithm='tail', gate_param={'tail_fraction': 0.001, 'threshold_factor': 1.5}),
    dict(gate_name='cd4cd8', gate_number=5, gate_type='2dsep', gate_marker_x='CD4',
         gate_marker_y='CD8', parent_gate='cd3', parent_popul='T',
         populations={'CD4': pop('CD4', 'x+y-'), 'CD8': pop('CD8', 'x-y+'),
                      'DP': pop('DP', 'x+y+'), 'DN': pop('DN', 'x-y-')},
         algorithm='mixture', gate_param={'algorithm_y': 'tail'},
         stats_parent={'CD4': 'live/live'}),
    dict(gate_name='notcd4', gate_number=6, gate_type='2dsep', gate_marker_x='CD4',
         gate_marker_y='CD8', parent_gate='cd3', parent_popul='T',
         populations={'nonCD4': pop('Non-CD4', 'x-y-|x-y+|x+y+')}, algorithm='replicate',
         gate_param={'replicate_x': {'gate': 'cd4cd8'}, 'replicate_y': {'gate': 'cd4cd8'}}),
    dict(gate_name='cd4ref', gate_number=7, gate_type='1dsep', gate_marker_x='CD4',
         gate_marker_y=None, parent_gate='cd3', parent_popul='nonT',
         populations={'CD4lo': pop('CD4 low', 'below'), 'CD4hi': pop('CD4 high', 'above')},
         algorithm='replicate',
         gate_param={'replicate_x': {'population': 'live/live', 'algorithm': 'kde_min'}}),
    dict(gate_name='cd8otsu', gate_number=8, gate_type='1dsep', gate_marker_x='CD8',
         gate_marker_y=None, parent_gate='cd3', parent_popul='T',
         populations={'lo': pop('lo', 'below'), 'hi': pop('hi', 'above')},
         algorithm='otsu', gate_param={}),
    dict(gate_name='cd8rep', gate_number=9, gate_type='replicate', origin_gate='cd8otsu',
         parent_gate='live', parent_popul='live', populations={}, gate_param={}),
    dict(gate_name='hull', gate_number=10, gate_type='free', gate_marker_x='CD3',
         gate_marker_y='CD4', parent_gate='cd3', parent_popul='T',
         populations={'in': pop('All T', 'inside')}, algorithm='all', gate_param={}),
]
with open(OUT / 'gates.json', 'w') as fh:
    json.dump(gate_defs, fh, indent=1)

transforms = {'CD3': {'id': 1, 'scale_t': 262144, 'logicle_w': 0.5, 'logicle_m': 4.5,
                      'logicle_a': 0}}

train = ag_core.train_boundaries(gate_defs, data_a, index, transforms=transforms)
applied = ag_core.run_gating(gate_defs, data_b, index, reference=train.boundaries,
                             modes={'cd3': 'recalculate', 'live': 'recalculate'},
                             transforms=transforms)


def result_rows(run, result):
    out = []
    for gate, pops in result.boundaries.items():
        st = result.status.get(gate, {})
        for p, e in pops.items():
            b = np.asarray(e.get('boundary') or [[np.nan, np.nan]], dtype=float)
            out.append([run, gate, p, int(e.get('n_events', -1)), num(e.get('threshold_x')),
                        num(e.get('threshold_y')), num(np.nanmin(b[:, 0])), num(np.nanmax(b[:, 0])),
                        num(np.nanmin(b[:, 1])), num(np.nanmax(b[:, 1])), st.get('source') or '',
                        st.get('flag') or '', num(st.get('drift'))])
    for gate, st in result.status.items():
        if gate not in result.boundaries:
            out.append([run, gate, '', -1, '', '', '', '', '', '', st.get('source') or '',
                        st.get('flag') or '', ''])
    return out


write_csv('gating.csv', ['run', 'gate', 'population', 'n_events', 'threshold_x', 'threshold_y',
                         'x_lo', 'x_hi', 'y_lo', 'y_hi', 'source', 'flag', 'drift'],
          result_rows('train', train) + result_rows('apply', applied))

table = ag_core.population_table(gate_defs, applied, n_total=len(data_b))
write_csv('population_table.csv',
          ['key', 'n_events', 'n_parent', 'fraction_of_parent', 'n_stats_parent',
           'fraction_of_stats_parent'],
          [[r['key'], r['n_events'], r['n_parent'], num(r['fraction_of_parent']),
            r.get('n_stats_parent', ''), num(r.get('fraction_of_stats_parent'))] for r in table])

with open(OUT / 'anchor.txt', 'w') as fh:
    fh.write(str(ag_core.anchor_population(gate_defs)) + '\n')
print(f"Reference written to {OUT}")
