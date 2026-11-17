"""Honeychrome side of the .agmodel round-trip check (scripts/agmodel_check_test.R).

Uses the ag_core.py snapshot in AutoGating-sandbox/dev/reference and the
synthetic events written by engine_reference.py. Run in the Honeychrome
virtual environment:

    python agmodel_reference.py make
        trains the engine_reference gates on sample A, saves the model as
        engine_reference/model_py.agmodel and applies it to sample B
        (engine_reference/counts_py.csv)

    python agmodel_reference.py apply engine_reference/model_r.agmodel
        opens a model written by AutoGating as Honeychrome does and applies
        it to sample B (engine_reference/counts_r_model.csv)
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'dev' / 'reference'))
import ag_core  # noqa: E402

REF = HERE / 'engine_reference'
CHANNELS = ['FSC-A', 'FSC-H', 'SSC-A', 'Viab', 'CD3', 'CD4', 'CD8']
INDEX = {c: i for i, c in enumerate(CHANNELS)}
TRANSFORMS = {'CD3': {'id': 1, 'scale_t': 262144, 'logicle_w': 0.5, 'logicle_m': 4.5,
                      'logicle_a': 0}}
RECALCULATE = ('live', 'cd4cd8')


def read_events(name):
    with open(REF / name, newline='') as fh:
        rows = list(csv.reader(fh))
    assert rows[0] == CHANNELS
    return np.array(rows[1:], dtype=float)


def apply_model(model, data):
    app = model['application']
    return ag_core.run_gating(model['gate_definitions'], data, INDEX,
                              reference=model['trained_boundaries'],
                              modes=app['mode_per_gate'],
                              default_mode=app['default_mode'],
                              drift_limit=app['drift_limit'],
                              drift_action=app['drift_action'],
                              transforms=model['transforms'])


def write_counts(result, name):
    with open(REF / name, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['gate', 'population', 'n_events', 'threshold_x', 'threshold_y', 'source', 'flag'])
        for gate, pops in result.boundaries.items():
            st = result.status.get(gate, {})
            for p, e in pops.items():
                w.writerow([gate, p, int(e.get('n_events', -1)),
                            '' if e.get('threshold_x') is None else repr(float(e['threshold_x'])),
                            '' if e.get('threshold_y') is None else repr(float(e['threshold_y'])),
                            st.get('source') or '', st.get('flag') or ''])


def make():
    with open(REF / 'gates.json') as fh:
        gate_defs = json.load(fh)
    data_a = read_events('data_a.csv')
    train = ag_core.train_boundaries(gate_defs, data_a, INDEX, transforms=TRANSFORMS)
    modes = {g['gate_name']: ('recalculate' if g['gate_name'] in RECALCULATE else 'fixed')
             for g in gate_defs}
    model = ag_core.new_model(
        gate_defs, ag_core.serialisable_boundaries(train.boundaries),
        channel_map={c: c for c in CHANNELS}, transforms=TRANSFORMS,
        training_summary={'training_sample_names': ['A'], 'events_cap_per_sample': 50000},
        metadata={'created': 'reference'},
        application={'mode_per_gate': modes})
    ag_core.write_model(REF / 'model_py.agmodel', model)
    model, warnings = ag_core.read_model(REF / 'model_py.agmodel')
    write_counts(apply_model(model, read_events('data_b.csv')), 'counts_py.csv')
    print(f"Model and counts written to {REF}")


def apply(path):
    model, warnings = ag_core.read_model(path)
    for w in warnings:
        print('warning:', w)
    write_counts(apply_model(model, read_events('data_b.csv')), 'counts_r_model.csv')
    print(f"Counts written to {REF / 'counts_r_model.csv'}")


if __name__ == '__main__':
    if sys.argv[1:2] == ['make']:
        make()
    elif sys.argv[1:2] == ['apply'] and len(sys.argv) == 3:
        apply(sys.argv[2])
    else:
        sys.exit(__doc__)
