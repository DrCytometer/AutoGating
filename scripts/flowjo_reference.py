"""Honeychrome side of the FlowJo import check (scripts/flowjo_check_test.R).

Uses the ag_core / ag_flowjo / ag_review snapshots in
AutoGating-sandbox/dev/reference and Honeychrome's label matcher with
AutoSpectral's fluorophore and marker databases. Run in the Honeychrome
virtual environment:

    python flowjo_reference.py <workspace> [<workspace> ...]

For each workspace (tag = file name without extension) it writes to
scripts/flowjo_reference/:
  <tag>_experiment.csv   the experiment both sides import into: channel
                         (the FlowJo parameter without its matrix prefix),
                         antigen, type
  <tag>_alignment.csv    align_parameters against channels named by the
                         fluorophore alone (no '-A'), antigens from the labels
  <tag>_gates.json       build_gate_definitions with every parameter mapped
                         to its own channel
  <tag>_events.csv       synthetic events (display units) for the review
  <tag>_review.csv       conversion proposals of every imported gate
"""
import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REF = HERE.parent / 'dev' / 'reference'
sys.path.insert(0, str(REF))
import ag_core  # noqa: E402
import ag_flowjo  # noqa: E402
import ag_review  # noqa: E402

AS_DIR = Path(os.environ.get('AUTOSPECTRAL_EXTDATA',
                             HERE.parents[1] / 'AutoSpectral-sandbox' / 'inst' / 'extdata'))
HC_LABELS = Path(os.environ.get('HONEYCHROME_LABELS', HERE.parents[1] / 'honeychrome-sandbox' / 'src'
                                / 'honeychrome' / 'controller_components' / 'label_matching.py'))
OUT = HERE / 'flowjo_reference'
OUT.mkdir(exist_ok=True)

# Display transforms of the experiment: Aurora ranges (22-bit) for Aurora
# workspaces, else 18-bit.
def transforms_for(cytometer):
    top = 4194304 if 'aurora' in (cytometer or '').lower() else 262144
    m = 5.6226599046 if top == 4194304 else 4.418539922
    return ({'id': 1, 'scale_t': top, 'logicle_w': 1, 'logicle_m': m, 'logicle_a': 0},
            {'id': 0, 'scale_t': top, 'linear_a': 0})


def label_matchers():
    import importlib.util
    spec = importlib.util.spec_from_file_location('label_matching', HC_LABELS)
    lm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(lm)
    fdb = lm._load_csv(AS_DIR / 'fluorophore_database.csv')
    mdb = lm._load_csv(AS_DIR / 'marker_database.csv')
    return (lambda n: lm.match_fluorophore(n, fdb)), (lambda n: lm.match_marker(n, mdb))


def is_scatter(stem):
    return stem.split('-')[0] in ('FSC', 'SSC')


def experiment(params, affixes, cytometer):
    stems = []
    for p in params:
        s = ag_flowjo.parameter_stem(p, affixes)
        if s not in stems:
            stems.append(s)
    scatter = [s for s in stems if is_scatter(s)]
    fluor = [s for s in stems if not is_scatter(s)]
    fluor_t, scatter_t = transforms_for(cytometer)
    transforms = {s: dict(scatter_t if s in scatter else fluor_t) for s in stems}
    pair = ('FSC-A', 'FSC-H') if {'FSC-A', 'FSC-H'} <= set(stems) else None
    return stems, ag_flowjo.ExperimentChannels(
        channels=stems, fluorescence=fluor, scatter=scatter,
        antigens={s: s for s in stems}, transforms=transforms, singlet_pair=pair)


def jsonable(obj):
    return ag_core.to_jsonable(obj)


def write_csv(path, header, rows):
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def synthetic_events(exp, n=30000, seed=7):
    """Scatter around a cell cloud, each fluorescence channel a negative
    population at its raw zero plus a positive one."""
    rng = np.random.default_rng(seed)
    cols = []
    for ch in exp.channels:
        if ch in exp.scatter:
            centre = 0.45 if ch.startswith('FSC') else 0.25
            v = rng.normal(centre, 0.08, n)
            if ch.endswith('-H'):
                v = 0.8 * v + rng.normal(0, 0.01, n)
        else:
            zero = float(ag_core.make_transform(exp.transforms[ch]).apply(np.zeros(1))[0])
            pos = rng.random(n) < rng.uniform(0.1, 0.5)
            v = np.where(pos, rng.normal(zero + 0.35, 0.05, n), rng.normal(zero, 0.04, n))
        cols.append(v)
    return np.column_stack(cols)


def review(gate_defs, data, exp):
    index = {c: i for i, c in enumerate(exp.channels)}
    result = ag_core.run_gating(gate_defs, data, index, reference=None, transforms=exp.transforms)
    zeros = {ch: float(ag_core.make_transform(t).apply(np.zeros(1))[0]) for ch, t in exp.transforms.items()}
    anchor = ag_core.anchor_population(gate_defs)
    rows = []
    for g in ag_core.ordered_gate_defs(gate_defs):
        if g.get('algorithm') != ag_core.ALGORITHM_IMPORTED:
            continue
        ev = ag_review.conversion_evidence(g, gate_defs, data, index, result, None, zeros, anchor)
        d = ag_review.proposed_definition(g, None, None, ev)
        gp = d.get('gate_param') or {}
        if d.get('gate_type') == '2dsep':
            ax = (gp.get('algorithm_x') or d['algorithm'], gp.get('algorithm_y') or d['algorithm'])
        else:
            ax = (d.get('algorithm'), '')
        fit = ev.get('fit') or {}
        rows.append([g['gate_name'], g.get('gate_type'), d.get('gate_type'), ax[0], ax[1],
                     ev.get('drop_axis') or '', ev.get('n', ''),
                     (ev.get('below_zero') or {}).get('x', ''), (ev.get('below_zero') or {}).get('y', ''),
                     (ev.get('positive') or {}).get('x', ''), (ev.get('positive') or {}).get('y', ''),
                     max(fit.values()) if fit else '',
                     json.dumps(jsonable({k: v for k, v in gp.items()}), sort_keys=True)])
    return rows


def main(paths):
    fluor_match, marker_match = label_matchers()
    for path in paths:
        tag = Path(path).stem.replace(' ', '_')
        ws = ag_flowjo.read_workspace(path)
        src = ws.default_source()
        params = ag_flowjo.gated_parameters(src)
        stems, exp = experiment(params, ws.affixes, ws.cytometer)
        keys = ['id', 'scale_t', 'linear_a', 'logicle_w', 'logicle_m', 'logicle_a']
        write_csv(OUT / f'{tag}_experiment.csv', ['channel', 'antigen', 'type'] + keys,
                  [[s, s, 'scatter' if s in exp.scatter else 'lineage']
                   + [exp.transforms[s].get(k, '') for k in keys] for s in stems])

        named = {s: (s[:-2] if s.endswith('-A') and s not in exp.scatter else s) for s in stems}
        exp_align = ag_flowjo.ExperimentChannels(
            channels=[named[s] for s in stems], fluorescence=[named[s] for s in exp.fluorescence],
            scatter=list(exp.scatter),
            antigens={named[s]: ((src.labels.get(p) or '').split('_')[0] or named[s])
                      for p in params for s in [ag_flowjo.parameter_stem(p, ws.affixes)]},
            transforms={}, singlet_pair=exp.singlet_pair)
        al = ag_flowjo.align_parameters(params, src.labels, exp_align, affixes=ws.affixes,
                                        fluor_match=fluor_match, marker_match=marker_match)
        write_csv(OUT / f'{tag}_alignment.csv',
                  ['parameter', 'channel', 'antigen', 'method', 'label'],
                  [[p, a['channel'] or '', exp_align.antigens.get(a['channel'], '') if a['channel'] else '',
                    a['method'] or '', a['label']] for p, a in al.items()])

        mapping = {p: ag_flowjo.parameter_stem(p, ws.affixes) for p in params}
        gate_defs, report = ag_flowjo.build_gate_definitions(src, mapping, exp,
                                                              workspace_name=Path(path).name)
        with open(OUT / f'{tag}_gates.json', 'w') as fh:
            json.dump(jsonable({'gate_definitions': gate_defs,
                                'dropped': report['plan'].messages(),
                                'populations': report['populations']}), fh, indent=1)

        data = synthetic_events(exp)
        write_csv(OUT / f'{tag}_events.csv', exp.channels, [[repr(float(v)) for v in r] for r in data])
        write_csv(OUT / f'{tag}_review.csv',
                  ['gate', 'imported_type', 'type', 'algorithm_x', 'algorithm_y', 'drop_axis', 'n',
                   'below_zero_x', 'below_zero_y', 'positive_x', 'positive_y', 'agreement',
                   'gate_param'],
                  review(gate_defs, data, exp))
        print(f"{tag}: {len(gate_defs)} gates")


if __name__ == '__main__':
    main(sys.argv[1:])
