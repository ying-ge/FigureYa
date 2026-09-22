"""Local FigureYa agent tools: profile, retrieve, inspect, fetch selected files, audit."""
import argparse
import collections
import csv
import datetime as dt
import hashlib
import json
import re
import sys
import urllib.parse
from pathlib import Path, PurePosixPath

from build_catalog import ROOT, blob_sha, read_url
from data_profile import profile


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def emit(value, output=None):
    data = json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)
    if output:
        target = Path(output)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(data + '\n', encoding='utf-8')
    else:
        print(data)


def catalog():
    return load(ROOT / 'catalog/modules.json'), load(ROOT / 'catalog/routing.json')['families']


def module_families(module, families, card):
    if 'families' in card:
        return [f for f in families if f['id'] in card['families']]
    return [f for f in families if module['number'] in f['numbers']]


def curated_cards(cat):
    data = load(ROOT / 'catalog/curated.json')
    current = data.get('reviewed_commit') == cat['commit']
    return data['modules'], current


def terms(query):
    # Preserve CJK phrases; English words require boundaries (e.g. GO != prognosis).
    return re.findall(r'[a-z0-9][a-z0-9_.+-]*|[\u3400-\u9fff]+', query.lower())


def matches(term, text):
    if re.search(r'[\u3400-\u9fff]', term):
        return term in text
    return bool(re.search(r'(?<![a-z0-9])' + re.escape(term) + r'(?![a-z0-9])', text))


def intent_families(query, families):
    return {f['id']: [k for k in f['keywords'] if matches(k.lower(), query.lower())] for f in families
            if any(matches(k.lower(), query.lower()) for k in f['keywords'])}


def contract_check(contract, p, explicit=None):
    if p is None:
        return {'fit': 0, 'status': 'unassessed', 'blockers': ['No data profile supplied'], 'evidence': []}
    roles = {k: list(v) for k, v in p.get('role_candidates', {}).items()}
    for k, v in (explicit or {}).items():
        roles[k] = v if isinstance(v, list) else [v]
    cols = {c['name']: c for c in p['columns']}
    for role, names in (explicit or {}).items():
        for name in names if isinstance(names, list) else [names]:
            if name not in cols:
                raise ValueError(f'Role {role}: unknown column {name}')
    numeric = [x for x in p['numeric_columns'] if x not in roles.get('id', []) and x not in roles.get('donor', [])]
    has = lambda k: bool(roles.get(k))
    evidence, blockers = [], []
    reqs = {
        'de_result': [('logfc',), ('padj', 'pvalue')],
        'survival': [('time',), ('event',)],
        'survival_score': [('time',), ('event',), ('score',)],
        'effect_ci': [('effect',), ('lower',), ('upper',)],
        'binary_score': [('label', 'event'), ('score',)],
        'enrichment': [('term',), ('padj', 'pvalue')],
        'edges': [('source',), ('target',)],
        'genomic': [('chromosome',), ('position',)],
    }
    fit = 0
    if contract in reqs:
        for choices in reqs[contract]:
            found = [k for k in choices if has(k)]
            if found:
                evidence.append('role candidate: ' + '/'.join(found))
            else:
                blockers.append('missing role: ' + '|'.join(choices))
        fit = round(25 * len(evidence) / len(reqs[contract]))
    elif contract == 'group_value':
        if numeric and (has('group') or p['category_candidates']):
            evidence.append('numeric measurement and candidate grouping exist'); fit = 25
        else:
            blockers.append('Need numeric measurements and grouping; summary-only tables are insufficient')
    elif contract == 'two_numeric':
        if len(numeric) >= 2:
            evidence.append('at least two numeric columns'); fit = 25
        else:
            blockers.append('Need two numeric measurements')
    elif contract == 'matrix':
        if len(numeric) >= 3:
            evidence.append('wide numeric matrix candidate; orientation unconfirmed'); fit = 25
        else:
            blockers.append('Need explicit matrix/coordinates or reshape plan; numeric-looking IDs excluded')
    elif contract == 'flow':
        if has('source') and has('target') or len(p['category_candidates']) >= 2:
            evidence.append('edge roles or multiple categorical columns'); fit = 25
        else:
            blockers.append('Need source/target or multiple categorical stages')
    elif contract == 'composition':
        if numeric and p['category_candidates']:
            evidence.append('category and numeric quantity candidates; denominator unconfirmed'); fit = 20
        else:
            blockers.append('Need categories and quantities/denominators')
    elif contract == 'sets':
        if has('feature') or any(c.get('binary_01') for c in cols.values()):
            evidence.append('feature IDs or binary membership candidate'); fit = 15
        else:
            blockers.append('Need explicit set membership; IDs alone do not establish sets')
    else:
        blockers.append('Specialized input object/analysis contract must be checked in source')
    if contract in {'de_result', 'enrichment'}:
        for role in ('padj', 'pvalue'):
            for name in roles.get(role, []):
                c = cols[name]
                if c['kind'] != 'numeric' or c.get('min', -1) < 0 or c.get('max', 2) > 1 or c['nonfinite']:
                    blockers.append(f'{name}: invalid P-value range/type')
                    fit = 0
    numeric_roles = {'de_result': ['logfc'], 'effect_ci': ['effect', 'lower', 'upper'],
                     'binary_score': ['score'], 'survival_score': ['score']}.get(contract, [])
    for role in numeric_roles:
        for name in roles.get(role, []):
            if cols[name]['kind'] != 'numeric' or cols[name]['nonfinite']:
                blockers.append(f'{name}: {role} must contain finite numeric values'); fit = 0
    if contract in {'survival', 'survival_score'}:
        for name in roles.get('time', []):
            c = cols[name]
            if c['kind'] != 'numeric' or c.get('min', -1) < 0 or c['nonfinite']:
                blockers.append(f'{name}: invalid follow-up time'); fit = 0
        for name in roles.get('event', []):
            if not cols[name].get('binary_01'):
                blockers.append(f'{name}: event coding must be mapped explicitly (including competing events)')
    if contract == 'binary_score':
        for name in roles.get('label', roles.get('event', [])):
            if cols[name]['distinct_observed'] != 2:
                blockers.append(f'{name}: outcome is not observed as binary'); fit = 0
    if contract in {'edges', 'flow', 'composition'}:
        for name in roles.get('weight', []):
            if cols[name]['kind'] != 'numeric' or cols[name].get('min', -1) < 0:
                blockers.append(f'{name}: invalid nonnegative weight'); fit = 0
    return {'fit': fit, 'status': 'blocked_or_unknown' if blockers else 'structural_candidate_only',
            'blockers': blockers, 'evidence': evidence}


def recommend(query, p=None, explicit=None, limit=8, design='unknown', domain=None):
    cat, families = catalog()
    cards, curation_current = curated_cards(cat)
    intents = intent_families(query, families)
    ranked = []
    for m in cat['modules']:
        card = cards.get(m['id'], {})
        mf = module_families(m, families, card)
        active = [f for f in mf if f['id'] in intents]
        text = ' '.join([m['id']] + [s.get('scenario_excerpt', '') + ' ' + s.get('input_excerpt', '') + ' ' + ' '.join(s.get('packages', [])) for s in m['scripts']]).lower()
        hits = [t for t in terms(query) if matches(t, text) or (len(t) >= 4 and t in m['id'].lower())]
        exact = m['id'].lower() in query.lower()
        number_ref = m['number'] in {int(n) for n in re.findall(r'figureya(\d+)', query.lower())}
        if not active and not hits and not exact and not number_ref:
            continue
        checks = [(f, contract_check(f['contract'], p, explicit)) for f in active or mf]
        best, check = max(checks, key=lambda item: item[1]['fit'], default=(None, {'fit': 0, 'status': 'unassessed', 'blockers': ['Unclassified; inspect script'], 'evidence': []}))
        # Heuristic retrieval only. No popularity, download status, or ease-of-coding bonus.
        purpose = 40 if active else 0
        lexical = min(15, 3 * len(hits))
        verified = 5 if m['evidence'] == 'source_metadata_extracted' else 0
        domain_fit = 0
        warnings = list(check['blockers'])
        if not curation_current:
            warnings.append('Curated cards refer to an older/unknown commit; re-read source before applying their distinctions')
        single_cell = any(f['id'] == 'single_cell' for f in mf) and not card.get('general_tabular')
        if domain in {'single-cell', 'spatial'}:
            domain_fit = 10 if single_cell else 0
        elif single_cell and domain not in {'single-cell', 'spatial'} and 'single_cell' not in intents:
            domain_fit = -15
            warnings.append('This module has single-cell prerequisites; domain not established')
        if design in {'paired', 'repeated'}:
            warnings.append('Preserve subject matching; adapt inference to paired/repeated observations')
        if best:
            warnings.append(best['cautions'])
        specific_hits = [k for k in card.get('prefer', []) if matches(k.lower(), query.lower())]
        specific = min(20, len(specific_hits) * 10) if curation_current else 0
        if card.get('note'):
            warnings.append(card['note'])
        id_score = 80 if exact else (40 if number_ref else 0)
        score = purpose + lexical + check['fit'] + verified + domain_fit + specific + id_score
        ranked.append({'id': m['id'], 'retrieval_score': score,
                       'score_components': {'purpose': purpose, 'text_match': lexical, 'data_structure': check['fit'], 'source_metadata': verified, 'domain': domain_fit, 'specific_semantics': specific, 'id_reference': id_score},
                       'families': [f['id'] for f in mf], 'matched_intent': [f['id'] for f in active],
                       'keyword_hits': hits, 'specific_hits': specific_hits, 'data_check': check,
                       'input_contract': best['needs'] if best else 'Inspect source',
                       'warnings': warnings, 'source_evidence': m['evidence'], 'source_url': m['source_url'],
                       'scripts': [s['path'] for s in m['scripts']],
                       'summary': ' '.join(s.get('scenario_excerpt', '') for s in m['scripts'])[:420]})
    ranked.sort(key=lambda x: (-x['retrieval_score'], x['id']))
    return {'query': query, 'commit': cat['commit'], 'curation_current': curation_current, 'matched_families': intents,
            'profile_used': p is not None, 'design': design, 'domain': domain,
            'decision': 'CANDIDATES_ONLY: agent must compare semantics, full source and visual reference; scores are not confidence or a final selection.',
            'needs': ['Confirm observational unit, roles, units and missingness', 'Compare at least 3 genuinely relevant modules when available; document rejection reasons', 'Inspect selected source and helper dependencies before fetching or running'],
            'candidates': ranked[:limit], 'matching_count': len(ranked)}


def inspect(module):
    cat, families = catalog()
    m = next((m for m in cat['modules'] if m['id'] == module), None)
    if m is None:
        raise ValueError('Unknown full module ID: ' + module)
    tree = load(ROOT / 'catalog/upstream-tree.json')
    cards, curation_current = curated_cards(cat)
    card = cards.get(module, {})
    return {'commit': cat['commit'], **m,
            'routing_families': module_families(m, families, card), 'curated_card': card, 'curation_current': curation_current,
            'files': [f for f in tree['files'] if f['path'].startswith(module + '/')],
            'note': 'excerpts/packages are mechanically extracted; absence in source metadata is not proof of absence. Preview/report URLs may change; pinned source is authoritative.'}


def safe_destination(base, relative):
    rel = PurePosixPath(relative)
    if rel.is_absolute() or '..' in rel.parts or '\\' in relative or ':' in relative:
        raise ValueError('Unsafe relative path')
    base = base.resolve()
    result = (base / relative).resolve()
    if not result.is_relative_to(base) or result == base:
        raise ValueError('Path escapes cache')
    return result


def fetch(module, paths, budget_mb=5):
    """Explicit allowlist; no archives, data, installers, recursive directory fetch or execution."""
    info = inspect(module)
    if not 0 < budget_mb <= 10:
        raise ValueError('Per-call budget must be >0 and <=10 MiB')
    if not paths or len(paths) > 6:
        raise ValueError('Select 1 to 6 exact code/document paths; inspect first')
    files = {f['path']: f for f in info['files']}
    chosen = []
    for path in dict.fromkeys(paths):
        if path not in files:
            raise ValueError('Path not in selected module inventory: ' + path)
        if PurePosixPath(path).suffix.lower() not in {'.rmd', '.r', '.py', '.md'}:
            raise ValueError('Only selected code/docs are supported; no datasets/reports/archives')
        chosen.append(files[path])
    if sum(f['size'] for f in chosen) > budget_mb * 1024**2:
        raise ValueError('Selected files exceed per-call byte budget')
    cache = ROOT / 'cache'
    cache.mkdir(exist_ok=True)
    cached_size = sum(f.stat().st_size for f in cache.rglob('*') if f.is_file())
    if cached_size + sum(f['size'] for f in chosen) > 32 * 1024**2:
        raise ValueError('32 MiB cache limit reached; review/delete obsolete selected files explicitly')
    results = []
    for f in chosen:
        target = safe_destination(cache, info['commit'] + '/' + f['path'])
        url = 'https://raw.githubusercontent.com/ying-ge/FigureYa/' + info['commit'] + '/' + urllib.parse.quote(f['path'])
        data = target.read_bytes() if target.exists() else read_url(url, min(512_000, int(budget_mb * 1024**2)))
        if len(data) != f['size'] or blob_sha(data) != f['sha']:
            raise ValueError('Source hash/size mismatch (or locally edited cache): ' + f['path'])
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            target.write_bytes(data)
        results.append({'path': str(target), 'upstream_path': f['path'], 'url': url, 'blob_sha': f['sha'],
                        'sha256': hashlib.sha256(data).hexdigest(), 'bytes': len(data)})
    manifest = {'module': module, 'commit': info['commit'], 'fetched_at': dt.datetime.now(dt.timezone.utc).isoformat(),
                'license': 'CC-BY-NC-SA-4.0 (upstream README; verify license for specific files)', 'executed': False, 'files': results}
    emit(manifest, cache / info['commit'] / module / ('fetch-' + dt.datetime.now().strftime('%Y%m%dT%H%M%S%f') + '.json'))
    return manifest


def audit(path):
    patterns = {
        'automatic_install': r'install[._](packages|github|version)|BiocManager::install|source\s*\([^\n]*install_dependencies',
        'working_directory': r'\bsetwd\s*\(',
        'absolute_path': r'[A-Z]:[/\\]|/Users/|/home/',
        'network_or_remote_source': r'https?://|download.file|GDCdownload',
        'external_command': r'\b(system2?|shell)\s*\(',
        'deletion': r'\b(unlink|file.remove)\s*\(',
        'workspace_clear': r'rm\s*\(\s*list\s*=',
        'helper_source': r'\bsource\s*\(',
        'demo_or_random_data': r'\b(rnorm|runif|rbinom|sample)\s*\(',
        'positional_or_fixed_slice': r'\[\s*\d+\s*:\s*\d+',
    }
    findings = []
    for i, line in enumerate(Path(path).read_text(encoding='utf-8-sig').splitlines(), 1):
        for label, pattern in patterns.items():
            if re.search(pattern, line, re.I):
                findings.append({'line': i, 'kind': label, 'excerpt': line.strip()[:240]})
    return {'path': str(Path(path).resolve()), 'findings': findings,
            'status': 'Static hints only; no execution. Review control flow, helpers, data semantics and formulas manually.'}


def export_index():
    cat, families = catalog()
    cards = load(ROOT / 'catalog/curated.json')['modules']
    lines = ['# FigureYa 模板索引', '', f"快照：`{cat['commit']}`；模块 {len(cat['modules'])} 个。", '',
             '本表为检索入口。分类按模块编号人工整理，说明/依赖由远程源码自动提取；不表示已经运行验证。',
             '先按数据契约筛选，再用 `inspect` 阅读完整条目。Plus/更新版必须使用完整 ID。', '',
             '| 模块 | 检索分组 | 依赖示例 | 源码证据 |', '|---|---|---|---|']
    for m in cat['modules']:
        cats = '、'.join(f['label'] for f in module_families(m, families, cards.get(m['id'], {}))) or '待语义复核；可全文检索'
        packages = ', '.join(sorted({p for s in m['scripts'] for p in s.get('packages', [])})[:5])
        lines.append(f"| [{m['id']}]({m['source_url']}) | {cats} | {packages} | {m['evidence']} |")
    (ROOT / 'catalog/INDEX.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return {'path': str(ROOT / 'catalog/INDEX.md'), 'modules': len(cat['modules'])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('profile'); p.add_argument('file'); p.add_argument('--out'); p.add_argument('--limit', type=int, default=10000)
    p.add_argument('--delimiter'); p.add_argument('--encoding', default='utf-8-sig'); p.add_argument('--sheet')
    p = sub.add_parser('recommend'); p.add_argument('--intent', required=True); p.add_argument('--profile'); p.add_argument('--roles')
    p.add_argument('--limit', type=int, default=8); p.add_argument('--design', choices=['unknown', 'independent', 'paired', 'repeated'], default='unknown')
    p.add_argument('--domain', choices=['general', 'bulk', 'single-cell', 'spatial', 'clinical']); p.add_argument('--out')
    p = sub.add_parser('inspect'); p.add_argument('module'); p.add_argument('--out')
    p = sub.add_parser('fetch'); p.add_argument('module'); p.add_argument('--path', action='append', required=True); p.add_argument('--budget-mb', type=float, default=5)
    p = sub.add_parser('audit'); p.add_argument('file'); p.add_argument('--out')
    sub.add_parser('export-index')
    args = parser.parse_args()
    if args.command == 'profile':
        emit(profile(args.file, args.limit, args.delimiter, args.encoding, args.sheet), args.out)
    elif args.command == 'recommend':
        if args.roles and not args.profile:
            raise ValueError('--roles requires --profile')
        emit(recommend(args.intent, load(args.profile) if args.profile else None, load(args.roles) if args.roles else None,
                       args.limit, args.design, args.domain), args.out)
    elif args.command == 'inspect':
        emit(inspect(args.module), args.out)
    elif args.command == 'fetch':
        emit(fetch(args.module, args.path, args.budget_mb))
    elif args.command == 'audit':
        emit(audit(args.file), args.out)
    else:
        emit(export_index())


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError, csv.Error) as e:
        print('ERROR: ' + str(e), file=sys.stderr)
        sys.exit(2)
