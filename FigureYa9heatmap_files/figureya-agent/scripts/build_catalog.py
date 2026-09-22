"""Browse upstream Rmd in memory; retain searchable metadata, never source copies."""
import argparse
import concurrent.futures
import datetime as dt
import hashlib
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPO = 'ying-ge/FigureYa'


def read_url(url, limit=4_000_000):
    for attempt in range(3):
        try:
            request = urllib.request.Request(url, headers={'User-Agent': 'FigureYa-local-catalog/1.0'})
            with urllib.request.urlopen(request, timeout=35) as response:
                data = response.read(limit + 1)
            if len(data) > limit:
                raise ValueError('Remote file exceeds byte limit')
            return data
        except Exception:
            if attempt == 2:
                raise
            time.sleep(attempt + 1)


def blob_sha(data):
    return hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()


def extract(source):
    # Markdown prose is separated from executable chunks before extraction.
    prose = re.sub(r'```.*?```', '', source, flags=re.S)
    sections = re.split(r'^#{1,6}\s+', prose, flags=re.M)
    selected = []
    inputs = []
    for section in sections:
        heading, _, body = section.partition('\n')
        body = re.sub(r'!\[[^\]]*\]\([^)]*\)', '', body)
        body = re.sub(r'https?://\S+', '[link]', body)
        body = re.sub(r'\s+', ' ', body).strip()
        if re.search('需求|应用场景|使用场景|requirement|scenario', heading, re.I):
            selected.append(body[:650])
        if re.search('输入|input', heading, re.I):
            inputs.append(body[:700])
    title = re.search(r'^title:\s*(.+)$', source, re.M)
    packages = sorted(set(re.findall(r'(?:library|require)\(\s*[\'"]?([\w.]+)', source)))
    packages += sorted(set(re.findall(r'\b([A-Za-z][\w.]*)::', source)) - set(packages))
    read_calls = re.findall(r'(?:read[.\w]*|fread|load)\s*\(\s*(?:file\s*=\s*)?[\'"]([^\'"\n]+)', source)
    source_calls = re.findall(r'source\s*\(\s*[\'"]([^\'"\n]+)', source)
    functions = sorted(set(re.findall(r'\b((?:geom_|stat_)[A-Za-z_]+|pheatmap|Heatmap|ggsurvplot|ggforest|forestplot|corrplot|ggalluvial|oncoPrint|ggVennDiagram|ggvenn|upset|plotROC|ggroc|BlandAltman|survfit|coxph|prcomp|cor|cor.test)\s*\(', source)))
    patterns = {
        'auto_install': r'install[._](?:packages|github)|BiocManager::install|source\([\'"]install_dependencies',
        'working_directory': r'\bsetwd\s*\(',
        'network': r'download.file|https?://|GDCdownload',
        'shell': r'\b(?:system2?|shell)\s*\(',
        'delete': r'\b(?:unlink|file.remove)\s*\(',
        'workspace_clear': r'rm\s*\(\s*list\s*=',
        'demo_simulation': r'\b(?:rnorm|runif|rbinom|sample)\s*\(',
    }
    return {
        'title': title.group(1).strip('"\' ') if title else '',
        'scenario_excerpt': ' / '.join(selected)[:1800],
        'input_excerpt': ' / '.join(inputs)[:1800],
        'packages': packages, 'input_literals': sorted(set(read_calls)),
        'source_literals': sorted(set(source_calls)), 'plot_functions': functions,
        'review_flags': [k for k, v in patterns.items() if re.search(v, source)],
        'evidence': 'source_metadata_extracted',
    }


def build(refresh=False, workers=6):
    directory = ROOT / 'catalog'
    directory.mkdir(exist_ok=True)
    tree_path = directory / 'upstream-tree.json'
    if refresh or not tree_path.exists():
        api = 'https://api.github.com/repos/' + REPO
        repo = json.loads(read_url(api))
        commit = json.loads(read_url(api + '/commits/' + repo['default_branch']))
        tree = json.loads(read_url(api + '/git/trees/' + commit['sha'] + '?recursive=1'))
        if tree.get('truncated'):
            raise ValueError('Truncated tree; refuse to publish an incomplete catalog')
        meta = {'repository': REPO, 'commit': commit['sha'], 'commit_date': commit['commit']['committer']['date'],
                'branch': repo['default_branch'], 'truncated': False,
                'files': [{k: x[k] for k in ('path', 'sha', 'size')} for x in tree['tree'] if x['type'] == 'blob']}
    else:
        meta = json.loads(tree_path.read_text(encoding='utf-8'))
    if meta.get('truncated'):
        raise ValueError('Incomplete tree')
    modules = {}
    for f in meta['files']:
        folder = f['path'].split('/')[0]
        if '/' in f['path'] and re.match(r'^FigureYa\d', folder):
            modules.setdefault(folder, []).append(f)
    def scan(item):
        folder, files = item
        rmds = [f for f in files if f['path'].lower().endswith('.rmd')]
        # Index every Rmd, including alternative scripts, but do not retain its code.
        records = []
        for f in rmds:
            url = 'https://raw.githubusercontent.com/' + REPO + '/' + meta['commit'] + '/' + urllib.parse.quote(f['path'])
            try:
                data = read_url(url, 512_000)
                if blob_sha(data) != f['sha']:
                    raise ValueError('Blob hash mismatch')
                record = extract(data.decode('utf-8-sig', errors='replace'))
                record.update(path=f['path'], blob_sha=f['sha'], bytes=f['size'], url=url)
            except Exception as e:
                record = {'path': f['path'], 'evidence': 'unavailable', 'error': str(e), 'url': url}
            records.append(record)
        return {'id': folder, 'number': int(re.search(r'FigureYa(\d+)', folder).group(1)),
                'scripts': records, 'file_count': len(files), 'upstream_bytes': sum(f['size'] for f in files),
                'source_url': 'https://github.com/' + REPO + '/tree/' + meta['commit'] + '/' + urllib.parse.quote(folder),
                'report_paths': [f['path'] for f in files if f['path'].lower().endswith('.html')],
                'evidence': 'source_metadata_extracted' if records and all(r['evidence'] != 'unavailable' for r in records) else 'inventory_only'}
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for index, item in enumerate(pool.map(scan, modules.items()), 1):
            results.append(item)
            if index % 40 == 0:
                print(f'Browsed {index}/{len(modules)} modules (source retained: 0)', flush=True)
    results.sort(key=lambda r: (r['number'], r['id']))
    catalog = {'schema_version': 1, 'repository': REPO, 'commit': meta['commit'], 'commit_date': meta['commit_date'],
               'indexed_at': dt.datetime.now(dt.timezone.utc).isoformat(), 'source_code_retained': False,
               'modules': results}
    # Publish only after browsing completes. Failures remain explicitly marked.
    for path, value in [(tree_path, meta), (directory / 'modules.json', catalog)]:
        tmp = path.with_suffix('.tmp')
        tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding='utf-8')
        tmp.replace(path)
    print(json.dumps({'modules': len(results), 'scripts': sum(len(x['scripts']) for x in results),
                      'unavailable': [s['path'] for m in results for s in m['scripts'] if s['evidence'] == 'unavailable']}, ensure_ascii=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refresh', action='store_true', help='Refresh metadata and browse current upstream; no clone/download archive')
    parser.add_argument('--workers', type=int, choices=range(1, 9), default=6)
    args = parser.parse_args()
    build(args.refresh, args.workers)
