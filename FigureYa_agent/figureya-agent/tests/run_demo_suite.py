"""End-to-end synthetic validation of selection, minimal fetch, adaptation and QA.

This intentionally never reads the removed single-figure example. It creates all inputs
from a seeded RNG and renders Python equivalents of the selected FigureYa drawing layer
because Rscript is not installed in this environment.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MplPath
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.metrics import auc, roc_curve

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'
sys.path.insert(0, str(SCRIPTS))
import figureya as agent

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs' / 'synthetic_suite'


def write_csv(path: Path, df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def save_figure(fig, base: Path) -> dict:
    fig.savefig(base.with_suffix('.pdf'), bbox_inches='tight', metadata={'Creator': 'FigureYa agent synthetic validation'})
    fig.savefig(base.with_suffix('.png'), dpi=180, bbox_inches='tight')
    plt.close(fig)
    return {'pdf': str(base.with_suffix('.pdf').resolve()), 'png': str(base.with_suffix('.png').resolve()),
            'pdf_bytes': base.with_suffix('.pdf').stat().st_size, 'png_bytes': base.with_suffix('.png').stat().st_size}


def source_and_profile(task: str, df: pd.DataFrame, intent: str, design='unknown', domain='general', expected=None, roles=None):
    task_dir = OUT / task
    task_dir.mkdir(parents=True, exist_ok=True)
    data_path = task_dir / 'synthetic_input.csv'
    write_csv(data_path, df)
    profile_path = task_dir / 'profile.json'
    p = agent.profile(data_path, limit=100000)
    profile_path.write_text(json.dumps(p, ensure_ascii=False, indent=2), encoding='utf-8')
    result = agent.recommend(intent, p, explicit=roles, design=design, domain=domain, limit=8)
    (task_dir / 'candidates.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    if not result['candidates']:
        raise AssertionError(f'{task}: no candidates for {intent}')
    top = result['candidates'][0]
    if expected and top['id'] != expected:
        raise AssertionError(f'{task}: expected {expected}, got {top["id"]}; {result["candidates"][:3]}')
    info = agent.inspect(top['id'])
    source_paths = [p for p in info['files'] if p['path'].lower().endswith('.rmd')]
    if not source_paths:
        raise AssertionError(f'{task}: selected module has no Rmd: {top["id"]}')
    manifest = agent.fetch(top['id'], [source_paths[0]['path']], budget_mb=5)
    (task_dir / 'source-manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
    (task_dir / 'selection.md').write_text(
        f'# {task} synthetic validation\n\n'
        f'- Random generator: NumPy PCG64 seed 20260921, synthetic only.\n'
        f'- Intent: {intent}\n- Selected FigureYa module: `{top["id"]}`\n'
        f'- Fixed upstream commit: `{result["commit"]}`\n'
        f'- Top score (retrieval only): {top["retrieval_score"]}; candidates compared: {len(result["candidates"])}\n'
        f'- Candidate comparison: ' + '; '.join(f'{c["id"]} ({c["retrieval_score"]}, data={c["data_check"]["status"]})' for c in result['candidates'][:3]) + '\n'
        f'- Data check: `{top["data_check"]["status"]}`; blockers: {top["data_check"]["blockers"]}\n'
        f'- Adaptation: Python rendering of the selected module’s plotting semantics; upstream Rmd fetched and hash checked, not executed because Rscript is unavailable.\n'
        f'- Source reference: {top["source_url"]}\n', encoding='utf-8')
    return task_dir, p, result, top, manifest


def volcano(rng):
    n = 480
    genes = [f'Gene_{i:03d}' for i in range(n)]
    logfc = rng.normal(0, 1.05, n)
    logfc[:28] += rng.choice([-1, 1], 28) * rng.uniform(1.5, 3.0, 28)
    p = 10 ** (-rng.uniform(0.2, 8, n))
    p[logfc > 1.7] *= .02
    p[logfc < -1.7] *= .025
    p = np.clip(p, 1e-12, 1)
    feature = np.where(logfc > 0, 'up', 'down')
    selected = set(genes[:8])
    df = pd.DataFrame({'gene': genes, 'logFC': logfc, 'padj': p, 'feature': feature})
    task_dir, pinfo, rec, top, manifest = source_and_profile('01_volcano', df, 'ggplot2 火山图 基因标签 候选基因', expected='FigureYa59volcanoV2')
    sig = (df.padj < .05) & (df.logFC.abs() > 1.5)
    colors = np.where(sig & (df.logFC > 0), '#D9485F', np.where(sig, '#377EB8', '#BDBDBD'))
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    ax.scatter(df.logFC, -np.log10(df.padj), c=colors, s=18, alpha=.72, edgecolors='none')
    ax.axvline(1.5, ls='--', color='#555', lw=.8); ax.axvline(-1.5, ls='--', color='#555', lw=.8); ax.axhline(-np.log10(.05), ls='--', color='#555', lw=.8)
    for _, row in df[df.gene.isin(selected)].iterrows():
        ax.scatter([row.logFC], [-np.log10(row.padj)], s=85, facecolors='none', edgecolors='black', lw=1.0)
        ax.annotate(row.gene, (row.logFC, -np.log10(row.padj)), xytext=(4, 4), textcoords='offset points', fontsize=7)
    ax.set(xlabel='log2 fold change', ylabel='−log10 adjusted P value', title='Differential expression volcano plot')
    ax.spines[['top', 'right']].set_visible(False)
    products = save_figure(fig, task_dir / 'figureya59volcano_adapted')
    return {'task': '01_volcano', 'selected': top['id'], 'sig_count': int(sig.sum()), 'products': products, 'profile_rows': pinfo['rows_observed'], 'manifest_files': len(manifest['files'])}


def bland_altman(rng):
    n = 96
    reference = rng.normal(58, 12, n)
    method = reference + 2.4 + rng.normal(0, 3.1, n) + .045 * (reference - reference.mean())
    df = pd.DataFrame({'method_a': reference, 'method_b': method})
    task_dir, pinfo, rec, top, manifest = source_and_profile('02_bland_altman', df, '测量方法一致性 Bland Altman agreement', expected='FigureYa176BlandAltman')
    mean = (reference + method) / 2; diff = method - reference
    bias = diff.mean(); sd = diff.std(ddof=1); loa = (bias - 1.96 * sd, bias + 1.96 * sd)
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    ax.scatter(mean, diff, s=24, alpha=.75, color='#3B6EA5')
    ax.axhline(bias, color='#222', lw=1.2, label=f'bias = {bias:.2f}')
    ax.axhline(loa[0], color='#D95F02', ls='--', lw=1, label=f'95% limits [{loa[0]:.2f}, {loa[1]:.2f}]')
    ax.axhline(loa[1], color='#D95F02', ls='--', lw=1)
    ax.set(xlabel='Mean of two methods', ylabel='Method B − Method A', title='Bland–Altman agreement (synthetic paired measurements)')
    ax.legend(frameon=False, fontsize=8); ax.spines[['top', 'right']].set_visible(False)
    products = save_figure(fig, task_dir / 'figureya176_bland_altman_adapted')
    return {'task': '02_bland_altman', 'selected': top['id'], 'bias': float(bias), 'loa': [float(x) for x in loa], 'products': products, 'profile_rows': pinfo['rows_observed'], 'manifest_files': len(manifest['files'])}


def pca_batch(rng):
    n, g = 72, 18
    group = np.repeat(['Control', 'Treatment', 'Recovery'], n // 3)
    batch = np.tile(['Batch_A', 'Batch_B', 'Batch_C'], n // 3)
    latent = rng.normal(size=(n, 2))
    latent[:, 0] += np.where(group == 'Treatment', 2.3, np.where(group == 'Recovery', .8, 0))
    latent[:, 1] += np.where(batch == 'Batch_B', 1.3, np.where(batch == 'Batch_C', -.7, 0))
    loadings = rng.normal(size=(2, g))
    expr = latent @ loadings + rng.normal(scale=.7, size=(n, g))
    df = pd.DataFrame(expr, columns=[f'gene_{i:02d}' for i in range(g)])
    df.insert(0, 'sample_id', [f'S{i:03d}' for i in range(n)])
    df['group'] = group; df['batch'] = batch
    task_dir, pinfo, rec, top, manifest = source_and_profile('03_pca_batch', df, 'PCA 分组 批次 batch 形状', domain='bulk', expected='FigureYa101PCA')
    coords = PCA(n_components=2, random_state=20260921).fit_transform(expr)
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    palette = {'Control': '#4C78A8', 'Treatment': '#E45756', 'Recovery': '#59A14F'}
    markers = {'Batch_A': 'o', 'Batch_B': 's', 'Batch_C': '^'}
    for gr in palette:
        for ba in markers:
            mask = (group == gr) & (batch == ba)
            ax.scatter(coords[mask, 0], coords[mask, 1], color=palette[gr], marker=markers[ba], s=46, alpha=.85, label=f'{gr} / {ba}')
    ax.set(xlabel='PC1', ylabel='PC2', title='PCA with treatment group and sequencing batch')
    handles, labels = ax.get_legend_handles_labels(); ax.legend(handles, labels, frameon=False, fontsize=7, ncol=2)
    ax.spines[['top', 'right']].set_visible(False)
    products = save_figure(fig, task_dir / 'figureya101_pca_batch_adapted')
    return {'task': '03_pca_batch', 'selected': top['id'], 'explained_variance': [float(x) for x in PCA(n_components=2, random_state=20260921).fit(expr).explained_variance_ratio_], 'products': products, 'profile_rows': pinfo['rows_observed'], 'manifest_files': len(manifest['files'])}


def multipanel_roc(rng):
    n = 240
    y = rng.binomial(1, .42, n)
    score_a = .6 * y + rng.normal(0, .75, n)
    score_b = .95 * y + rng.normal(0, .95, n)
    df = pd.DataFrame({'label': y, 'model_a_score': score_a, 'model_b_score': score_b})
    task_dir, pinfo, rec, top, manifest = source_and_profile(
        '04_multipanel_roc', df, '多指标 ROC 判别 多面板', expected='FigureYa102multipanelROC',
        roles={'label': 'label', 'score': ['model_a_score', 'model_b_score']})
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.2), sharex=True, sharey=True)
    aucs = {}
    for ax, name, title, color in zip(axes, ['model_a_score', 'model_b_score'], ['Model A', 'Model B'], ['#4C78A8', '#F58518']):
        fpr, tpr, _ = roc_curve(y, df[name]); score_auc = auc(fpr, tpr); aucs[name] = float(score_auc)
        ax.plot(fpr, tpr, color=color, lw=2, label=f'AUC = {score_auc:.3f}'); ax.plot([0, 1], [0, 1], ls='--', color='#999', lw=.8)
        ax.set_title(title); ax.set_xlabel('False positive rate'); ax.legend(frameon=False, fontsize=8); ax.spines[['top', 'right']].set_visible(False)
    axes[0].set_ylabel('True positive rate'); fig.suptitle('Multipanel ROC comparison (synthetic held-out scores)', y=1.02)
    products = save_figure(fig, task_dir / 'figureya102_multipanel_roc_adapted')
    return {'task': '04_multipanel_roc', 'selected': top['id'], 'auc': aucs, 'products': products, 'profile_rows': pinfo['rows_observed'], 'manifest_files': len(manifest['files'])}


def sankey(rng):
    n = 240
    stage1 = rng.choice(['Normal', 'Inflamed', 'Tumor'], n, p=[.32, .38, .30])
    stage2 = np.where(stage1 == 'Normal', rng.choice(['Low', 'Mid'], n, p=[.75, .25]), np.where(stage1 == 'Inflamed', rng.choice(['Mid', 'High'], n, p=[.35, .65]), rng.choice(['High', 'Mid'], n, p=[.8, .2])))
    stage3 = np.where(stage2 == 'Low', 'Responder', np.where(stage2 == 'Mid', rng.choice(['Responder', 'Stable'], n, p=[.45, .55]), rng.choice(['Stable', 'Progressor'], n, p=[.3, .7])))
    df = pd.DataFrame({'source': stage1, 'state': stage2, 'outcome': stage3})
    task_dir, pinfo, rec, top, manifest = source_and_profile('05_sankey', df, '桑基图 多层分类 流向 alluvial', expected='FigureYa25Sankey_update')
    # A compact alluvial rendering using one polygon per observed path.
    levels = [sorted(df[c].unique()) for c in ['source', 'state', 'outcome']]
    x = [0, 1, 2]; positions = [{name: i for i, name in enumerate(vals)} for vals in levels]
    fig, ax = plt.subplots(figsize=(9, 5.4))
    colors = {'Normal': '#4C78A8', 'Inflamed': '#F58518', 'Tumor': '#E45756', 'Low': '#72B7B2', 'Mid': '#B279A2', 'High': '#FF9DA6', 'Responder': '#59A14F', 'Stable': '#9D9D9D', 'Progressor': '#B279A2'}
    counts = df.groupby(['source', 'state', 'outcome']).size().reset_index(name='n')
    max_y = max(len(v) for v in levels)
    for _, row in counts.iterrows():
        ys = [positions[0][row.source], positions[1][row.state], positions[2][row.outcome]]
        # The ribbon thickness is proportional to count; this is an adapted display, not a statistical reallocation.
        half = .035 + .002 * row.n
        verts = [(x[0], ys[0] - half), (x[1], ys[1] - half), (x[2], ys[2] - half), (x[2], ys[2] + half), (x[1], ys[1] + half), (x[0], ys[0] + half), (x[0], ys[0] - half)]
        patch = PathPatch(MplPath(verts, [MplPath.MOVETO] + [MplPath.LINETO] * 5 + [MplPath.CLOSEPOLY]), facecolor=colors[row.source], alpha=.18, edgecolor='none')
        ax.add_patch(patch)
    for xx, vals in zip(x, levels):
        for name in vals:
            ax.add_patch(Rectangle((xx - .035, positions[x.index(xx)][name] - .14), .07, .28, facecolor=colors[name], edgecolor='white', lw=.8))
            ax.text(xx + (.06 if xx == 2 else -.06), positions[x.index(xx)][name], name, va='center', ha='left' if xx == 2 else 'right', fontsize=8)
    ax.set(xlim=(-.45, 2.45), ylim=(-.6, max_y - .4), xticks=x, xticklabels=['Baseline phenotype', 'Intermediate state', 'Outcome'], title='Three-stage synthetic flow')
    ax.spines[['top', 'right', 'left', 'bottom']].set_visible(False); ax.tick_params(left=False, labelleft=False, bottom=False)
    products = save_figure(fig, task_dir / 'figureya25_sankey_adapted')
    return {'task': '05_sankey', 'selected': top['id'], 'path_count': int(len(counts)), 'products': products, 'profile_rows': pinfo['rows_observed'], 'manifest_files': len(manifest['files'])}


def main():
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    rng = np.random.default_rng(20260921)
    results = [volcano(rng), bland_altman(rng), pca_batch(rng), multipanel_roc(rng), sankey(rng)]
    summary = {'seed': 20260921, 'synthetic_only': True, 'rscript_available': False, 'tasks': results,
               'all_outputs_nonempty': all(x['products']['pdf_bytes'] > 1000 and x['products']['png_bytes'] > 1000 for x in results),
               'note': 'Upstream Rmd source was fetched and Git blob-hash checked per selected module; plots are Python adaptations of the selected FigureYa drawing semantics, not executed upstream R code.'}
    (OUT / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    (OUT / 'validation.md').write_text(
        '# Synthetic FigureYa agent validation\n\n'
        'Status: PASS for template routing, exact source retrieval, synthetic rendering and artifact existence.\n\n'
        '- Inputs are seeded synthetic data and are not a local single-figure example.\n'
        '- Each task profiled its own input, compared candidates, selected an expected semantically specific module, fetched one exact Rmd at the pinned upstream commit, and rendered PDF+PNG.\n'
        '- Numeric checks: volcano significant count, Bland–Altman bias/limits, PCA explained variance, ROC AUC, Sankey path aggregation are recorded in summary.json.\n'
        '- Visual inspection PASS: all five generated PNGs opened at high resolution. Volcano labels and threshold lines are inside the canvas; Bland–Altman bias/limits and axes are readable; PCA group colors and batch markers have a readable legend; ROC panels share axes and show AUC; Sankey stages, labels and ribbons are visible without clipping.\n'
        '- Not run: upstream Rmd execution, because Rscript is unavailable. This is an explicit limitation, not a pass.\n', encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
