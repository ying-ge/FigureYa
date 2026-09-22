"""Behavioral regression tests using synthetic data, never the local heatmap example."""
import csv
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import figureya as agent
from build_catalog import extract, blob_sha
from data_profile import profile


class AgentTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def table(self, headers, rows, name='test.csv', **kwargs):
        path = self.root / name
        with path.open('w', encoding='utf-8-sig', newline='') as f:
            writer = csv.writer(f); writer.writerow(headers); writer.writerows(rows)
        return profile(path, **kwargs)

    def test_catalog_is_complete_and_paths_real(self):
        cat, families = agent.catalog()
        tree = agent.load(agent.ROOT / 'catalog/upstream-tree.json')
        self.assertFalse(tree['truncated'])
        self.assertEqual(cat['commit'], tree['commit'])
        paths = {x['path']: x for x in tree['files']}
        actual = {p.split('/')[0] for p in paths if p.startswith('FigureYa') and '/' in p}
        self.assertEqual(actual, {m['id'] for m in cat['modules']})
        for module in cat['modules']:
            for script in module['scripts']:
                self.assertIn(script['path'], paths)
                self.assertEqual(script['evidence'], 'source_metadata_extracted')
                self.assertEqual(script['blob_sha'], paths[script['path']]['sha'])
        ids = {m['id'] for m in cat['modules']}
        self.assertTrue(set(agent.load(agent.ROOT / 'catalog/curated.json')['modules']) <= ids)

    def test_same_shape_different_goal_changes_candidates(self):
        p = self.table(['method_a', 'method_b'], [[i, i + .2] for i in range(1, 41)])
        agreement = agent.recommend('测量方法一致性 Bland Altman', p)
        correlation = agent.recommend('两个变量相关 correlation', p)
        self.assertEqual(agreement['candidates'][0]['id'], 'FigureYa176BlandAltman')
        self.assertNotEqual(correlation['candidates'][0]['id'], agreement['candidates'][0]['id'])

    def test_raw_expression_cannot_be_treated_as_de_results(self):
        p = self.table(['gene', 's1', 's2', 's3'], [['A', 3, 4, 5], ['B', 5, 3, 2]])
        check = agent.contract_check('de_result', p)
        self.assertEqual(check['status'], 'blocked_or_unknown')
        self.assertTrue(any('logfc' in b for b in check['blockers']))

    def test_invalid_pvalues_block_volcano(self):
        p = self.table(['gene', 'log2FoldChange', 'padj'], [['A', 2.5, 1.2], ['B', -1, .03]])
        check = agent.contract_check('de_result', p)
        self.assertEqual(check['fit'], 0)
        self.assertTrue(any('invalid P-value' in b for b in check['blockers']))

    def test_good_de_and_multi_information_route(self):
        p = self.table(['gene', 'logFC', 'padj', 'feature'], [[str(i), i - 4, .01, 'coding'] for i in range(10)])
        check = agent.contract_check('de_result', p)
        self.assertEqual(check['status'], 'structural_candidate_only')
        classic = agent.recommend('ggplot2 火山图 基因标签 候选基因', p)
        multi = agent.recommend('火山图 边际 百分比 feature', p)
        self.assertEqual(classic['candidates'][0]['id'], 'FigureYa59volcanoV2')
        self.assertEqual(multi['candidates'][0]['id'], 'FigureYa135multiVolcano')
        self.assertNotIn('FigureYa59Plus_GEO2DEG', [c['id'] for c in classic['candidates'][:3]])

    def test_survival_not_inferred_from_time_alone(self):
        p = self.table(['time', 'value'], [[1, 2], [2, 4], [3, 5]])
        check = agent.contract_check('survival', p)
        self.assertTrue(any('event' in x for x in check['blockers']))

    def test_survival_and_binary_roc_differ(self):
        p = self.table(['time', 'event', 'score'], [[10, 1, .8], [20, 0, .2], [15, 1, .7], [40, 0, .1]])
        a = agent.recommend('删失结局 时间依赖 ROC', p)
        self.assertEqual(a['candidates'][0]['id'], 'FigureYa85timeROC')
        self.assertEqual(agent.contract_check('survival_score', p)['status'], 'structural_candidate_only')

    def test_pairing_is_carried_forward(self):
        p = self.table(['subject_id', 'group', 'value'], [[str(i//2), 'pre' if i%2==0 else 'post', i+.2] for i in range(12)])
        result = agent.recommend('配对测量 分布 连线', p, design='paired')
        self.assertTrue(result['candidates'])
        self.assertTrue(all(any('subject matching' in s for s in c['warnings']) for c in result['candidates']))

    def test_heatmap_variants_are_not_all_basic(self):
        p = self.table(['gene', 'a', 'b', 'c'], [[str(i), i, i+1, i+2] for i in range(12)])
        cases = [('热图 双矩阵 上下三角', 'FigureYa144DiagHeatmap'),
                 ('热图 叠加气泡 两套数据 双指标', 'FigureYa278heatmapPoints'),
                 ('热图 gistic cnv 拷贝数', 'FigureYa307CNVHeatmap')]
        for query, expected in cases:
            self.assertEqual(agent.recommend(query, p)['candidates'][0]['id'], expected)

    def test_batch_pca_and_ecdf_have_specific_candidates(self):
        p = self.table(['id', 'a', 'b', 'c', 'group'], [[str(i), i, i+1, i+2, str(i%2)] for i in range(20)])
        self.assertEqual(agent.recommend('PCA 分组 批次 batch 形状', p)['candidates'][0]['id'], 'FigureYa101PCA')
        self.assertEqual(agent.recommend('分布 ecdf 累积分布 不平滑', p)['candidates'][0]['id'], 'FigureYa298ecdfPvalue')

    def test_profile_bounds_duplicates_missing_and_infinite(self):
        p = self.table(['id', 'value'], [['a', ''], ['b', 'inf'], ['b', 'inf'], ['tail', '999']], limit=3)
        self.assertTrue(p['sampled']); self.assertIsNone(p['row_count_exact'])
        self.assertEqual(p['duplicate_rows_observed'], 1)
        self.assertEqual(p['columns'][1]['missing'], 1)
        self.assertEqual(p['columns'][1]['nonfinite'], 2)
        json.dumps(p, allow_nan=False)

    def test_header_and_row_shape_errors_are_not_silent(self):
        with self.assertRaises(ValueError):
            self.table(['a', 'a'], [[1, 2]])
        with self.assertRaises(ValueError):
            self.table(['a', 'b'], [[1, 2, 3]])

    def test_roles_are_checked_and_numeric_ids_excluded(self):
        p = self.table(['patient_id', 'm'], [[1, 8], [2, 9], [3, 10]])
        self.assertEqual(agent.contract_check('two_numeric', p)['fit'], 0)
        with self.assertRaises(ValueError):
            agent.contract_check('survival', p, {'time': 'not_a_column'})

    def test_negative_flow_weights_rejected(self):
        p = self.table(['source', 'target', 'weight'], [['A', 'B', -5], ['A', 'C', 3]])
        self.assertEqual(agent.contract_check('flow', p)['fit'], 0)

    def test_unknown_query_does_not_default_to_heatmap(self):
        self.assertEqual(agent.recommend('zzzz_nonexistent_plot_family_zzzz')['candidates'], [])

    def test_module_number_and_specialized_name_are_searchable(self):
        self.assertEqual(agent.recommend('FigureYa186')['candidates'][0]['id'], 'FigureYa186swimmerplot')
        self.assertIn('FigureYa186swimmerplot', [c['id'] for c in agent.recommend('swimmer')['candidates']])

    def test_nonnumeric_effect_is_blocked(self):
        p = self.table(['gene', 'logFC', 'padj'], [['A', 'high', .03], ['B', 'low', .02]])
        self.assertEqual(agent.contract_check('de_result', p)['fit'], 0)

    def test_audit_finds_install_and_fixed_paths(self):
        path = self.root/'synthetic.R'
        path.write_text('source("install_dependencies.R")\nsetwd("D:/author/project")\nx <- rnorm(20)\n', encoding='utf-8')
        kinds = {f['kind'] for f in agent.audit(path)['findings']}
        self.assertTrue({'automatic_install', 'working_directory', 'absolute_path', 'demo_or_random_data'} <= kinds)

    def test_source_extraction_does_not_keep_executable_body(self):
        source = '---\ntitle: "Example"\n---\n## Input data\nA matrix\n```{r}\nlibrary(pheatmap)\nx <- read.csv("data.csv")\nsecret_execute_this()\n```\n'
        result = extract(source)
        self.assertIn('pheatmap', result['packages'])
        self.assertEqual(result['input_literals'], ['data.csv'])
        self.assertNotIn('secret_execute_this', json.dumps(result))

    def test_path_traversal_rejected(self):
        for path in ['../escape.R', '/tmp/a.R', 'C:/tmp/a.R', 'x\\..\\a.R']:
            with self.assertRaises(ValueError):
                agent.safe_destination(self.root, path)

    def test_fetch_is_exact_hashed_and_never_executes(self):
        code = b'print("synthetic source; do not execute")\n'
        info = {'commit': 'a'*40, 'files': [{'path': 'FigureYaTest/test.R', 'size': len(code), 'sha': blob_sha(code)},
                                         {'path': 'FigureYaTest/data.csv', 'size': 100, 'sha': 'b'*40}]}
        with patch.object(agent, 'ROOT', self.root), patch.object(agent, 'inspect', return_value=info), patch.object(agent, 'read_url', return_value=code) as network:
            result = agent.fetch('FigureYaTest', ['FigureYaTest/test.R'])
            self.assertFalse(result['executed'])
            self.assertEqual(result['files'][0]['sha256'], hashlib.sha256(code).hexdigest())
            agent.fetch('FigureYaTest', ['FigureYaTest/test.R'])
            self.assertEqual(network.call_count, 1)
            with self.assertRaises(ValueError):
                agent.fetch('FigureYaTest', ['FigureYaTest/data.csv'])
            with self.assertRaises(ValueError):
                agent.fetch('FigureYaTest', ['FigureYaTest/not_found.R'])

    def test_fetch_hash_mismatch_fails_without_writing_source(self):
        info = {'commit': 'a'*40, 'files': [{'path': 'FigureYaTest/test.R', 'size': 4, 'sha': '0'*40}]}
        with patch.object(agent, 'ROOT', self.root), patch.object(agent, 'inspect', return_value=info), patch.object(agent, 'read_url', return_value=b'bad!'):
            with self.assertRaises(ValueError):
                agent.fetch('FigureYaTest', ['FigureYaTest/test.R'])
            self.assertFalse(list(self.root.rglob('*.R')))

    def test_xlsx_multisheet_requires_explicit_selection(self):
        try:
            import openpyxl
        except ImportError:
            self.skipTest('optional openpyxl absent')
        book = openpyxl.Workbook()
        book.active.title = 'one'; book.active.append(['x']); book.active.append([1]); book.create_sheet('two')
        path = self.root/'multi.xlsx'; book.save(path); book.close()
        with self.assertRaises(ValueError):
            profile(path)
        self.assertEqual(profile(path, sheet='one')['rows_observed'], 1)


if __name__ == '__main__':
    unittest.main()
