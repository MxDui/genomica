import math
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import run_analysis  # noqa: E402


class MetadataParsingTests(unittest.TestCase):
    def test_title_to_short_id_for_covid_and_control_samples(self):
        self.assertEqual(run_analysis.title_to_short_id("COVID_08_78y_male_ICU"), "C8")
        self.assertEqual(run_analysis.title_to_short_id("NONCOVID_12_70y_female_ICU"), "NC12")

    def test_title_to_short_id_rejects_unexpected_format(self):
        with self.assertRaises(ValueError):
            run_analysis.title_to_short_id("sample_without_expected_prefix")

    def test_normalize_key(self):
        self.assertEqual(run_analysis.normalize_key("Charlson score"), "charlson_score")
        self.assertEqual(run_analysis.normalize_key("Age (years)"), "age_years")

    def test_parse_age(self):
        self.assertEqual(run_analysis.parse_age("63y"), 63)
        self.assertEqual(run_analysis.parse_age(">89"), 89.0)
        self.assertTrue(math.isnan(run_analysis.parse_age(None)))


class MetricTests(unittest.TestCase):
    def test_cluster_accuracy_is_label_permutation_invariant(self):
        labels = np.array([0, 0, 1, 1])
        truth = np.array([1, 1, 0, 0])
        self.assertEqual(run_analysis.cluster_accuracy_for_binary(labels, truth), 1.0)

    def test_cluster_accuracy_requires_binary_labels(self):
        labels = np.array([0, 1, 2, 2])
        truth = np.array([0, 1, 1, 0])
        self.assertTrue(math.isnan(run_analysis.cluster_accuracy_for_binary(labels, truth)))


if __name__ == "__main__":
    unittest.main()

