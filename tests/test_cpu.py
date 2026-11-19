
import unittest

from src.cpu_engine import find_matches


class TestCPUEngine(unittest.TestCase):

    def setUp(self):
        self.genome = "ATCGATCGATCGGATCCGATCGATCGATCG"
        self.guide = "ATCGAT"

    def test_exact_matches(self):
        results, elapsed = find_matches(
            self.genome, self.guide, 0
        )
        positions = [r["position"] for r in results]

        self.assertEqual(positions, [0, 4, 18, 22])

    def test_mismatch_tolerance(self):
        results, elapsed = find_matches(
            self.genome, self.guide, 1
        )

        self.assertTrue(all(r["mismatches"] <= 1 for r in results))

    def test_empty_genome(self):
        with self.assertRaises(ValueError):
            find_matches("", self.guide)

    def test_empty_guide(self):
        with self.assertRaises(ValueError):
            find_matches(self.genome, "")

    def test_guide_longer_than_genome(self):
        results, elapsed = find_matches("ATCG", "ATCGAT")
        self.assertEqual(results, [])


if __name__ == "__main__":
    unittest.main()
