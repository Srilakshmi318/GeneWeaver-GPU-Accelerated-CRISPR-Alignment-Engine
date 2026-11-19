
import unittest
from numba import cuda

from src.gpu_engine import scan_crispr_gpu


@unittest.skipUnless(cuda.is_available(), "CUDA GPU not available")
class TestGPUEngine(unittest.TestCase):

    def test_gpu_finds_forward_and_reverse_matches(self):
        genome = "ATCGATCGATCGGATCCGATCGATCGATCG"
        guide = "ATCGAT"

        results, elapsed = scan_crispr_gpu(
            genome,
            guide,
            max_mismatches=2
        )

        strands = {result["strand"] for result in results}

        self.assertIn("+", strands)
        self.assertIn("-", strands)

    def test_gpu_returns_match_details(self):
        genome = "ATCGATCGATCGGATCCGATCGATCGATCG"
        guide = "ATCGAT"

        results, elapsed = scan_crispr_gpu(
            genome,
            guide,
            max_mismatches=2
        )

        for result in results:
            self.assertIn("position", result)
            self.assertIn("target", result)
            self.assertIn("pam", result)
            self.assertIn("mismatches", result)
            self.assertIn("penalty", result)


if __name__ == "__main__":
    unittest.main()
