import unittest
from subprocess import PIPE, run
import os


class ResFinderRunTest(unittest.TestCase):
    test_data = {
        "test1": "data/test_isolate_01.fa",  # Test published resistance
        "test3": "data/test_isolate_03.fa",  # Test no resistance
    }

    def setUp(self):
        os.makedirs("running_test", exist_ok=True)

    def test_on_data_with_just_acquired_resgene(self):
        # Maria has an E. coli isolate, with unknown resistance.
        # At first, she just wants to know which acquired resistance genes are
        # found in the genome.
        # She therefore runs resfinder cmd line

        # First Maria checks out the documentation
        procs = run("python3 ../run_resfinder.py -h", shell=True, stdout=PIPE,
                    check=True)
        output = procs.stdout.decode()
        self.assertIn("--help", output)

        # Maria goes on to run ResFinder for acquired genes with her E. coli
        # isolate
#        cmd_acquired = ("python3 ../run_resfinder.py"
#                        " -i " + test_data["test1"]
#                        " ")


if __name__ == "__main__":
    unittest.main()
