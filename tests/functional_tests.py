import unittest
from subprocess import Popen, PIPE, STDOUT, run


class ResFinderRunTest(unittest.TestCase):

    def test_on_data_with_just_acquired_resgene(self):
        # Maria has an E. coli isolate, with unknown resistance.
        # At first, she just wants to know which acquired resistance genes are
        # found in the genome.
        # She therefore runs resfinder cmd line

        # First Maria checks out the documentation
        procs = run("python3 ../run_resfinder.py -h", shell=True, stdout=PIPE, stderr=PIPE)

        # output = procs.stdout.read().decode()
        # outs, errs = procs.communicate()

        print("ERR: \n" + procs.stderr.decode())


if __name__ == "__main__":
    unittest.main()