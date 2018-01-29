import unittest
from subprocess import Popen, PIPE, STDOUT


class ResFinderRunTest(unittest.TestCase):

    def test_on_data_with_just_acquired_resgene(self):
        # Maria has an E. coli isolate, with unknown resistance.
        # At first, she just wants to know which acquired resistance genes are
        # found in the genome.
        # She therefore runs resfinder cmd line

        # First Maria checks out the documentation
        procs = Popen("python3 ../run_resfinder.py -h", shell=True) #stdin=PIPE,
                      # stdout=PIPE, stderr=PIPE)

        # output = procs.stdout.read().decode()
        outs, errs = procs.communicate()

        print("ERR: \n" + errs.decode())


if __name__ == "__main__":
    unittest.main()
