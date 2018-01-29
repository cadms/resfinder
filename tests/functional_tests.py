import unittest
from subprocess import Popen, PIPE, STDOUT


class ResFinderRunTest(unittest.TestCase):

    def test_on_data_with_just_acquired_resgene:
        # Maria has an E. coli isolate, with unknown resistance.
        # At first, she just wants to know which acquired resistance genes are
        # found in the genome.
        # She therefore runs resfinder cmd line

        # First Maria checks out the documentation
        procs = Popen("../run_resfinder.py -h", shell=True, stdin=PIPE,
                      stdout=PIPE, stderr=STDOUT, close_fds=True)

        output = p.stdout.read().decode()

        print("OUT: \n" + output)
