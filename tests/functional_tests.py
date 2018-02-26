import unittest
from subprocess import PIPE, run
import os
import shutil


test_names = ["test1", "test3"]
test_data = {
    test_names[0]: "data/test_isolate_01.fa",  # Test published resistance
    test_names[1]: "data/test_isolate_03.fa",  # Test no resistance
}
run_test_dir = "running_test"
blast_out_dir = "blast_out"


class ResFinderRunTest(unittest.TestCase):

    def setUp(self):
        # Change working dir to test dir
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        # Does not allow running two tests in parallel
        os.makedirs(run_test_dir, exist_ok=False)

    def tearDown(self):
        # shutil.rmtree(run_test_dir)
        pass

    def test_on_data_with_just_acquired_resgene(self):
        # Maria has an E. coli isolate, with unknown resistance.
        # At first, she just wants to know which acquired resistance genes are
        # found in the genome.
        # She therefore runs resfinder cmd line.

        # First Maria checks out the documentation.
        procs = run("python3 ../run_resfinder.py -h", shell=True, stdout=PIPE,
                    check=True)
        output = procs.stdout.decode()
        self.assertIn("--help", output)

        # Maria goes on to run ResFinder for acquired genes with her E. coli
        # isolate.

        # First she creates a few directories to store her output.
        test1_dir = run_test_dir + "/" + test_names[0]
        test1_blast_dir = test1_dir + "/" + blast_out_dir
        os.makedirs(test1_dir)
        os.makedirs(test1_blast_dir)

        # Then she runs run_resfinder with her first isolate.
        cmd_acquired = ("python3 ../run_resfinder.py"
                        + " -i " + test_data[test_names[0]]
                        + " -o " + test1_blast_dir
                        + " -s e.coli"
                        + " --min_cov 0.6"
                        + " -t 0.8"
                        + " --acquired"
                        + " --databasePath_res ../database")

        procs = run(cmd_acquired, shell=True, stdout=PIPE, stderr=PIPE,
                    check=True)

        # Expected output files
        res_out = test1_blast_dir + "/resfinder_out"

        fsa_hit = res_out + "/Hit_in_genome_seq.fsa"
        fsa_res = res_out + "/Resistance_gene_seq.fsa"
        res_table = res_out + "/results_table.txt"
        res_tab = res_out + "/results_tab.txt"
        results = res_out + "/results.txt"

        with open(fsa_hit, "r") as fh:
            check_result = fh.readline()
        self.assertIn("blaB-2_1_AF189300", check_result)

        with open(fsa_res, "r") as fh:
            check_result = fh.readline()
        self.assertIn("blaB-2_AF189300", check_result)

        with open(res_table, "r") as fh:
            for i, line in enumerate(fh):
                check_result = fh.readline()
                if(i == 17):
                    print("RES: " + check_result)
                    break
        self.assertIn("blaB-2_1_AF189300", check_result)

        with open(res_tab, "r") as fh:
            fh.readline()
            check_result = fh.readline()
        self.assertIn("blaB-2_1_AF189300", check_result)

        with open(results, "r") as fh:
            fh.readline()
            fh.readline()
            fh.readline()
            fh.readline()
            fh.readline()
            check_result = fh.readline()
        self.assertIn("blaB-2_1_AF189300", check_result)

        # self.fail("Finish the test for pointmutations")
        # self.fail("Finish the test for phenotypes")


if __name__ == "__main__":
    unittest.main()
