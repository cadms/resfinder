import unittest
from subprocess import PIPE, run
import os
import shutil
import sys


# This is not best practice but for testing, this is the best I could
# come up with
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))


test_names = ["test1", "test2", "test3"]
test_data = {
    # Test published acquired resistance
    test_names[0]: "data/test_isolate_01.fa",
    test_names[1]: "data/test_isolate_01_1.fq data/test_isolate_01_2.fq",
    # Test published point mut resistance
    test_names[2]: "data/test_isolate_05.fa",
}
run_test_dir = "running_test"


class ResFinderRunTest(unittest.TestCase):

    def setUp(self):
        # Change working dir to test dir
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        # Does not allow running two tests in parallel
        os.makedirs(run_test_dir, exist_ok=False)

    def tearDown(self):
        # shutil.rmtree(run_test_dir)
        pass

    def p_test_on_data_with_just_acquired_resgene_using_blast(self):
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
        os.makedirs(test1_dir)

        # Then she runs run_resfinder with her first isolate.
        cmd_acquired = ("python3 ../run_resfinder.py"
                        + " -ifa " + test_data[test_names[0]]
                        + " -o " + test1_dir
                        + " -s e.coli"
                        + " --min_cov 0.6"
                        + " -t 0.8"
                        + " --acquired"
                        + " --databasePath_res ../database")

        procs = run(cmd_acquired, shell=True, stdout=PIPE, stderr=PIPE,
                    check=True)

        # Expected output files
        res_out = test1_dir + "/resfinder_blast"

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
            for line in fh:
                if(line.startswith("blaB-2")):
                    check_result = line
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

    def p_test_on_data_with_just_acquired_resgene_using_kma(self):
        # Maria has another E. coli isolate, with unknown resistance.
        # This time she does not have an assembly, but only raw data.
        # She therefore runs resfinder cmd line using KMA.

        # First she creates a few directories to store her output.
        test2_dir = run_test_dir + "/" + test_names[1]
        os.makedirs(test2_dir, exist_ok=False)

        # Then she runs run_resfinder with her first isolate.
        cmd_acquired = ("python3 ../run_resfinder.py"
                        + " -ifq " + test_data[test_names[1]]
                        + " -o " + test2_dir
                        + " -s e.coli"
                        + " --min_cov 0.6"
                        + " -t 0.8"
                        + " --acquired"
                        + " --databasePath_res ../database"
                        + " --kmaPath ../cge/kma/kma")

        procs = run(cmd_acquired, shell=True, stdout=PIPE, stderr=PIPE,
                    check=True)

        # Expected output files
        res_out = test2_dir + "/resfinder_kma"

        fsa_hit = res_out + "/Hit_in_genome_seq.fsa"
        fsa_res = res_out + "/Resistance_gene_seq.fsa"
        res_table = res_out + "/results_table.txt"
        res_tab = res_out + "/results_tab.txt"
        results = res_out + "/results.txt"

        with open(fsa_hit, "r") as fh:
            check_result = fh.readline()
        self.assertIn("blaB-2", check_result)

        with open(fsa_res, "r") as fh:
            check_result = fh.readline()
        self.assertIn("blaB-2_AF189300", check_result)

        with open(res_table, "r") as fh:
            for line in fh:
                if(line.startswith("blaB-2")):
                    check_result = line
                    break
        self.assertIn("blaB-2", check_result)

        with open(res_tab, "r") as fh:
            fh.readline()
            check_result = fh.readline()
        self.assertIn("blaB-2", check_result)

        with open(results, "r") as fh:
            fh.readline()
            fh.readline()
            fh.readline()
            fh.readline()
            fh.readline()
            check_result = fh.readline()
        self.assertIn("blaB-2", check_result)

    def test_on_data_with_just_point_mut_using_blast(self):
        # Maria also wants to check her assembled E. coli isolate for
        # resistance caused by point mutations.

        # First she creates a few directories to store her output.
        test3_dir = run_test_dir + "/" + test_names[2]
        os.makedirs(test3_dir)

        # Then she runs run_resfinder with her first isolate.
        cmd_point = ("python3 ../run_resfinder.py"
                     + " -ifa " + test_data[test_names[0]]
                     + " -o " + test3_dir
                     + " -s e.coli"
                     + " --min_cov 0.6"
                     + " --threshold 0.8"
                     + " --point"
                     + " --databasePath_point ../database_pointfinder")

        procs = run(cmd_point, shell=True, stdout=PIPE, stderr=PIPE,
                    check=True)

        # Expected output files
        pf_pred = test3_dir + "/PointFinder_prediction.txt"
        pf_res = test3_dir + "/PointFinder_results.txt"
        pf_table = test3_dir + "/PointFinder_table.txt"

        with open(pf_pred, "r") as fh:
            fh.readline()
            fh.readline()
            fh.readline()
            pred_line = fh.readline()
        pred_lst = pred_line.split()
        print("LIST: " + str(pred_lst))
        self.assertEqual(1 == pred_lst[14])
        self.assertEqual(1 == pred_lst[15])

        with open(pf_res, "r") as fh:
            fh.readline()
            check_result = fh.readline()
        self.assertIn("gyrA", check_result)
        self.assertIn("p.S83A", check_result)

        point_mut_found = False
        with open(res_table, "r") as fh:
            for line in fh:
                if(line.startswith("gyrA p.S83A")):
                    check_result = line
                    point_mut_found = True
                    break
        self.assertEqual(point_mut_found is True)


if __name__ == "__main__":
    unittest.main()
