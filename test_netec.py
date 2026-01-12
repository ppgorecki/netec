#!/usr/bin/env python3
"""
Test suite for netec.py based on README.md examples

This comprehensive test suite validates netec.py functionality by running
tests based on the examples provided in README.md.

Usage:
    python3 test_netec.py

Requirements:
    - netec.py must be in the current directory
    - data_sim/ directory with test data must exist
    - data_yeast/ directory (optional, test skipped if missing)

Test Coverage:
    - Basic usage and output options
    - Verbose modes and statistics
    - Episode analysis features
    - Network with reticulations
    - Distribution maps
    - Performance optimizations
    - Various output file formats
    - GSE-style output
    - Combined options

Exit Codes:
    0 - All tests passed
    1 - One or more tests failed
"""

import os
import sys
import tempfile
import shutil
import subprocess
from pathlib import Path

class Colors:
    """ANSI color codes for terminal output"""
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    RESET = '\033[0m'
    BOLD = '\033[1m'

class NetECTester:
    def __init__(self):
        self.passed = 0
        self.failed = 0
        self.test_dir = None

    def setup(self):
        """Create temporary test directory"""
        self.test_dir = tempfile.mkdtemp(prefix="netec_test_")
        print(f"{Colors.BLUE}Test directory: {self.test_dir}{Colors.RESET}\n")

    def cleanup(self):
        """Remove temporary test directory"""
        if self.test_dir and os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    def run_command(self, cmd, description=""):
        """Run a shell command and return result"""
        try:
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                timeout=60
            )
            return result.returncode == 0, result.stdout, result.stderr
        except subprocess.TimeoutExpired:
            return False, "", "Command timed out"
        except Exception as e:
            return False, "", str(e)

    def test(self, name, cmd, check_output=None, check_files=None, should_fail=False):
        """Run a single test"""
        print(f"{Colors.BOLD}Test: {name}{Colors.RESET}")
        print(f"Command: {cmd}")

        success, stdout, stderr = self.run_command(cmd)

        # Check if command succeeded/failed as expected
        if should_fail:
            if success:
                print(f"{Colors.RED}✗ FAILED: Expected command to fail but it succeeded{Colors.RESET}\n")
                self.failed += 1
                return False
        else:
            if not success:
                print(f"{Colors.RED}✗ FAILED: Command returned non-zero exit code{Colors.RESET}")
                if stderr:
                    print(f"Error: {stderr[:200]}")
                print()
                self.failed += 1
                return False

        # Check output contains expected strings
        if check_output:
            for expected in check_output:
                if expected not in stdout:
                    print(f"{Colors.RED}✗ FAILED: Expected '{expected}' in output{Colors.RESET}")
                    print(f"Got: {stdout[:300]}...")
                    print()
                    self.failed += 1
                    return False

        # Check that expected files were created
        if check_files:
            for filepath in check_files:
                if not os.path.exists(filepath):
                    print(f"{Colors.RED}✗ FAILED: Expected file not created: {filepath}{Colors.RESET}\n")
                    self.failed += 1
                    return False

        print(f"{Colors.GREEN}✓ PASSED{Colors.RESET}\n")
        self.passed += 1
        return True

    def run_all_tests(self):
        """Run all test cases"""
        print(f"{Colors.BOLD}{Colors.BLUE}{'='*60}")
        print("NetEC Test Suite")
        print(f"{'='*60}{Colors.RESET}\n")

        self.setup()

        try:
            # Test 1: Basic usage - minimal example
            self.test(
                "Basic usage - minimal example",
                "python3 netec.py --gene_trees data_sim/wgd-1-gene-trees --network data_sim/s_tree",
                check_output=["Cost:", "Exact:"]
            )

            # Test 2: With output file
            out_file = os.path.join(self.test_dir, "output.log")
            self.test(
                "With output file",
                f"python3 netec.py --gene_trees data_sim/wgd-1-gene-trees --network data_sim/s_tree --out_file {out_file}",
                check_files=[out_file]
            )

            # Test 3: With output directory
            out_dir = os.path.join(self.test_dir, "results")
            os.makedirs(out_dir, exist_ok=True)
            self.test(
                "With output directory",
                f"python3 netec.py --gene_trees data_sim/wgd-1-gene-trees --network data_sim/s_tree --out_file {out_dir}",
                check_files=[os.path.join(out_dir, "netec.log")]
            )

            # Test 4: Verbose output
            self.test(
                "Verbose output (level 1)",
                "python3 netec.py --gene_trees data_sim/wgd-1-gene-trees --network data_sim/s_tree --verbose 1"
            )

            # Test 5: Verbose output level 2
            self.test(
                "Verbose output (level 2)",
                "python3 netec.py --gene_trees data_sim/wgd-1-gene-trees --network data_sim/s_tree --verbose 2"
            )

            # Test 6: Print duplication statistics
            self.test(
                "Print duplication statistics",
                "python3 netec.py --gene_trees data_sim/wgd-1-gene-trees --network data_sim/s_tree --print_dup_stats",
                check_output=["Duplication statistics:"]
            )

            # Test 7: Minimum duplication filter
            self.test(
                "Filter trees with --mindup",
                "python3 netec.py --gene_trees data_sim/wgd-1-gene-trees --network data_sim/s_tree --mindup 2"
            )

            # Test 8: Episode summary file
            out_dir2 = os.path.join(self.test_dir, "episode_test")
            os.makedirs(out_dir2, exist_ok=True)
            self.test(
                "Save episode summary",
                f"python3 netec.py --gene_trees data_sim/wgd-1-gene-trees --network data_sim/s_tree --out_file {out_dir2} --episummaryfile",
                check_files=[os.path.join(out_dir2, "episummary")]
            )

            # Test 9: Save embedding
            out_dir3 = os.path.join(self.test_dir, "embedding_test")
            os.makedirs(out_dir3, exist_ok=True)
            self.test(
                "Save embedding",
                f"python3 netec.py --gene_trees data_sim/wgd-1-gene-trees --network data_sim/s_tree --out_file {out_dir3} --save_embedding",
                check_files=[os.path.join(out_dir3, "embedding")]
            )

            # Test 10: Locked episode support (use directory for simpler file naming)
            out_dir4 = os.path.join(self.test_dir, "locked_epi_test")
            os.makedirs(out_dir4, exist_ok=True)
            self.test(
                "Locked episode support analysis",
                f"python3 netec.py --gene_trees data_sim/wgd-1-gene-trees --network data_sim/s_tree --locked_epi_support --out_file {out_dir4}",
                check_files=[
                    os.path.join(out_dir4, "locked_epi"),
                    os.path.join(out_dir4, "locked_epi_net")
                ]
            )

            # Test 11: User-defined episodes (all)
            self.test(
                "User episodes - all nodes",
                "python3 netec.py --gene_trees data_sim/wgd-1-gene-trees --network data_sim/s_tree --user_episodes all --verbose 1"
            )

            # Test 12: Randomization options
            self.test(
                "Performance optimization with randomization",
                "python3 netec.py --gene_trees data_sim/wgd-1-gene-trees --network data_sim/s_tree --randomize_from 1000 --noimprovement_stop 10"
            )

            # Test 13: Reversed climb search
            self.test(
                "Reversed climb search strategy",
                "python3 netec.py --gene_trees data_sim/wgd-1-gene-trees --network data_sim/s_tree --reversed_climb"
            )

            # Test 14: GSE-style output
            out_dir5 = os.path.join(self.test_dir, "gse_test")
            os.makedirs(out_dir5, exist_ok=True)
            self.test(
                "GSE-style output format",
                f"python3 netec.py --gene_trees data_sim/wgd-1-gene-trees --network data_sim/s_tree --out_file {out_dir5} --gse"
            )

            # Test 15: Create network with reticulations and test
            network_file = os.path.join(self.test_dir, "network.txt")
            gtrees_file = os.path.join(self.test_dir, "gtrees.txt")

            with open(network_file, 'w') as f:
                f.write("((c)#A,((#A,b),a))\n")

            with open(gtrees_file, 'w') as f:
                f.write("((a,a),(b,c))\n")
                f.write("((a,b),(c,c))\n")

            self.test(
                "Network with reticulations",
                f"python3 netec.py --gene_trees {gtrees_file} --network {network_file}"
            )

            # Test 16: Distribution maps (if data_sim is a tree without reticulations)
            self.test(
                "Distribution maps",
                "python3 netec.py --gene_trees data_sim/wgd-1-gene-trees --network data_sim/s_tree --distribution_maps --print_distr_maps"
            )

            # Test 17: Test with yeast data if available
            if os.path.exists("data_yeast/gtrees") and os.path.exists("data_yeast/s_tree"):
                out_dir6 = os.path.join(self.test_dir, "yeast_test")
                os.makedirs(out_dir6, exist_ok=True)
                self.test(
                    "Real data - Yeast dataset",
                    f"python3 netec.py --gene_trees data_yeast/gtrees --network data_yeast/s_tree --out_file {out_dir6} --verbose 1"
                )
            else:
                print(f"{Colors.YELLOW}Skipping yeast test - data not available{Colors.RESET}\n")

            # Test 18: Test with small gene tree file
            small_gtrees = os.path.join(self.test_dir, "small_gtrees.txt")
            with open(small_gtrees, 'w') as f:
                f.write("((1,2),(3,4))\n")
                f.write("((1,1),(2,3))\n")
                f.write("(((1,1),2),(3,4))\n")

            small_stree = os.path.join(self.test_dir, "small_stree.txt")
            with open(small_stree, 'w') as f:
                f.write("((1,2),(3,4))\n")

            self.test(
                "Small custom dataset",
                f"python3 netec.py --gene_trees {small_gtrees} --network {small_stree}"
            )

            # Test 19: WGD debug option
            self.test(
                "WGD debug mode",
                f"python3 netec.py --gene_trees {small_gtrees} --network {small_stree} --wgddebug"
            )

            # Test 20: Combined options (note: print_dup_stats causes early return, so excluded)
            out_dir7 = os.path.join(self.test_dir, "combined_test")
            os.makedirs(out_dir7, exist_ok=True)
            self.test(
                "Combined options test",
                f"python3 netec.py --gene_trees data_sim/wgd-1-gene-trees --network data_sim/s_tree "
                f"--out_file {out_dir7} --verbose 1 --episummaryfile --save_embedding",
                check_files=[
                    os.path.join(out_dir7, "netec.log"),
                    os.path.join(out_dir7, "episummary"),
                    os.path.join(out_dir7, "embedding")
                ]
            )

        finally:
            self.cleanup()

        # Print summary
        print(f"{Colors.BOLD}{Colors.BLUE}{'='*60}")
        print("Test Summary")
        print(f"{'='*60}{Colors.RESET}")
        print(f"{Colors.GREEN}Passed: {self.passed}{Colors.RESET}")
        if self.failed > 0:
            print(f"{Colors.RED}Failed: {self.failed}{Colors.RESET}")
        else:
            print(f"Failed: {self.failed}")

        total = self.passed + self.failed
        if total > 0:
            success_rate = (self.passed / total) * 100
            print(f"Success rate: {success_rate:.1f}%")

        print(f"{Colors.BLUE}{'='*60}{Colors.RESET}\n")

        return self.failed == 0

def main():
    # Check if we're in the right directory
    if not os.path.exists("netec.py"):
        print(f"{Colors.RED}Error: netec.py not found. Run this script from the netec directory.{Colors.RESET}")
        sys.exit(1)

    if not os.path.exists("data_sim"):
        print(f"{Colors.RED}Error: data_sim directory not found.{Colors.RESET}")
        sys.exit(1)

    tester = NetECTester()
    success = tester.run_all_tests()

    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
