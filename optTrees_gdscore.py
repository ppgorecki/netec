# optTrees => https://github.com/agnmyk/optTrees
# gdscore => https://github.com/j-paszek/genomicduplicationr

import os
from treeop import Tree
import subprocess


def opttrees_gdscore(st: Tree, gt: Tree, outfile="optTree.out"):
    opttrees_command = f"optTrees -s '{st}' -g '{gt}' -t --maxtrees=1"
    os.system(opttrees_command)
    with open(outfile, "r") as out_f:
        lines = out_f.readlines()
        opt_tree_num = int(lines[0])
        if opt_tree_num == 1:
            gdscore_command = f"gdscore -e {outfile}"
        else:
            gdscore_command = f"gdscore - e -g '{lines[1]}' -s '{lines[-1]}'"  # multiple optimal mappings, use only the first one
        process = subprocess.Popen(gdscore_command, shell=True, stdout=subprocess.PIPE)
        (output, error) = process.communicate()
        gdscore_output = output.decode().strip().split()
        if error:
            print(error)
            raise ChildProcessError
        assert gdscore_output[:1] == ["EC", "score"]
        upper_bound = int(gdscore_output[2])
        return upper_bound
