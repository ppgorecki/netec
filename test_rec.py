from treeop import Tree, str2tree
from rec import rec
import pytest

ST = Tree(str2tree("(((a,b),(c,d)),e)"))
GT1 = Tree(str2tree("(c,c)"))
GT2 = Tree(str2tree("(a,a)"))
GT3 = Tree(str2tree("((b,b),e)"))
LOCATION1 = "(((a,b),(c,d)),e)"  # intervals end at root so that node is returned (location2 has also the best EC score)
LOCATION2 = "((a,b),(c,d))"      # root(gt3) is a speciation, so root(st) is not reachable


@pytest.mark.parametrize("gtrees, st, res, loc",
                         [([GT1, GT2], ST, 1, LOCATION1), ([GT1, GT3], ST, 1, LOCATION2)])
def test_rec(gtrees, st, res, loc):
    score, nodes = rec(gtrees, st)
    assert score == res
    assert len(nodes) == score
    assert str(nodes[0]) == loc
