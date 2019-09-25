from lib.blister import *
import random as rnd

Blister.Logo("Test!")

def sf(a, b): return a > b
def mf(a, b): return a == b

tree = Blister.BST(sf, mf)

for _ in range(0, 100000):
	tree.Insert(rnd.randrange(0, 10000))

print(tree.Search(100))
