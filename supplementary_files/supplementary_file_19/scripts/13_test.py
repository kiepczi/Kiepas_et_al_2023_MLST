from ete3 import Tree

# Load your phylogenetic tree
tree = Tree("../output/tree/04_tbe_mlst.raxml.support")

# Initialize a counter for branches with bootstrap values below 50
count_bs_below_50 = 0

# Traverse the tree and check bootstrap values
for node in tree.traverse():
    if node.support == 1.0:
        print(node.support)
        count_bs_below_50 += 1

# Print the result
print(f"Number of branches with bootstrap values below 50: {count_bs_below_50}")
