import ipyrad.analysis as ipa
import toytree

def visualize_tree(path_to_file):
    """
    """
    rax = ipa.raxml(data=path_to_file, T=4, N=10)
    rax.run(block=True, force=True)
    tree = toytree.tree(rax.trees.bipartitions)
    draw = tree.root(wildcard="prz")
    draw.draw(tip_labels_align=True, node_labels="support")
