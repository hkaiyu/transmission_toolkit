import toytree
import toyplot
import numpy as np


"""
Visualization classes for viewing data as charts, trees, etc..
"""

class Phylo:
    """
    Methods that the user can use to visualize phylogenetic trees.
    """
    def __init__(self, newick):
        self._newick = newick
        self.tree = toytree.tree(newick)
        
    def get_clades(self):
        pass

    def color_clades(self, cladeColors=''):
        pass

    def _colored_heatmap_matrix(self, vcf_files):
        """
        Private static method for generating the matrix needed to make the heatmap.
        """
        pass
    
    def draw_tree(self):
        pass

    def draw_heatmap_tree(self, scaleLength=2, w=500, h=350, heatmapColors="BlueRed"): #eventually add cladeColors
        """
        Method for drawing tree with heatmap.
        """
        # Create tree object
        tree = self.tree
        tree = tree.root(wildcard="prz")

        # Create canvas
        canvas = toyplot.Canvas(width=w, height=h)

        # Add tree
        axes = canvas.cartesian(bounds=(50, 150, 70, 250))
        tree.draw(
            axes=axes,
            tip_labels=False,
            tips_labels_align=True
        )

        # Add matrix
        #table = canvas.table(
            #rows= number of positions containing lfvs across all data in newick file
            #columns= tree.ntips,
            #margin=0,
            #bounds=(175,250,65,255))

        colormap = toyplot.color.brewer.map("BlueRed")

        # Apply colors to heatmap matrix (need to make it based off of allele frequencies)
        
class BB_bottleneck:
    """
    Methods to create figures regarding bottleneck size and allele frequencies.
    """
    pass

