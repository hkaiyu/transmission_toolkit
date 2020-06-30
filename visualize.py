import toytree
import toyplot
import numpy as np
import os
import toyplot.svg

"""
Visualization classes for viewing data as charts, trees, etc..
"""

class Phylo:
    """
    Methods that the user can use to visualize phylogenetic trees.
    """
    def __init__(self, newick):
        self._newick = newick
        self.tree = toytree.tree(newick, tree_format=0)
        
    def get_clades(self):
        pass

    def color_clades(self, cladeColors=''):
        pass

    def _colored_heatmap_matrix(self, vcf_files):
        """
        Method for generating the matrix needed to make the heatmap.
        """
        pass
    
    def draw_tree(self, outputDir="TreeFigures", name="tree-plot", customStyle={}):
        """
        Creates standard phylogenetic tree visual using ToyTree and writes the image to outputDir.

        Args:
            orient (str, optional): Orientation of tree (ex. tree pointing up, down, right, etc.). Defaults to 'down'.
            render_type ([type], optional): The type of rendering for the image. Defaults to html.
            outputDir (str, optional): The name of the directory that this method will write the image to. Defaults to "TreeFigures".
            name (str, optional): The name of the file. Defaults to "tree-plot".
        """

        filename = name + '.html'

        # If custom style specified, use it
        if customStyle:
            style = customStyle

        # Otherwise, use default style (figure out good default style, also be able to update default style with custom style arguments)
        else:
            style = {
                "height": 875,
                "width": 1400,
                "layout": 'd',
                "edge_type": 'p',
                "edge_style": {
                    "stroke": toytree.colors[2],
                    "stroke-width": 2.5},
                "tip_labels_align": True,
                "tip_labels_style": {
                    "font-size": "10px"},
                "node_labels": False,
                "node_sizes": 7,
                "node_colors": toytree.colors[1]}
        

        # Create tree drawing
        canvas = self.tree.draw(**style)[0]

        # Make directory if does not already exist
        if not os.path.exists(outputDir):
            os.mkdir(outputDir)
        
        # Overwrite file if file with same name is already in the directory
        path = os.path.join(outputDir, filename)
        if os.path.exists(path):
            os.remove(path)

        # Save tree image to outputDir
        toyplot.html.render(canvas, path)


    def draw_heatmap_tree(self, scaleLength=2, w=500, h=350, heatmapColors="BlueRed"): #eventually add cladeColors
        """
        Method for drawing tree with heatmap.
        """
        # Create tree object
        tree = self.tree

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

#tree = Phylo("RAxML_bestTree.raxml")
#tree.draw_tree()
