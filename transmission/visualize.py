"""
Visualization classes for viewing data as charts, trees, etc..
"""
#Standard library imports
import os

#Third party imports
import numpy as np
import toytree
import toyplot
import toyplot.svg

class Phylo:
    """
    Methods that the user can use to visualize phylogenetic trees.
    """
    def __init__(self, newick):
        self._newick = newick
        self.tree = toytree.tree(newick)

    def get_clades(self):
        pass

    def color_clades(self):
        pass

    def _colored_heatmap_matrix(self):
        """
        Method for generating the matrix needed to make the heatmap.
        """
        pass

    def draw_tree(self, output_dir="TreeFigures", name="tree-plot", custom_style=dict()):
        """
        Draws a phylogenetic tree and saves the figure in output_dir

        Args:
            output_dir (str, optional): [description]. Defaults to "TreeFigures".
            name (str, optional): [description]. Defaults to "tree-plot".
            custom_style (dict, optional): [description]. Defaults to {}.
        """

        filename = name + '.html'

        # If custom style specified, use it
        if custom_style:
            style = custom_style

        # Otherwise, use default style
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
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        # Overwrite file if file with same name is already in the directory
        path = os.path.join(output_dir, filename)
        if os.path.exists(path):
            os.remove(path)

        # Save tree image to output_dir
        toyplot.html.render(canvas, path)

    def draw_heatmap_tree(self): 
        """
        Method for drawing tree with heatmap.
        """
        pass

#tree = Phylo("RAxML_bestTree.raxml")
#tree.draw_tree()
