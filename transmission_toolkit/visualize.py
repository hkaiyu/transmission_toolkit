#Standard library imports
import os
import glob
import shutil

#Third party imports
import numpy as np
import toytree
import toyplot

# Local imports
from transmission_toolkit.FASTAtools import FastaAligner, MultiFastaParser, write_fasta
from transmission_toolkit.utils import getpathleaf
from transmission_toolkit.VCFtools import extract_lfv

COLORS = [
    '#800000', '#FF0000', '#FF4500', '#FF8C00', '#FFFF00', '#7FFF00', '#008000',
    '#00FFFF', '#4169E1', '#000080', '#0000FF', '#FF00FF', '#FF1493', '#FFFAFA'
]

#vcfpath = 'test_data'
#path = 'testTransmission/parsnp/parsnp.tree'
#multifasta = 'TransmissionViz/parsnp/parsnp.mfa'

class Viz:
    """
    Class for visualizing phylogenies using toytree
    """
    def __init__(self, newick):
        tree = toytree.tree(newick, tree_format=1)
        self.tree = tree
        default_style = {
            'layout': 'r',
            'edge_type': 'p',
            'tip_labels_align': True,
            'node_labels': False,
            'node_sizes': [0 if i else 8 for i in tree.get_node_values(None, 1, 0)],
            'node_style': {
                'stroke': 'black'
        }
    }
        self.style = default_style
    
    def update(self, **kwargs):
        self.style.update(kwargs)
        
    def shorten_names(self):
        for node in self.tree.treenode.traverse():
            node.name = node.name.split('.')[0]
    
    def color_groups(self, multifasta, colors=COLORS):
        # Parse data and group records with same sequence
        data = MultiFastaParser(multifasta)
        assign_color = dict() #{name: color}
        max_idx = len(colors) - 2
        idx = 0
        for group in data.get_groups():
            if len(group) >= 2:
                for name in group:
                    assign_color[name.split('.')[0]] = colors[idx%max_idx]
                idx += 1
        
        # Traverse tree coloring nodes with same sequence
        for node in self.tree.treenode.traverse('postorder'):
            node.add_feature('color', colors[-1])
            node_name = node.name.split('.')[0]
            if node_name in assign_color:
                node.add_feature('color', assign_color[node_name])
        color_data = self.tree.get_node_values('color', show_root=1, show_tips=1)
        self.update(node_colors=color_data)

    def add_heatmap(self, vcf_dir, colors=COLORS, max_positions=15, width=1000, height=450, max_AF=1):

        # Create heatmap matrix from low frequency variants
        variant_pos = set() #keeps track of var. positions in all files
        matrix_data = dict() # maps filename to LFV data

        for fname in os.listdir(vcf_dir):
            path = os.path.join(vcf_dir, fname)
            data = extract_lfv(path, store_ref=False, parse_type='biallelic', max_AF=max_AF)
            for pos in data:
                variant_pos.add(pos)
            matrix_data[fname.split('.')[0]] = data

        variant_pos = sorted(variant_pos)

        # Create matrix and keep track of which row carries which file's data
        rows, cols = len(matrix_data), len(variant_pos)
        matrix = [[0 for j in range(cols)] for i in range(rows)]
        row_map = {idx: name.split('.')[0] for idx, name in enumerate(self.tree.get_tip_labels()[::-1])}

        # Populate matrix with data
        for i, data in enumerate(matrix_data):
            name = row_map[i]
            for j, pos in enumerate(variant_pos):
                if pos in matrix_data[name]:
                    for var in matrix_data[name][pos]:
                        matrix[i][j] = matrix_data[name][pos][var][0]
                else:
                    matrix[i][j] = 0.0

        # Get rid of columns with < 2 nonzero entries
        bad_columns = set()
        for j in range(cols):
            count = 0 # Counts nonzero entries in a column
            for i in range(rows):
                if matrix[i][j] != 0:
                    count += 1
            if count < 2:
                bad_columns.add(j)

        for i in range(rows):
            new_row = []
            for j in range(cols):
                if not j in bad_columns:
                    new_row.append(matrix[i][j])
            matrix[i] = new_row
        
        # Keep track of which column is which variant position
        position_labels = []
        for j in range(cols):
            if not j in bad_columns:
                position_labels.append(variant_pos[j])
        
        rows = len(matrix)
        matrix = [matrix[i][:max_positions] for i in range(rows)]
        cols = len(matrix[0])
        
        # create a canvas
        self.canvas = toyplot.Canvas(width=width, height=height);

        # add tree 
        tree_bounds = ('1%', '20%', '15%', '75%')
        axes = self.canvas.cartesian(bounds=tree_bounds)
        self.tree.draw(axes=axes, tip_labels=False, **self.style)
        
        
        colormap = toyplot.color.brewer.map("BlueRed", domain_min=0.0, domain_max=max_AF)

        # add matrix
        matrix_bounds = ('21%', '90%', '8%', '82%')
        tlocator = toyplot.locator.Explicit(range(cols), position_labels[:max_positions])
        table = self.canvas.matrix(
            (matrix, colormap),
            tshow=True,
            tlabel='Variant Positions',
            lshow=False,
            bounds=matrix_bounds,
            rshow=True,
            rlocator=toyplot.locator.Explicit(range(rows), self.tree.get_tip_labels()[::-1]),
            tlocator=tlocator
        )
        
        xmax_range = .95 * width
        ymax_range = .8 * height
        ymin_range = .1 * height
        
        scale = self.canvas.color_scale(
                colormap=colormap,
                x1=xmax_range,
                y1=ymax_range,
                x2=xmax_range,
                y2=ymin_range,
                width=10,
                padding=10,
                show=True,
                label="Variant Frequency",
                scale="linear",
            )

        # hide axes coordinates
        axes.show = False

    def subtree_phylogeny(self):
        pass

    def draw(self):
        self.tree.draw(**self.style)
    
    def save(self, output='', render_type='html'):
        if render_type=='svg':
            import toyplot.svg
            toyplot.svg.render(self.canvas, output)
        elif render_type=='pdf':
            import toyplot.pdf
            toyplot.pdf.render(self.canvas, output)
        elif render_type=='html':
            toyplot.html.render(self.canvas, output)
        else:
            raise ValueError('Invalid render type.')

def visualize(
    vcfdir, 
    ref, 
    colors=COLORS, 
    output_dir='TransmissionViz', 
    render_type='html', 
    threads=2,
    **filters,
):

    # Check if path exists
    if not os.path.exists(vcfdir):
        raise FileNotFoundError(f"Directory does not exist: {vcfdir}")
    if not os.path.exists(ref):
        raise FileNotFoundError(f"File does not exist: {ref}")

    # Make new directory
    os.mkdir(output_dir)

    # Generate fasta file and put those files in tmp folder
    tmp_path = os.path.join(output_dir, "tmp")
    os.mkdir(tmp_path)
    for fn in os.listdir(vcfdir):
        fp = os.path.join(vcfdir, fn)
        write_fasta(fp, ref, output_dir=tmp_path)

    # Align fasta files using parsnp and put them it in parsnp folder
    parsnp = os.path.join(output_dir, 'parsnp')
    tmpdata = FastaAligner(tmp_path)
    tmpdata.align(ref, output_dir=parsnp)

    # Get rid of tmp directory
    shutil.rmtree(tmp_path, ignore_errors=True)

    # Create Viz object that uses toytree
    newick = os.path.join(parsnp, "parsnp.tree")
    tree = Viz(newick)

    # Color nodes on tree by group
    
    # Create a heatmap matrix for the ToyTree


#visualize('example_data/test_data', 'example_data/sequence.fasta', output_dir='testTransmission')