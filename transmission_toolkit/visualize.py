"""
Module for visualizing trees using toytrees.
"""

#Standard library imports
import os
import glob
import shutil

#Third party imports
import toytree
import toyplot

# Local imports
from transmission_toolkit.FASTAtools import FastaAligner, MultiFastaParser
from transmission_toolkit.convertFile import vcf2fasta, root_newick
from transmission_toolkit.utils import getpathleaf
from transmission_toolkit.VCFtools import extract_lfv

COLORS = [
    '#DC050C', '#E8601C', '#F1932D', '#F6C141', '#F7F056', '#CAE0AB',
    '#90C987', '#4EB265', '#7BAFDE', '#5289C7', '#1965B0', '#882E72'
]

class PhyloTree:
    """
    Class for visualizing phylogenies using toytree
    """
    def __init__(self, newick, colors=COLORS):
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
        self.colors = colors
    
    def update(self, **kwargs):
        static_attributes = {'layout', 'node_sizes'}
        for key in kwargs.keys():
            if key in static_attributes:
                raise ValueError(f'Cannot change attribute: {key}.')
        self.style.update(kwargs)
    
    def color_groups(self, multifasta):
        # Parse data and group records with same sequence
        data = MultiFastaParser(multifasta)
        assign_color = dict() #{name: color}
        ref = data.records[0].name
        max_idx = len(self.colors) - 1
        idx = 0
        for group in data.get_groups():
            if ref in group:
                group.remove(ref)
            if len(group) >= 2:
                for name in group:
                    assign_color[name.split('.')[0]] = self.colors[idx%max_idx]
                idx += 1
        self.group_colors = assign_color
        
        # Traverse tree coloring nodes with same sequence
        for node in self.tree.treenode.traverse('postorder'):
            node.add_feature('color', '#FFFAFA') # default color white
            node_name = node.name.split('.')[0]
            if node_name in assign_color:
                node.add_feature('color', assign_color[node_name])
        color_data = self.tree.get_node_values('color', show_root=1, show_tips=1)
        self.update(node_colors=color_data)

    def add_heatmap(
        self, 
        vcf_dir,  
        position_range=None, 
        width=1000, 
        height=450, 
        min_read_depth=1, 
        store_ref=True, 
        masks=None, 
        mask_status='hide',
        filter_columns=True
    ):
        if position_range and not isinstance(position_range, tuple):
            raise TypeError("position_range parameter must be a tuple.")

        # Create heatmap matrix from low frequency variants
        variant_pos = set() #keeps track of var. positions in all files
        matrix_data = dict() # maps filename to LFV data

        for fname in os.listdir(vcf_dir):
            path = os.path.join(vcf_dir, fname)
            data = extract_lfv(
                path, 
                min_read_depth=min_read_depth, 
                max_AF=1,
                parse_type='biallelic',
                store_ref=store_ref, 
                masks=masks, 
                mask_status=mask_status
            )
            for pos in data:
                variant_pos.add(pos)
            matrix_data[fname.split('.')[0]] = data

        variant_pos = sorted(variant_pos)

        # Create matrix and keep track of which row carries which file's data
        rows, cols = len(matrix_data), len(variant_pos)
        matrix = [[0 for j in range(cols)] for i in range(rows)]
        row_map = {idx: name.split('.')[0] for idx, name in enumerate(self.tree.get_tip_labels()[::-1])}
        column2position = {idx: pos for idx, pos in enumerate(variant_pos)}

        # Populate matrix with data
        for i, data in enumerate(matrix_data):
            name = row_map[i]
            for j, pos in column2position.items():
                if pos in matrix_data[name]:
                    for var in matrix_data[name][pos]:
                        matrix[i][j] = matrix_data[name][pos][var][0]
                else:
                    matrix[i][j] = 0.0
        
        # Keep track of columns that should be removed
        bad_columns = set()
        for j in range(cols):
            
            # filters out columns with less than 2 nonzero frequencies
            if filter_columns:
                count = 0 
                for i in range(rows):
                    if matrix[i][j] != 0:
                        count += 1
                if count < 2:
                    bad_columns.add(j)
                    
            # filters out columns not within range
            if position_range:
                if not position_range[0] <= column2position[j] <= position_range[1]:
                    bad_columns.add(j)
        
        # Label positions on heatmap
        position_labels = []
        for j in range(cols):
            if not j in bad_columns:
                position_labels.append(column2position[j])
        
        # Modify matrix to not include filtered columns
        for i in range(rows):
            new_row = []
            for j in range(cols):
                if not j in bad_columns:
                    new_row.append(matrix[i][j])

            matrix[i] = new_row
    
        rows = len(matrix)
        cols = len(matrix[0])
        if cols == 0:
            raise IOError('There are no variants to make a heatmap in given VCF files.')
        
        # create a canvas
        self.canvas = toyplot.Canvas(width=width, height=height);

        # add tree 
        tree_bounds = ('1%', '20%', '15%', '75%')
        axes = self.canvas.cartesian(bounds=tree_bounds)
        self.tree.draw(axes=axes, tip_labels=False, **self.style)
        
        
        colormap = toyplot.color.brewer.map("BlueRed", domain_min=0, domain_max=1)

        # add matrix
        matrix_bounds = ('21%', '88%', '6%', '84%')
        tlocator = toyplot.locator.Explicit(range(cols), position_labels)
        rlocator = toyplot.locator.Explicit(range(rows), self.tree.get_tip_labels()[::-1])
        table = self.canvas.matrix(
            (matrix, colormap),
            tshow=True,
            tlabel='Variant Positions',
            lshow=False,
            bounds=matrix_bounds,
            rshow=True,
            rlocator=rlocator,
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

    def draw(self):
        self.tree.draw(**self.style)
    
    def save(self, output, render_type='html'):
        if render_type=='svg':
            import toyplot.svg
            toyplot.svg.render(self.canvas, output)
        elif render_type=='pdf':
            import toyplot.pdf
            toyplot.pdf.render(self.canvas, output)
        elif render_type=='html':
            import toyplot
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
    min_read_depth=1,
    parse_type='biallelic',
    store_ref=True,
    masks=None,
    mask_status='hide'
    ):
    '''
    Function that generates figures on the full tree and subtree phylogenies given
    a directory of VCF files and a reference genome.

    It should be noted that running this provides less customizablity than using the
    PhyloTree class methods themselves.
    '''

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
        vcf2fasta(
            fp, 
            ref, 
            output_dir=tmp_path,
            min_read_depth=min_read_depth,
            max_AF=1,
            parse_type=parse_type,
            store_ref=store_ref,
            masks=masks,
            mask_status=mask_status,
        )

    # Align fasta files using parsnp and put them it in parsnp folder
    parsnp = os.path.join(output_dir, 'parsnp')
    tmpdata = FastaAligner(tmp_path)
    tmpdata.align(ref, output_dir=parsnp)
    newick = os.path.join(parsnp, 'parsnp.tree')
    root_newick(newick, os.path.join(parsnp, 'rparsnp.tree'))
    rnewick = os.path.join(parsnp, 'rparsnp.tree')

    # Run normal tree stuff or whatever
    multifasta = os.path.join(parsnp, 'parsnp.mfa')
    tree = PhyloTree(rnewick)
    tree.color_groups(multifasta)
    tree.draw()
    tree.add_heatmap(vcfdir, height=400, width=1000)
    tree.save(os.path.join(output_dir, 'full_phylogeny.' + render_type))

    # store colors assigned
    assigned_colors = tree.group_colors

    # Group vcf files into list of groups
    data = MultiFastaParser(multifasta)
    old_groups = list(data.get_groups())
    groups = list()
    refname = data.records[0].name
    for group in old_groups:
        if refname in group:
            group.remove(refname)
        if len(group) >= 2:
            groups.append(group)
            
    dir_map = {name.split('.')[0]: name for name in os.listdir(vcfdir)}
    i = 1
    for group in groups:
        group_dir = os.path.join(output_dir, f'group_{i}')
        group_tmp = os.path.join(group_dir, 'tmp')
        tmp_vcf = os.path.join(group_tmp, 'vcf')
        os.mkdir(group_dir)
        os.mkdir(group_tmp)
        os.mkdir(tmp_vcf)
        for name in group:
            splitname = name.split('.')[0]
            color = assigned_colors[splitname]
            vcf = dir_map[splitname]
            vcfpath = os.path.join(vcfdir, vcf)
            shutil.copy(vcfpath, tmp_vcf)
            vcf2fasta(
                vcfpath, 
                ref, 
                output_dir=group_tmp, 
                min_read_depth=min_read_depth, 
                max_AF=1, 
                parse_type=parse_type, 
                store_ref=store_ref, 
                masks=masks, 
                mask_status=mask_status,
                consensus_type='minor'
            )

        parsnp = os.path.join(group_dir, 'parsnp')
        tmpdata = FastaAligner(group_tmp)
        tmpdata.align(ref, output_dir=parsnp)

        newick = os.path.join(parsnp, 'parsnp.tree')
        multifasta = os.path.join(parsnp, 'parsnp.mfa')

        tree = PhyloTree(newick)
        tree.update(node_colors=[color])
        tree.draw()
        i += 1
        try:
            tree.add_heatmap(tmp_vcf, height=400, width=1000, filter_columns=False)
        except:
            print('Cannot make heatmap... continue!')
            continue
        tree.save(os.path.join(output_dir, f'subtree_heatmap{i}.' + render_type))
        shutil.rmtree(group_tmp, ignore_errors=True)

    # Get rid of tmp directory
    shutil.rmtree(tmp_path, ignore_errors=True)

#visualize('example_data/Denmark_WGS_47', 'example_data/sequence.fasta', output_dir='demo')