import sys, os
from ete3 import Tree, TreeStyle, random_color, PhyloTree


if len(sys.argv) < 3:
	print("Missing arguments: correct usage: " + sys.argv[0] + "clust2.newick out_dir")
	exit(1)

if os.path.isfile(sys.argv[1]) == False:
	print('FileNotFoundError: No such file or directory:' + sys.argv[1])
	exit(1)

if os.path.isdir(sys.argv[2]) == False:
	print('FileNotFoundError: No such file or directory:' + sys.argv[2])
	exit(1)

def clean_ly(node):
    node.img_style['size'] = 0

t = Tree(sys.argv[1])

ts = TreeStyle()

# set a clean view for the tree
ts.show_leaf_name = True
ts.mode = 'c'
ts.layout_fn = clean_ly


# fix a reasonable size for the image, so tree is scaled to it
t.render(sys.argv[2]+'/PhyloTree.png', w=3000, units='px', tree_style=ts)


#################################################################
#           Obtain all the evolutionary events                  #
#################################################################

#pt = PhyloTree('clust2.newick')
#events = pt.get_descendant_evol_events()
#print("Events detected from the root of the tree")
#for ev in events:
#    if ev.etype == "S":
#        print('   ORTHOLOGY RELATIONSHIP:', ','.join(ev.in_seqs), "<====>", ','.join(ev.out_seqs))
#    elif ev.etype == "D":
#        print('   PARALOGY RELATIONSHIP:', ','.join(ev.in_seqs), "<====>", ','.join(ev.out_seqs))
