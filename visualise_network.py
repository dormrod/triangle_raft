import sys
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib._color_data as mcd
import matplotlib.colors as colors
import matplotlib.pylab as pylab
import numpy as np

def main():

    # Get options for visualisation and unpack
    vis_options=get_options()
    vis_prefix=vis_options.get("prefix")
    vis_tri=vis_options.get("triangle",False)
    vis_x=vis_options.get("x_atom",False)
    vis_m=vis_options.get("m_atom",False)
    vis_ring=vis_options.get("ring",False)
    vis_atom_label=vis_options.get("atom_label",False)
    vis_tri_label=vis_options.get("triangle_label",False)
    # vis_atom_label=vis_options.get("atom_label",False)
    # vis_ring_label=vis_options.get("ring_label",False)
    vis_save_pdf=vis_options.get("save_pdf",False)
    vis_save_png=vis_options.get("save_png",False)

    # Set up axes
    updateParams()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.axis('off')

    # Get data and plot
    data=get_data(vis_prefix)
    if vis_tri:
        plot_triangle_network(data,vis_x,vis_m,vis_atom_label,vis_tri_label,fig,ax)
    if vis_ring:
        plot_ring_network(data,fig,ax)

    # Save and Show plot
    if vis_save_pdf: savePlot(vis_prefix,fmt="pdf")
    if vis_save_png: savePlot(vis_prefix,fmt="png")
    plt.show()

def get_options():
    options={}
    options["prefix"]=sys.argv[1]
    if len(sys.argv)==2:
        options["triangle"]=True
        options["x_atom"]=True
        options["ring"]=False
    else:
        if("t" in sys.argv[2]): options["triangle"]=True
        else: options["triangle"]=False
        if("x" in sys.argv[2]): options["x_atom"]=True
        else: options["x_atom"]=False
        if("m" in sys.argv[2]): options["m_atom"]=True
        else: options["m_atom"]=False
        if("r" in sys.argv[2]): options["ring"]=True
        else: options["ring"]=False
        if("A" in sys.argv[2]): options["atom_label"]=True
        else: options["atom_label"]=False
        if("T" in sys.argv[2]): options["triangle_label"]=True
        else: options["triangle_label"]=False
        if("s" in sys.argv[2]): options["save_pdf"]=True
        if("S" in sys.argv[2]): options["save_png"]=True
    return options

def updateParams():
    params = {'legend.fontsize': 8,
              'font.size' : 8,
              'figure.figsize': (6.5, 6.5),
              'axes.labelsize': 8,
              'axes.titlesize': 8,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8}
    pylab.rcParams.update(params)
    return

def get_data(prefix):
    atom_filename = "{0}_atoms.out".format(prefix)
    unit_filename = "{0}_units.out".format(prefix)
    ring_filename = "{0}_rings.out".format(prefix)

    data={}
    atoms=np.genfromtxt(atom_filename,dtype=float)
    elements=np.array(atoms[:,0],dtype=int)
    coordination=np.array(atoms[:,1],dtype=int)
    crds=np.array(atoms[:,2:],dtype=float)
    units=np.genfromtxt(unit_filename,dtype=int)
    file=open(ring_filename,"r")
    rings=[]
    for line in file:
        ring=np.array([int(a) for a in line.split()])
        rings.append(ring)
    file.close()
    ring_sizes=[]
    for r in rings: ring_sizes.append(r.size)
    ring_sizes=np.array(ring_sizes)

    data["elements"]=elements
    data["coordination"]=coordination
    data["crds"]=crds
    data["unit_m"]=units[:,0]
    data["unit_x"]=units[:,1:]
    data["rings"]=rings
    data["ring_sizes"]=ring_sizes

    return data

def plot_triangle_network(data,show_x,show_m,atom_label,tri_label,fig,ax):
    crds=data.get("crds")
    unit_x=data.get("unit_x")
    unit_m=data.get("unit_m")

    # Generate patch drawing commands
    polygon_cmds=generate_polygon_commands(3,3)

    # Loop over x atom triangles and draw
    for tri in unit_x:
        tri_crds=np.array([crds[a] for a in tri])
        tri_crds=np.append(tri_crds, [tri_crds[0]], axis=0)
        path=Path(tri_crds, polygon_cmds[0])
        patch = patches.PathPatch(path, facecolor="black", lw=1.0, alpha=0.5, zorder=1)
        ax.add_patch(patch)

    # Plot x atoms
    if show_x:
        x_crds=crds[unit_x.flatten()]
        ax.scatter(x_crds[:,0],x_crds[:,1], facecolor="red", edgecolor="black", alpha=0.8, s=10, zorder=2)
    # Plot m atoms
    if show_m:
        m_crds=crds[unit_m]
        ax.scatter(m_crds[:,0],m_crds[:,1], facecolor="yellow", edgecolor="black", alpha=0.8, s=10, zorder=2)

    # Make labels
    if atom_label:
        for label, crd in enumerate(crds):
            plt.text(crd[0], crd[1], label, zorder=2)

    if tri_label:
        for label,tri in enumerate(unit_x):
            tri_crd_x=np.average(np.array([crds[a,0] for a in tri]))
            tri_crd_y=np.average(np.array([crds[a,1] for a in tri]))
            plt.text(tri_crd_x, tri_crd_y, label, zorder=2)

    # Set axes limits
    limLb=np.amax([np.amin(crds),np.amax(crds)])*1.1
    limUb=-limLb
    ax.set_xlim(limLb,limUb)
    ax.set_ylim(limLb,limUb)

def plot_ring_network(data,fig,ax):

    # Unpack dictionary
    crds=data.get("crds")
    unit_m=data.get("unit_m")
    rings=data.get("rings")
    ring_sizes=data.get("ring_sizes")

    # Generate patch drawing commands and colours
    polygon_cmds=generate_polygon_commands(3,np.max(ring_sizes));
    ring_colours=generate_colours(ring_sizes)

    # rings=rings[::-1]
    # ring_sizes=ring_sizes[::-1]
    # ring_colours=ring_colours[::-1]
    #
    # Loop over rings and draw polygons
    for i, ring in enumerate(rings):
        ring_crds=np.array([crds[unit_m[a]] for a in ring])
        ring_crds=np.append(ring_crds, [crds[0]], axis=0)
    # depth=np.average(np.array([depth_cueing[a] for a in ring]))
        path=Path(ring_crds, polygon_cmds[ring_sizes[i]-3])
        patch = patches.PathPatch(path, facecolor=ring_colours[i], lw=1.0, alpha=1.0, zorder=0)
        # patch = patches.PathPatch(path, facecolor="white", lw=1.0, alpha=1.0, zorder=0)
        ax.add_patch(patch)

    # if(label):
    #     for i, c in enumerate(crds):
    #         plt.text(c[0],c[1],i)
    #
    # Set axes limits
    limLb=np.amax([np.amin(crds),np.amax(crds)])*1.1
    limUb=-limLb
    ax.set_xlim(limLb,limUb)
    ax.set_ylim(limLb,limUb)
    return

def generate_polygon_commands(min, max):
    polygon_code = [1, 79]
    for i in range(1, min):
        polygon_code.insert(1, 2)
    all_poly_codes = []
    all_poly_codes.append(polygon_code[:])
    for i in range(min, max):
        polygon_code.insert(1, 2)
        all_poly_codes.append(polygon_code[:])
    return all_poly_codes

def generate_colours(ring_sizes):
    n_rings=ring_sizes.size
    colours=[]
    colormap_greens=plt.cm.get_cmap("Greens")
    colormap_blues=plt.cm.get_cmap("Blues")
    colormap_greys=plt.cm.get_cmap("Greys")
    colormap_reds=plt.cm.get_cmap("Reds")
    colormap_oranges=plt.cm.get_cmap("YlOrBr")
    colormap_purples=plt.cm.get_cmap("PuRd")
    colormap_pinks=plt.cm.get_cmap("RdPu")

    green=colormap_greens(100)
    blue=colormap_blues(150)
    grey=colormap_greys(90)
    red=colormap_reds(105)
    orange=colormap_oranges(100)
    purple=colormap_purples(100)
    pink=colormap_pinks(80)

    colourList=[green,blue,grey,red,orange,purple,pink]

    for i in range(n_rings):
        size=int(ring_sizes[i])
        if(size<4 or size>10): colours.append("white")
        else: colours.append(colourList[size-4])

    return colours

# def plot_ring_network(ring_data,label,fig,ax):
#
#     # Unpack dictionary
#     crds = ring_data.get("crds")
#     cnxs = ring_data.get("cnxs")
#     ring_sizes = ring_data.get("ring_sizes")
#
#     for cnxList in cnxs:
#         cnx0 = cnxList[0]
#         for cnx1 in cnxList[1:]:
#             if (cnx0 < cnx1): plt.plot([crds[cnx0][0], crds[cnx1][0]], [crds[cnx0][1], crds[cnx1][1]], color="k",
#                                        lw=0.5, zorder=2)
#     plt.scatter(crds[:, 0], crds[:, 1], c='k', s=2, zorder=3)
#
#     if label:
#         for i, crd in enumerate(crds):
#             plt.text(crd[0],crd[1],i)

def savePlot(prefix,fmt="pdf"):
    filename="{0}.{1}".format(prefix,fmt)
    plt.savefig(filename, dpi=400, bbox_inches="tight")

if __name__=="__main__":
    main()
