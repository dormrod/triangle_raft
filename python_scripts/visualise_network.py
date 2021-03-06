import sys
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib._color_data as mcd
import matplotlib.colors as colors
import matplotlib.pylab as pylab
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize
import numpy as np
import os

def main():

    # Get options for visualisation and unpack
    vis_options=get_options()
    vis_prefix=vis_options.get("prefix")
    vis_tri=vis_options.get("triangle",False)
    vis_x=vis_options.get("x_atom",False)
    vis_m=vis_options.get("m_atom",False)
    vis_ring=vis_options.get("ring",False)
    vis_ring_colour=vis_options.get("ring_colour_style",0)
    vis_ring_filter=vis_options.get("ring_size_filter",None)
    vis_atom_label=vis_options.get("atom_label",False)
    vis_tri_label=vis_options.get("triangle_label",False)
    vis_boundary=vis_options.get("boundary",False)
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
        plot_ring_network(data,fig,ax,vis_ring_colour,vis_ring_filter)
    if vis_boundary:
        plot_boundary(data,fig,ax)

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
        if("r" in sys.argv[2]):
            options["ring"]=True
            options["ring_colour_style"]=3
        else: options["ring"]=False
        if("A" in sys.argv[2]): options["atom_label"]=True
        else: options["atom_label"]=False
        if("T" in sys.argv[2]): options["triangle_label"]=True
        else: options["triangle_label"]=False
        if("b" in sys.argv[2]): options["boundary"]=True
        else: options["boundary"]=False
        if("e" in sys.argv[2]): options["ring_colour_style"]=1
        if("c" in sys.argv[2]): options["ring_colour_style"]=2
        if("s" in sys.argv[2]): options["save_pdf"]=True
        if("S" in sys.argv[2]): options["save_png"]=True
    if len(sys.argv)==4:
        options["ring_size_filter"]=int(sys.argv[3][1:])
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
    boundary_filename = "{0}_boundary.out".format(prefix)
    colour_filename="{0}_vis.out".format(prefix)

    data={}
    atoms=np.genfromtxt(atom_filename,dtype=float)
    elements=np.array(atoms[:,0],dtype=int)
    coordination=np.array(atoms[:,1],dtype=int)
    crds=np.array(atoms[:,2:],dtype=float)
    units=np.genfromtxt(unit_filename,dtype=int)
    dimensionality=crds.shape[1]
    file=open(ring_filename,"r")
    rings=[]
    for line in file:
        ring=np.array([int(a) for a in line.split()])
        rings.append(ring)
    file.close()
    ring_sizes=[]
    for r in rings: ring_sizes.append(r.size)
    ring_sizes=np.array(ring_sizes)
    if os.path.isfile(boundary_filename):
        boundary=np.genfromtxt(boundary_filename)
    else: boundary=[]
    colours=np.genfromtxt(colour_filename,dtype="int")
    try:
        ring_edges=colours[:,1]
        ring_clusters=colours[:,2]
    except:
        ring_edges=[]
        ring_clusters=[]

    data["elements"]=elements
    data["coordination"]=coordination
    data["crds"]=crds
    data["unit_m"]=units[:,0]
    data["unit_x"]=units[:,1:]
    data["rings"]=rings
    data["ring_sizes"]=ring_sizes
    data["boundary"]=boundary
    data["ring_edges"]=ring_edges
    data["ring_clusters"]=ring_clusters
    data["dimensionality"]=dimensionality

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
        #ax.scatter(x_crds[:,0],x_crds[:,1], facecolor="red", edgecolor="black", alpha=1.0, s=50, zorder=2)
        ax.scatter(x_crds[:,0],x_crds[:,1], facecolor="red", edgecolor="red", alpha=1, s=10, zorder=2)
    # Plot m atoms
    if show_m:
        m_crds=crds[unit_m]
        #ax.scatter(m_crds[:,0],m_crds[:,1], facecolor="yellow", edgecolor="black", alpha=1.0, s=50, zorder=2)
        ax.scatter(m_crds[:,0],m_crds[:,1], facecolor="yellow", edgecolor="yellow", alpha=1, s=10, zorder=2)

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
    #limLb=6.4095262
    #limLb=10
    limUb=-limLb
    ax.set_xlim(limLb,limUb)
    ax.set_ylim(limLb,limUb)

def plot_ring_network(data,fig,ax,colour_style,size_filter):

    # Unpack dictionary
    crds=data.get("crds")
    unit_m=data.get("unit_m")
    rings=data.get("rings")
    ring_sizes=data.get("ring_sizes")
    ring_edges=data.get("ring_edges")
    ring_clusters=data.get("ring_clusters")
    dimensionality=data.get("dimensionality")

    # Generate patch drawing commands and colours
    polygon_cmds=generate_polygon_commands(3,np.max(ring_sizes));
    if colour_style==0: ring_colours=generate_colours(ring_sizes,colour_style)
    elif colour_style==1: ring_colours=generate_colours(ring_edges,colour_style)
    elif colour_style==2:
        if size_filter is not None:
            for i,s in enumerate(ring_sizes):
                if s!=size_filter: ring_clusters[i]=-1
        ring_colours=generate_colours(ring_clusters,colour_style)
    elif colour_style==3: ring_colours=generate_colours(ring_sizes,colour_style)

    # rings=rings[::-1]
    # ring_sizes=ring_sizes[::-1]
    # ring_colours=ring_colours[::-1]
    #

    if dimensionality==2:
        # Loop over rings and draw polygons
        for i, ring in enumerate(rings):
            ring_crds=np.array([crds[unit_m[a]] for a in ring])
            ring_crds=np.append(ring_crds, [crds[0]], axis=0)
        # depth=np.average(np.array([depth_cueing[a] for a in ring]))
            path=Path(ring_crds, polygon_cmds[ring_sizes[i]-3])
            patch = patches.PathPatch(path, facecolor=ring_colours[i], lw=0.6, alpha=1.0, zorder=0)
            # patch = patches.PathPatch(path, facecolor="white", lw=1.0, alpha=1.0, zorder=0)
            ax.add_patch(patch)
        #    plt.text(ring_crds[0,0],ring_crds[0,1],i)
    elif dimensionality==3:
        # Add depth cueing
        z=crds[:,2]
        crds=crds[:,:2]
        z=z/np.max(z)
        for i,a in enumerate(z):
            if a<0: z[i]=0
        for i, ring in enumerate(rings):
            ring_crds=np.array([crds[unit_m[a]] for a in ring])
            ring_crds=np.append(ring_crds, [crds[0]], axis=0)
            depth=np.average(np.array([z[a] for a in ring]))
            path=Path(ring_crds, polygon_cmds[ring_sizes[i]-3])
            patch = patches.PathPatch(path, facecolor=ring_colours[i], lw=1.0, alpha=depth, zorder=0)
            ax.add_patch(patch)
        pass

     #if(label):
     #    for i, c in enumerate(crds):
     #        plt.text(c[0],c[1],i)
    #
    # Set axes limits
    limLb=np.amax([np.amin(crds),np.amax(crds)])*1.1
    limLb=-15
    limUb=-limLb
    ax.set_xlim(limLb,limUb)
    ax.set_ylim(limLb,limUb)
    ax.set_xlim(limLb,limUb)
    return

def plot_boundary(data,fix,ax):
    #unpack
    boundary=data.get("boundary")
    ax.plot(boundary[:,0],boundary[:,1],lw=2)
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

def generate_colours(ring_codes, colour_set=0):
    n_rings=ring_codes.size
    colours=[]
    colormap_greens=plt.cm.get_cmap("Greens")
    colormap_blues=plt.cm.get_cmap("Blues")
    colormap_greys=plt.cm.get_cmap("Greys")
    colormap_reds=plt.cm.get_cmap("Reds")
    colormap_oranges=plt.cm.get_cmap("YlOrBr")
    colormap_purples=plt.cm.get_cmap("PuRd")
    colormap_pinks=plt.cm.get_cmap("RdPu")
    colormap_rainbow=plt.cm.get_cmap("gist_rainbow")
    colormaps=[colormap_greens,colormap_blues,colormap_reds,colormap_oranges,colormap_purples,colormap_pinks]

    green=colormap_greens(100)
    blue=colormap_blues(150)
    grey=colormap_greys(90)
    red=colormap_reds(105)
    orange=colormap_oranges(100)
    purple=colormap_purples(100)
    pink=colormap_pinks(80)
        
    map_lower = cm.get_cmap('Blues_r', 128)
    map_upper = cm.get_cmap('Reds', 128)
    map_mean=cm.get_cmap("Greys")
    map_lower=ListedColormap(map_lower(np.arange(30,100)))
    map_upper=ListedColormap(map_upper(np.arange(30,100)))

    norm_lower=Normalize(vmin=3,vmax=6)
    norm_upper=Normalize(vmin=6,vmax=12)
    colour_mean=map_mean(50)
     #colour_mean='ivory'
     # map = np.vstack((map_upper(np.linspace(0, 1, 128)),map_lower(np.linspace(0, 1, 128))))
     # colormap = ListedColormap(map, name='custom')


    if colour_set==0:
        colourList=[green,blue,grey,red,orange,purple,pink]
        for i in range(n_rings):
            size=int(ring_codes[i])
            if(size<4 or size>10): colours.append("white")
            else: colours.append(colourList[size-4])
    elif colour_set==1:
        colourList=["red","blue","green","orange",grey]
        for code in ring_codes:
            colours.append(colourList[code])
    elif colour_set==2:
        colourList=[]
        np.random.seed(seed=0)
        for i in range(np.max(ring_codes)+1):
            lb=(i%6)*40
            ub=(i%6+1)*40
            rand=np.random.randint(low=lb,high=ub)
            colourList.append(colormap_rainbow(rand))
        #     rand1=np.random.randint(low=0,high=6)
        #     rand2=np.random.randint(low=50,high=200)
        #     colourList.append(colormaps[rand1](rand2))
        colourList.append(grey)
        for code in ring_codes:
            colours.append(colourList[code])
    elif colour_set==3:  
        colourList=[]
        for i in range(30):
            if i < 3:
                colourList.append("white")
            elif np.abs(i-6)<1e-6:
                colourList.append(colour_mean)
            elif i<6:
                colourList.append(map_lower(norm_lower(i)))
            else:
                colourList.append(map_upper(norm_upper(i)))
        for code in ring_codes:
            colours.append(colourList[code])

    return colours


def savePlot(prefix,fmt="pdf"):
    filename="{0}.{1}".format(prefix,fmt)
    plt.savefig(filename, dpi=800, bbox_inches="tight")

if __name__=="__main__":
    main()
