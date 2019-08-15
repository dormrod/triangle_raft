import numpy as np
import sys

def main():

    poly_size=int(sys.argv[1])

    # file names
    atom_crd_filename="{0}_atoms.out".format(poly_size)
    unit_filename="{0}_units.out".format(poly_size)
    ring_filename="{0}_rings.out".format(poly_size)
    cnx_filename="{0}_connections.out".format(poly_size)

    # atom coordinates
    file=open(atom_crd_filename,"w")
    # inner x atoms
    rot_ang=2.0*np.pi/poly_size
    l=np.sin(np.pi/3.0)
    r=l/(np.sin(0.5*rot_ang))
    for i in range(poly_size):
        x=r*np.cos(-0.5*rot_ang+i*rot_ang)
        y=r*np.sin(-0.5*rot_ang+i*rot_ang)
        file.write("{:}  {:}  {:8.6f}  {:8.6f}    \n".format(8,4,x,y))
    # outer x atoms
    rot_ang=2.0*np.pi/poly_size
    rr=np.sqrt(r**2-l**2)+2*l*np.sin(np.pi/3)
    for i in range(poly_size):
        x=rr*np.cos(i*rot_ang)
        y=rr*np.sin(i*rot_ang)
        file.write("{:}  {:}  {:8.6f}  {:8.6f}    \n".format(8,3,x,y))
    # si atoms
    rot_ang=2.0*np.pi/poly_size
    rr=np.sqrt(r**2-l**2)+2*l*np.sqrt(3.0)/6.0
    for i in range(poly_size):
        x=rr*np.cos(i*rot_ang)
        y=rr*np.sin(i*rot_ang)
        file.write("{:}  {:}  {:8.6f}  {:8.6f}    \n".format(14,3,x,y))
    file.close()

    # units
    file=open(unit_filename,"w")
    for i in range(poly_size):
        a=i
        b=(i+1)%poly_size
        c=i+poly_size
        d=i+2*poly_size
        file.write("{0} {1} {2} {3}    \n".format(d,a,b,c))
    file.close()

    # rings
    file=open(ring_filename,"w")
    for i in range(poly_size):
        file.write("{0}  ".format(i))
    file.close()

    # connections
    file=open(cnx_filename,"w")
    # # unit-atom m
    # file.write("{}  \n".format(poly_size))
    # for i in range(poly_size):
    #     file.write("{}  {}  \n".format(i,(i+12)))
    # # unit-atom x
    # file.write("{}  \n".format(poly_size*3))
    # for i in range(poly_size):
    #     a=i
    #     b=(i+1)%6
    #     c=i+6
    #     file.write("{}  {}  \n".format(i,a))
    #     file.write("{}  {}  \n".format(i,b))
    #     file.write("{}  {}  \n".format(i,c))
    # unit-unit
    uu=[]
    for i in range(poly_size):
        uu.append(np.array([i,(i+1)%poly_size]))
    uu=np.array(uu)
    file.write("{}  \n".format(uu[:,0].size))
    for p in uu:
        file.write("{}  {}  \n".format(p[0],p[1]))
    # unit-ring
    # file.write("{}  \n".format(poly_size))
    # for i in range(poly_size):
    #     file.write("{}  {}  \n".format(i,0))
    # ring-ring
    file.write("{}".format(0))


    file.close()



if __name__=="__main__":
    main()
