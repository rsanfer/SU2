from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.tri as mtri
import numpy as np


def plotmodes(X_mode, Y_mode, Z_mode,i):
    matplotlib.interactive(True)
    fig = plt.figure(i)   
    ax = fig.gca(projection='3d')
    #ax.set_aspect('equal')
    
    # Create cubic bounding box to simulate equal aspect ratio
    max_range = np.array([X_mode.max()-X_mode.min(), Y_mode.max()-Y_mode.min(), Z_mode.max()-Z_mode.min()]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X_mode.max()-X_mode.min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y_mode.max()-Y_mode.min())
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z_mode.max()-Z_mode.min())
    # Comment or uncomment following both lines to test the fake bounding box:
    for xb, yb, zb in zip(Xb, Yb, Zb):
        ax.plot([xb], [yb], [zb], 'w')
    
   

    #triang = mtri.Triangulation(Aero_Matrix[:,0], Aero_Matrix[:,1], Connectivity.tolist()) 
    ax.scatter(X_mode, Y_mode,Z_mode,  linewidth=0.2, color ='b') #triang,
    #ax.scatter(X_mode2, Y_mode2,Z_mode2,  linewidth=0.2, color ='r', marker='^') #triang,

    plt.grid()
    plt.pause(.1)
    plt.draw()   
    #print("something")
    #wait = input("PRESS ENTER TO CONTINUE.")

    #print("something")
    
    