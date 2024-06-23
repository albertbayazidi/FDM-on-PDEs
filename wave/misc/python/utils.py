import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

def plot_sol(data,N,jump = 10,T = 50):
    for i in range(T):
        if i % jump == 0:
            plt.plot(np.linspace(0,1,len(data(i))),data(i), label=f"{np.round(i/N,4)} sek")

    plt.ylabel(r"Density [$\frac{kg}{m^3}$]")
    plt.xlabel("Space [m]")
    plt.grid()
    plt.legend()
    plt.show()

def plot(obj1,obj2,title_x,title_y):
    plt.plot(obj1,obj2)
    plt.xlabel(title_x)
    plt.ylabel(title_y)
    plt.grid()
    plt.show()


def compute_integral(data,T = 50):
    integral = np.zeros(T)
    for i in range(T):
        integral[i] = np.trapz(data(i))

    return integral



def surf_animation(X, Y, Z, Size,  png_location,
                   gif_location_name, gif_length = 200,decimals = 3, color_map = 'plasma',skip = 1):
    
    xDim = X.shape[1]
    yDim = Y.shape[0]

    # The amount of frames in the animation
    frames = len(Z[0,:])

    # Max min solution
    maxZ = np.round(np.max(Z),decimals)
    minZ = np.round(np.min(Z),decimals)

    # Generate each frame
    for n in range(0,frames,skip):
        fig = plt.figure(figsize=Size)
        ax = fig.add_subplot(111, projection='3d')

        ax.plot_surface(X, Y, Z[:,n].reshape(xDim,yDim), cmap=color_map)
        
        # Set ticks
        ax.set_xticks([X[0,0],X[-1,-1]])
        ax.set_yticks([Y[0,0],Y[-1,-1]])
        ax.set_zticks([minZ,maxZ])
        
        # Set labels
        ax.set_xlim(X[0,0],X[-1,-1])
        ax.set_ylim(Y[0,0],Y[-1,-1])
        ax.set_zlim(minZ,maxZ)

        plt.savefig(f'{png_location + str(n)}.png')
        plt.close()

    images = [Image.open(f'{png_location + str(n)}.png') for n in range(0,frames,skip)]

    images[0].save(f'{gif_location_name}', save_all=True,
                    append_images=images[1:], duration=gif_length, loop=0)