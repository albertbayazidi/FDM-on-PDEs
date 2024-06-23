import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

def plot_sol(data,N,jumps ,T = 50):
    data0 = data(0) 
    size = len(data0)
    array = np.linspace(0,1,size)
    colors = ["tab:orange","tab:green","tab:red","tab:purple","tab:brown","tab:pink","tab:gray","tab:olive","tab:cyan"]

    for i,j in enumerate(jumps):
        plt.plot(array,data(j) ,label=f"{np.round(j/N,4)} sek",color = f"{colors[(i+8)%len(colors)]}")


    plt.ylabel(r"Density")
    plt.xlabel("Space ")
    plt.grid()
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=True, ncol=8)
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


def plot_animation(X, Y, N , Size,  png_location,
                   gif_location_name, gif_length = 200,skip = 2):
    
    # The amount of frames in the animation
    frames = N

    maxX = np.max(X) * 1.1
    minX = np.min(X) - 0.1
    maxY = np.max(Y(1)) * 1.1 
    minY = np.min(Y(1)) - 0.1

    # Generate each frame
    for n in range(0,frames,skip):
        fig = plt.figure(figsize=Size)
        ax = fig.add_subplot()

        ax.plot(X,Y(n))

        # Set labels
        ax.set_xlim(X[0],X[-1])
        ax.set_ylim(minY,maxY)

        plt.savefig(f'{png_location + str(n)}.png')
        plt.close()

    images = [Image.open(f'{png_location + str(n)}.png') for n in range(0,frames,skip)]

    images[0].save(f'{gif_location_name}', save_all=True,
                    append_images=images[1:], duration=gif_length, loop=0)