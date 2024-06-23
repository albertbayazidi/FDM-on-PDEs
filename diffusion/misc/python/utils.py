import numpy as np
import matplotlib.pyplot as plt
from scipy import special

def plot_sol(data,N,jumps ,T = 50):
    data0 = data(0) 
    size = len(data0)
    array = np.linspace(0,1,size)
    colors = ["tab:orange","tab:green","tab:red","tab:purple","tab:brown","tab:pink","tab:gray","tab:olive","tab:cyan"]
    plt.plot([0.5,0.5],[0,np.max(data0)],"-.", label=f"{np.round(0/N,4)} sek", color = "tab:blue")
    plt.plot([0,1],[np.min(data0),np.min(data0)],"-.",color = "tab:blue")
    
    for i,j in enumerate(jumps):
        plt.plot(array,data(j) ,label=f"{np.round(j/N,4)} sek",color = f"{colors[(i+8)%len(colors)]}")

    plt.plot(array,data(N-1), label=f"1 sek",color = f"{colors[(i+9)%len(colors)]}")

    
    plt.ylabel(r"Density [$\frac{kg}{m}$]")
    plt.xlabel("Space [m]")
    plt.ylim(-0.5, 1.1*np.max(data(jumps[0])))
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


def plot_sol_n_integral(N,jump_array,filename):
    # exstract data
    numConst  = np.loadtxt(filename ,delimiter=" ", dtype=str)
    splitNumConst = np.char.split(numConst, sep=",")
    data = lambda x: np.asarray(splitNumConst[x][:-1], dtype=float)

    #PLOT CONSTANTS DIFFUSION NUMERICAL SOLUTION
    
    plot_sol(data,N,jump_array, T = N)

    #PLOT INTEGRAL 
    integral = compute_integral(data,N)
    array = np.linspace(0,1,len(integral))
    plot(array,integral,"time","mass")

    print("Standard deviation of integral: ", np.std(integral))


def plot_diff_dist(filename= "../data/diffDist.txt"):
    # GET DATA
    diffDist  = np.loadtxt(filename ,delimiter=" ", dtype=str)
    new = np.char.split(diffDist, sep=",")
    diff = lambda x: np.asarray(new[x][:-1], dtype=float)

    # PLOT DIFFUSION DISTRIBUTION
    array = np.linspace(0,1,len(diff(0)))
    plot(array,diff(0),"x","Distribution")


def plot_UNBOUND_ANAL_NON_CONST(N,jump_array):
    file  = np.loadtxt("../data/anal_unbound_nonConstD.txt" ,delimiter=" ", dtype=str)
    newfile = np.char.split(file, sep=",")
    data = lambda x: np.asarray(newfile[x][10:-1], dtype=float)
    data0 = data(0) 
    size = len(data0)
    array = np.linspace(0,1,size)
    colors = ["tab:orange","tab:green","tab:red","tab:purple","tab:brown","tab:pink","tab:gray","tab:olive","tab:cyan"]
    plt.plot([0.5,0.5],[0,np.max(data0)],"-.", label=f"{np.round(0/N,4)} sek", color = "tab:blue")
    plt.plot([0,1],[np.min(data0),np.min(data0)],"-.",color = "tab:blue")
    
    for i,j in enumerate(jump_array):
        plt.plot(array,data(j) ,label=f"{np.round(j/N,4)} sek",color = f"{colors[(i+8)%len(colors)]}")

    plt.plot(array,data(N-1), label=f"1 sek",color = f"{colors[(i+9)%len(colors)]}")

    
    plt.ylabel(r"Density [$\frac{kg}{m}$]")
    plt.xlabel("Space [m]")
    plt.ylim(-0.5, 1.1*np.max(data(jump_array[0])))
    plt.grid()
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=True, ncol=8)
    plt.show()

def plt_anal_vs_num(N,jump_array, data1,data2,data3):
    x = np.linspace(0,1,len(data1(0)))
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=1, ncols=3,
                                        sharex=True,figsize=(18, 6))
    for j,i in enumerate(jump_array):
        ax0.set_title('Numerical solution with Neumann BC')
        ax0.plot(x, data1(i), label=f'{round(i/N,4)} sek')
        ax0.grid()
        ax0.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=True, ncol=9)

        ax1.set_title('Analytical solution with no BC')
        ax1.plot(x, data2(i), label=f'{round(i/N,4)} sek')
        ax1.grid()
        ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=True, ncol=9)

        ax2.set_title('Numerical solution with Dirichlet BC')
        ax2.plot(x, data3(i), label=f'{round(i/N,4)} sek')
        ax2.grid()
        ax2.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=True, ncol=9)

    fig.tight_layout(pad=5.0)
    fig.suptitle('Numerical vs analytical solutions')
    plt.show()


def plt_bound_anal_vs_num(N,jump_array, data1,data2,bc_condition = "no BC"):
    x = np.linspace(0,1,len(data1(0)))
    fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2,
                                        sharex=True,figsize=(18, 6))
    for _,i in enumerate(jump_array):
        ax0.set_title(f'Numerical solution with {bc_condition}')
        ax0.plot(x, data1(i), label=f'{round(i/N,4)} sek')
        ax0.grid()
        ax0.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=True, ncol=9)

        ax1.set_title(f'Analytical solution with {bc_condition}')
        ax1.plot(x, data2(i), label=f'{round(i/N,4)} sek')
        ax1.grid()
        ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=True, ncol=9)


    fig.tight_layout(pad=5.0)
    fig.suptitle('Numerical vs analytical solutions')
    plt.show()


    