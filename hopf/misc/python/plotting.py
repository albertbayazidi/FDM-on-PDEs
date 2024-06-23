#%% 
from utils import *
# %% ## const flux sols
sol  = np.loadtxt("../data/numsolConstFluxHat.txt" ,delimiter=" ", dtype=str)
new_sol = np.char.split(sol, sep=",")
dataConst = lambda x: np.asarray(new_sol[x][:-1], dtype=float)

# %% const flux sols IMAGES
jumps = np.array([0,2,4])*10
N = 100

plot_sol(dataConst,N,jumps)

# %% const flux sols GIF
N = 100
X = np.linspace(0,1,len(dataConst(0)))
Y = dataConst
Size = (10,10)
png_location = '../animation/pngs/'
gif_location_name = '../animation/ConstFlux.gif'

plot_animation(X, Y, N ,Size, png_location, gif_location_name)

# %% ## NONCONST FLUX SOLS
sol  = np.loadtxt("../data/numsolNonConstFluxWaveFront.txt" ,delimiter=" ", dtype=str)
new_sol = np.char.split(sol, sep=",")
dataNonConst = lambda x: np.asarray(new_sol[x][:-1], dtype=float)

x = np.linspace(0,1,len(dataNonConst(0)))
plt.plot(x ,dataNonConst(0))

# %%  ## NONCONST FLUX GIF
N = 400
X = np.linspace(0,1,len(dataNonConst(0)))
Y = dataNonConst
Size = (10,10)
png_location = '../animation/pngs/'
gif_location_name = '../animation/NonConstFluxWideWave.gif'

plot_animation(X, Y, N , Size, png_location,
                gif_location_name, skip=4)

# %%
