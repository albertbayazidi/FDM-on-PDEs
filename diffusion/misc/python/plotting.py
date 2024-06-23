#%% 
from utils import *

#%% CONST DIFF NUMERICAL SOL WITH DIRICHLET + INTEGRAL
filename = "../data/numsol_dirichlet.txt"
jump_array = np.array([1,2,4])*100
N = 10000
plot_sol_n_integral(N,jump_array,filename)

#%% CONST DIFF NUMERICAL SOL WITH NEUMANN + INTEGRAL
filename = "../data/numsol_neumann.txt"
jump_array = np.array([1,2,4])*100
N = 10000
plot_sol_n_integral(N,jump_array,filename)

#%% ANALYTICAL UNBOUND
AnalConst  = np.loadtxt("../data/anal_unbound.txt" ,delimiter=" ", dtype=str)
new = np.char.split(AnalConst, sep=",")
dataAnal = lambda x: np.asarray(new[x][:-1], dtype=float)

jump_array = np.array([1,2,4])*100
N = 10000

plot_sol_n_integral(N,jump_array,"../data/anal_unbound.txt")

#%% ANALYTIACL VS NUMERICAL SOL
NumConstDir  = np.loadtxt("../data/numsol_dirichlet.txt" ,delimiter=" ", dtype=str)
newDir = np.char.split(NumConstDir, sep=",")
dataNumDir = lambda x: np.asarray(newDir[x][:-1], dtype=float)

NumConstNeu  = np.loadtxt("../data/numsol_neumann.txt" ,delimiter=" ", dtype=str)
newNeu = np.char.split(NumConstNeu, sep=",")
dataNumNeu = lambda x: np.asarray(newNeu[x][:-1], dtype=float)

jump_array = np.array([1,2,4])*100
data1 = dataNumNeu
data2 = dataAnal
data3 = dataNumDir

plt_anal_vs_num(N,jump_array, data1,data2,data3)
#%% ANALYTICAL BOUND neumann (NO FLUX ON BOUNDRY)
N = 10000-1
jump_array = np.array([1,2,4])*100
filename = "../data/anal_bound_constD_neumann.txt"
plot_sol_n_integral(N,jump_array,filename)

#%% ANAL BOUND VS NUM NEUMANN
AnalBoundNeu  = np.loadtxt("../data/anal_bound_constD_neumann.txt" ,delimiter=" ", dtype=str)
newNeu = np.char.split(AnalBoundNeu, sep=",")
dataAnalNeu = lambda x: np.asarray(newNeu[x][:-1], dtype=float)

jump_array = np.array([1,2,4])*100
data1 = dataNumNeu
data2 = dataAnalNeu

plt_bound_anal_vs_num(N,jump_array, data1, data2, "Neumann BC")
#%% ANALYTICAL BOUND dirichlet (abssoulte absorption on boundry)
N = 10000-1
jump_array = np.array([1,2,4])*100
filename = "../data/anal_bound_constD_dirichlet.txt"
plot_sol_n_integral(N,jump_array,filename)

#%% ANAL BOUND VS NUM DIRICHLET

AnalBoundDir  = np.loadtxt("../data/anal_bound_constD_dirichlet.txt" ,delimiter=" ", dtype=str)
newNeu = np.char.split(AnalBoundDir, sep=",")
dataAnalDir = lambda x: np.asarray(newNeu[x][:-1], dtype=float)

jump_array = np.array([1,2,4])*100
data1 = dataNumDir
data2 = dataAnalDir
plt_bound_anal_vs_num(N,jump_array, data1, data2, "Dirichlet BC")

#%% DIFFUSION DISTRIBUTION  ### NON CONST DIFFUSION ####
plot_diff_dist()

#%% NON CONST DIFF NUMERICAL SOL WITH NEUMANN + INTEGRAL  
N = 10000
jump_array = np.array([1,2,4])*100
filename = "../data/numsolNonConst_neumann.txt"
plot_sol_n_integral(N,jump_array,filename)

#%% NON CONST DIFF NUMERICAL SOL WITH DIRICHLET + INTEGRAL 
N = 10000
jump_array = np.array([1,2,4])*100
filename = "../data/numsolNonConst_dirichlet.txt"
plot_sol_n_integral(N,jump_array,filename)

# %% UNBOUND ANAL SOL WITH NON CONST DIFFUSION
N = 10000
jump_array = np.array([10,15,20])*100
plot_UNBOUND_ANAL_NON_CONST(N,jump_array)# %%

# %% TEST OF THE ANALYTICAL UNBOUND NON CONST DIFFUSION

fileAnal  = np.loadtxt("../data/anal_unbound_nonConstD.txt" ,delimiter=" ", dtype=str)
newfileAnal = np.char.split(fileAnal, sep=",")
dataAnalNonConstDiff = lambda x: np.asarray(newfileAnal[x][:-1], dtype=float)

fileDir  = np.loadtxt("../data/numsolNonConst_dirichlet.txt" ,delimiter=" ", dtype=str)
newfileDir = np.char.split(fileDir, sep=",")
dataNumNonConstDiffDir = lambda x: np.asarray(newfileDir[x][:-1], dtype=float)

fileNeu  = np.loadtxt("../data/numsolNonConst_neumann.txt" ,delimiter=" ", dtype=str)
newfileNeu = np.char.split(fileNeu, sep=",")
dataNumNonConstDiffNeu = lambda x: np.asarray(newfileNeu[x][:-1], dtype=float)


data1 = dataNumNonConstDiffNeu
data2 = dataAnalNonConstDiff
data3 = dataNumNonConstDiffDir
N = 10000

jump_array = np.array([2,4,8])*100

plt_anal_vs_num(N,jump_array, data1,data2,data3)
# %%
#%% making images
N = 10000
plot_diff_dist()

fileNeu  = np.loadtxt("../data/numsolNonConst_neumann.txt" ,delimiter=" ", dtype=str)
newfileNeu = np.char.split(fileNeu, sep=",")
dataNumNonConstDiffNeu = lambda x: np.asarray(newfileNeu[x][:-1], dtype=float)

fileDir  = np.loadtxt("../data/numsolNonConst_dirichlet.txt" ,delimiter=" ", dtype=str)
newfileDir = np.char.split(fileDir, sep=",")
dataNumNonConstDiffDir = lambda x: np.asarray(newfileDir[x][:-1], dtype=float)

jump_array = np.array([1,2,4])*100

plot_sol(dataNumNonConstDiffDir,N,jump_array, N-1)


plot_sol(dataNumNonConstDiffNeu,N,jump_array, N-1)
# %%
