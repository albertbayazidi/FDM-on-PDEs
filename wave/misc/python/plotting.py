#%% 
from utils import *
# %%
K = 20
N = 50

a  = np.loadtxt("../data/numSolSinInit.txt" ,delimiter=" ", dtype=str)
new = np.char.split(a, sep=",")
data = lambda x,K: np.asarray(new[x][:-1], dtype=float).reshape(K,K)

storage = np.zeros((K**2,N))
for i in range(N):
    storage[:,i] = data(i,K).ravel()


# %%  NUMERICAL SOLUTION
Size = (10,10)               

X = np.linspace(0,1,K)
X, Y = np.meshgrid(X, X)
Z = storage

gif_length = 150
png_location = '../animation/pngs/numerical/'
gif_location_name = '../animation/wave_numerical.gif'

surf_animation(X, Y, Z, Size,  png_location,
               gif_location_name, gif_length,skip = 1)
#%%  ANALYTICAL SOLUTION
x = np.linspace(0,1,K)
t = np.linspace(0,1,N)
XX, YY = np.meshgrid(x, x)
Z = np.zeros((K**2,N))

def u(xx,yy,t):
    return np.cos(np.sqrt(5)*t)*np.sin(np.pi*xx)*np.sin(2*np.pi*yy)

for i in range(N):
    Z[:,i] = u(XX,YY,t[i]).ravel()


png_location = '../animation/pngs/analytical/'
gif_location_name = '../animation/wave_analytical.gif'

surf_animation(X, Y, Z, Size,  png_location,
               gif_location_name, gif_length,skip = 1)
# %% CHECKING IF DISSPERIVE OR NAH
from scipy import interpolate, integrate
x = np.linspace(0,1,K)
y = np.linspace(0,1,K)
sol = np.zeros(N)
for i in range(N):

    fun = interpolate.RectBivariateSpline(x,y,storage[:,0].reshape(K,K))
    sol[i] = integrate.dblquad(fun,0,1,0,1)[0]

#%%
plt.plot(sol)
plt.show()
