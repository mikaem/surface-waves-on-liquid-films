import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import PolyCollection
from shenfun import *
from mpl_toolkits.mplot3d import axes3d

N = 24
T = Basis(N, 'F', dtype='D') # Using complex initial data
#Tp = T
Tp = Basis(N, 'F', dtype='D', padding_factor=1.5)
#Tp = Basis(N, 'F', dtype='D', dealias_direct=True)
x = T.mesh()
u = TrialFunction(T)
v = TestFunction(T)
k = T.wavenumbers(scaled=True, eliminate_highest_freq=True)

u_ = Array(T)
Up = Array(Tp)
u_hat = Function(T)
mask = T.get_mask_nyquist()

Re = 5.0
We = 755.
theta = np.pi/2.
h0 = 0.000042
lmbda = 0.0098
alpha = 2*np.pi*h0/lmbda
B_1 = 8./15.*Re-2./3./np.tan(theta)
C_1 = 2./3.*alpha**2*We

def LinearRHS(self, **params):
    L = -2.*inner(Dx(u, 0, 1), v)[0]-alpha*(B_1*inner(Dx(u, 0, 2), v)[0]+C_1*inner(Dx(u, 0, 4), v)[0])
    return L

def NonlinearRHS(self, u, u_hat, rhs, **params):
    rhs.fill(0)
    Up[:] = Tp.backward(u_hat, Up)
    rhs = Tp.forward(-2*Up**2, rhs)
    rhs *= 1j*k
    #for m in range(N//2+1):
    #    for kk in range(N//2+1):
    #        l = m - kk
    #        rhs[m] -= 2j*(kk+l)*u_hat[kk]*u_hat[l]
    return rhs

def ZeroRHS(self, u, u_hat, rhs, **params):
    rhs.fill(0)
    return rhs

# initialize
#u_[:] = h0*0.1*np.exp(1j*x)
#u_hat = T.forward(u_, u_hat)
#u_hat.mask_nyquist(mask)
#u_hat = np.where(abs(u_hat) < 1e-15, 0, u_hat)
u_hat[1] = h0/2.
u_hat0 = u_hat.copy()

data = []
tdata = []
plt.figure()

def update(self, u, u_hat, t, tstep, plot_step, **params):
    if tstep % plot_step == 0 and plot_step > 0:
        u = T.backward(u_hat, u)
        plt.plot(x, u)
        plt.draw()
        plt.pause(1e-6)
        data.append(u_hat.copy())

lindata = []
def linupdate(self, u, u_hat, t, tstep, plot_step, **params):
    if tstep % plot_step == 0 and plot_step > 0:
        lindata.append(u_hat.copy())

dt = 0.002
end_time = 50.
#par = {'plot_step': int(end_time/25/dt)}
par = {'plot_step': 500}

#integrator = RK4(T, L=LinearRHS, N=NonlinearRHS, update=update, **par)
integrator = ETDRK4(T, L=LinearRHS, N=NonlinearRHS, update=update, **par)
linear = ETDRK4(T, L=LinearRHS, N=ZeroRHS, update=linupdate, **par)
integrator.setup(dt)
u_hat = integrator.solve(u_, u_hat, dt, (0, end_time))

linear.setup(dt)
u_hat[:] = u_hat0
u_hat = linear.solve(u_, u_hat, dt, (0, end_time))

# Get solution on denser mesh for nicer plots
Tplot = Basis(N, 'F', dtype='D', padding_factor=3.0)
padded_data = []
for d0 in data:
    dc = Function(Tplot, buffer=d0)
    padded_data.append(dc.backward())

padded_lindata = []
for d0 in lindata:
    dc = Function(Tplot, buffer=d0)
    padded_lindata.append(dc.backward())

xp = Tplot.mesh()

t = end_time
s = []
for d in padded_data:
    s.append(np.vstack((xp, d)).T)

N = len(padded_data)
tdata = np.linspace(0, end_time, N)
ddata = np.array(padded_data)
ldata = np.array(padded_lindata)

fig = plt.figure(figsize=(8, 3))
ax = axes3d.Axes3D(fig)
X, Y = np.meshgrid(xp, tdata)
ax.plot_wireframe(X, Y, ddata, cstride=1000)
ax.set_xlim(0, 2*np.pi)
ax.set_ylim(0, t)
#ax.set_zlim(0, 2000)
ax.view_init(65, -105)
#ax.set_zticks([0, 2000])
ax.grid()


fig2 = plt.figure(figsize=(8,3))
ax2 = fig2.gca(projection='3d')
poly = PolyCollection(s, facecolors=(1, 1, 1, 1), edgecolors='b')
ax2.add_collection3d(poly, zs=tdata, zdir='y')
ax2.set_xlim3d(0, 2*np.pi)
ax2.set_ylim3d(0, t)
ax2.set_zlim3d(0, h0)
ax2.view_init(65, -105)
#ax2.set_zticks([0, 2000])
ax2.grid()

fig3 = plt.figure(figsize=(8, 3))
ax3 = fig3.gca(projection='3d')
X, Y = np.meshgrid(xp, tdata)
ax3.plot_surface(X, Y, ddata, cstride=1, rstride=1, color='w')
ax3.set_xlim(0, 2*np.pi)
ax3.set_ylim(0, t)
ax3.set_zlim(0, h0)
#ax3.set_zlim(0, 2000)
ax3.view_init(65, -105)
#ax3.set_zticks([0, 2000])
ax3.grid()

fig4 = plt.figure(figsize=(8,3))
ax4 = fig4.gca(projection='3d')
for i in range(len(tdata)):
    ax4.plot(xp, ddata[i], tdata[i])
ax4.view_init(65, -105)
#ax4.set_zticks([0, 2000])
ax4.grid()

fig5 = plt.figure(facecolor='k')
ax5 = fig5.add_subplot(111, facecolor='k')
N = len(tdata)
for i in range(N):
    offset = (N-i-1)*h0
    ax5.plot(xp, ddata[N-i-1]+offset, 'w', lw=2, zorder=(i+1)*2)
    ax5.fill_between(xp, ddata[N-i-1]+offset, offset, facecolor='k', lw=0, zorder=(i+1)*2-1)
fig5.savefig('film.png')
#plt.show()
