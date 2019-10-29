from shenfun import *
import matplotlib.pyplot as plt

N = 24
T = Basis(N, 'F', dtype='D') # Using complex initial data
Tp = Basis(N, 'F', dtype='D', padding_factor=1.5)
x = T.mesh()
u = TrialFunction(T)
v = TestFunction(T)
k = T.wavenumbers(scaled=True, eliminate_highest_freq=True)

# Choose papameters
Re = 5.0
We = 755.
theta = np.pi/2.
h0 = 0.000042
lmbda = 0.0098
alpha = 2*np.pi*h0/lmbda
B_1 = 8./15.*Re-2./3./np.tan(theta)
C_1 = 2./3.*alpha**2*We

# Arrays to hold the solution
eta_ = Array(T)
etap = Array(Tp)      # For the nonlinear part
eta_hat = Function(T)
eta_lin_hat = Function(T)
eta_lin = Array(T)

# Initialize wavenumber w0
w0 = 1
eta_hat[w0] = h0

# The integrator needs one function to compute the linear part,
# and one term for the nonlinear.
def LinearRHS(self, **params):
    """Return linear operator for right hand side."""
    L = -2*1j*k - alpha*(-B_1*k**2 + C_1*k**4)
    return L

def NonlinearRHS(self, eta, eta_hat, rhs, **params):
    """Return nonlinear part of right hand side."""
    rhs.fill(0)
    etap[:] = Tp.backward(eta_hat, etap)
    rhs = Tp.forward(-2*etap**2, rhs)
    rhs *= 1j*k
    return rhs

# Linear solution at time t
def get_linear_solution(u_hat, t):
    u_hat[w0] = h0*np.exp(LinearRHS(None)[w0]*t)
    return u_hat.backward()

data = []
def update(self, eta, eta_hat, t, tstep, plot_step, fig, **par):
    """Function called at the end of each time step."""
    if tstep % plot_step == 0 and plot_step > 0:
        eta = eta_hat.backward(eta)
        eta_lin = get_linear_solution(eta_lin_hat, t)
        fig.add_subplot(121)
        fig.gca().plot(x, eta)
        fig.add_subplot(122)
        fig.gca().plot(x, eta-eta_lin)
        plt.pause(1e-6)
        data.append(eta_hat.copy())

# Specify time step size and end time.
dt = 0.002
end_time = 2.
fig = plt.figure(1, figsize=(10, 3))
fig.add_subplot(121)
plt.title('$\eta$')
fig.add_subplot(122)
plt.title('nonlinear contribution')
par = {'plot_step': 100, 'fig': fig}
integrator = ETDRK4(T, L=LinearRHS, N=NonlinearRHS, update=update, **par)
integrator.setup(dt)
eta_hat = integrator.solve(eta_, eta_hat, dt, (0, end_time))
plt.show()