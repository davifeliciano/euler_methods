import sys
from euler_methods import check_dom
import numpy as np
from math import sin, cos, pi, radians, sqrt
import matplotlib.pyplot as plt


def first_quad(deg):
    if abs(deg) == 90:
        return 90
    if deg >= 0:
        return deg % 90
    else:
        return - abs(deg) % 90


def euler_system(funcs, init_values, dom=(0, 1), h=0.1):
    '''
    solves a system of 3 odes in the form y[i]' = funcs[i](x, y[i])
    given initial conditions y[i](dom[0]) = init_values[i] using
    euler explicit method

    funcs: a tuple of functions of the form lambda i, t, x, y, z:
    init_values: a tuple of float initial values for x, y, z
    dom: an ordered 2 float tuple with the t domain to consider
    h: the step of euler explicit method  
    '''
    a, b = check_dom(dom)

    x = np.arange(a, b + h, h)
    ys = [np.array([init_value]) for init_value in init_values]

    for i in range(len(x) - 1):
        ys_i = [y[i] for y in ys]
        for j in range(len(ys)):
            ys[j] = np.append(ys[j], ys[j][i] + h * funcs[j](i, x[i], *ys_i))
    return x, ys


usg_str = f'Usage:\n{sys.argv[0]} [latitude] [velocity] [angle] [direction]\n'
opt_arg_str = f'Optional arguments: [mass] [alpha]'
exit_str = usg_str + opt_arg_str

# check if the cl arguments are valid, if not, show user the way
try:
    lat = radians(first_quad(float(sys.argv[1])))
    velocity = abs(float(sys.argv[2]))
    angle = radians(abs(first_quad(float(sys.argv[3]))))
    direction = float(sys.argv[4])
except ValueError:
    raise SystemExit(exit_str + '\nAll the arguments must be numbers')
except:
    raise SystemExit(exit_str)

try:
    drag = True     # this flag indicates a drag coeficcient was supplied
    alpha = abs(float(sys.argv[6]))
except ValueError:
    raise SystemExit(exit_str + '\nAll the arguments must be numbers')
except IndexError:
    drag = False
    alpha = 0.0
    pass

try:
    mass = abs(float(sys.argv[4]))
except ValueError:
    raise SystemExit(exit_str + '\nAll the arguments must be numbers')
except IndexError:
    mass = 1.0
    pass

arg_info = (
    'Latitude',
    'Speed',
    'Angle of Shot',
    'Direction',
    'Mass',
    'Drag Coefficient'
)

arg_unit = [
    '°',
    ' m/s',
    '°',
    '°',
    ' kg',
    ''
]

if lat >= 0:
    arg_unit[0] += ' N'
else:
    arg_unit[0] += ' S'

for i in range(len(sys.argv) - 1):
    print(arg_info[i] + f' = {float(sys.argv[i + 1])}' + arg_unit[i])

print('\nComputing solution...')

omega = 7.292e-5
g = 9.807


def abs_vec(e, n, z): return sqrt(e ** 2 + n ** 2 + z ** 2)


def e_dot_func(i, t, e, n, z):
    return (- alpha * abs_vec(e, n, z) * e + 2 * omega * (n * sin(lat) - z * cos(lat))) / mass


def n_dot_func(i, t, e, n, z):
    return (- alpha * abs_vec(e, n, z) * n + 2 * omega * e * sin(lat)) / mass


def z_dot_func(i, t, e, n, z):
    return (- alpha * abs_vec(e, n, z) * z + 2 * omega * cos(lat)) / mass - g


e_dot_init = velocity * cos(angle) * sin(direction)
n_dot_init = velocity * cos(angle) * cos(direction)
z_dot_init = velocity * sin(angle)

h = 0.01
dom = (0, 120)

vel_funcs = (e_dot_func, n_dot_func, z_dot_func)
vel_init_values = (e_dot_init, n_dot_init, z_dot_init)

t, vel = euler_system(vel_funcs, vel_init_values, dom, h)

if drag:
    '''
    if alpha is supplied, we're gonna plot two trajectories
    one with drag, another without drag 
    '''
    old_alpha = alpha
    alpha = 0.0
    t_nodrag, vel_nodrag = euler_system(vel_funcs, vel_init_values, dom, h)
    alpha = old_alpha

pos_funcs = (
    lambda i, t, e, n, z: vel[0][i],
    lambda i, t, e, n, z: vel[1][i],
    lambda i, t, e, n, z: vel[2][i]
)

pos_nodrag_funcs = (
    lambda i, t, e, n, z: vel_nodrag[0][i],
    lambda i, t, e, n, z: vel_nodrag[1][i],
    lambda i, t, e, n, z: vel_nodrag[2][i]
)

pos_init_values = (0, 0, 0)

t, pos = euler_system(pos_funcs, pos_init_values, dom, h)

if drag:
    '''
    if alpha is supplied, we're gonna plot two trajectories
    one with drag, another without drag 
    '''
    old_alpha = alpha
    alpha = 0.0
    t_nodrag, pos_nodrag = euler_system(
        pos_nodrag_funcs, pos_init_values, dom, h)
    alpha = old_alpha

# pick only the subdomains where z is positive
for i in range(len(t)):
    if pos[2][i] < 0:
        t = t[:i]
        pos = [p[:i] for p in pos]
        vel = [v[:i] for v in vel]
        break

for i in range(len(t_nodrag)):
    if pos_nodrag[2][i] < 0:
        t_nodrag = t_nodrag[:i]
        pos_nodrag = [p[:i] for p in pos_nodrag]
        vel_nodrag = [v[:i] for v in vel_nodrag]
        break

print('Ploting graphics...')

''' setting up tex in mpl
    you must install texlive-latex-extra texlive-fonts-recommended
    cm-super dvipng ghostscript packages for this '''
plt.rcParams.update({
    'text.usetex': True,
    'figure.figsize': [12.8, 4.8],
    'figure.autolayout': True,
    'axes.labelsize': 15
})

# creating figures and axes
vel_fig, vel_axs = plt.subplots(1, 3)
pos_fig, pos_axs = plt.subplots(1, 2)

fig_3d = plt.figure()
ax1 = fig_3d.add_subplot(1, 3, 1, projection='3d')
ax2 = fig_3d.add_subplot(1, 3, 2, projection='3d')
ax3 = fig_3d.add_subplot(1, 3, 3, projection='3d')
ax_3d = [ax1, ax2, ax3]

# setting up axes and plotting
axes_labels = ('e', 'n', 'z')

ax1.view_init(azim=45, elev=15)
ax2.view_init(azim=-20, elev=20)
ax3.view_init(azim=135, elev=20)

for i in range(3):
    vel_axs[i].grid(ls='--')
    vel_axs[i].set(
        xlabel=r'$t$',
        ylabel=r'$v_{sub}$'.format(sub=axes_labels[i])
    )

    ax_3d[i].set(
        xlabel=r'$e$',
        ylabel=r'$n$'
    )
    ax_3d[i].plot(*pos, label='Com resistência do ar', color='red')
    vel_axs[i].plot(t, vel[i], label='Com resistência do ar', color='red')

    if drag:
        ax_3d[i].plot(*pos_nodrag, label='Sem resistência do ar', color='blue')
        vel_axs[i].plot(t_nodrag, vel_nodrag[i],
                        label='Sem resistência do ar', color='blue')

pos_axs[0].grid(ls='--')
pos_axs[0].set(
    xlabel=r'$e$',
    ylabel=r'$n$',
)

pos_axs[1].grid(ls='--')
pos_axs[1].set(
    xlabel=r'$t$',
    ylabel=r'$z$'
)

pos_axs[0].plot(pos[0], pos[1], label='Com resistência do ar', color='red')
pos_axs[1].plot(t, pos[2], label='Com resistência do ar', color='red')

if drag:
    pos_axs[0].plot(pos_nodrag[0], pos_nodrag[1],
                    label='Sem resistência do ar', color='blue')
    pos_axs[1].plot(t_nodrag, pos_nodrag[2],
                    label='Sem resistência do ar', color='blue')

print('Saving images...')

pos_fig.savefig('images/problem2/position.png', dpi=300)
vel_fig.savefig('images/problem2/velocity.png', dpi=300)
fig_3d.savefig('images/problem2/trajectory.png', dpi=300)
