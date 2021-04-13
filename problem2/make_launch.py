import argparse
import numpy as np
from math import sin, cos, pi, sqrt, radians, degrees, modf
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def first_quad(deg):
    if abs(deg) <= 90:
        return deg
    else:
        if deg > 0:
            return deg % 90
        else:
            return - abs(deg) % 90


def check_dom(dom):
    '''
    Check if dom is a valid domain
    i. e. a 2 elem tuple with distinct number values
    Return the values in crescent order
    '''
    if not isinstance(dom, tuple) or len(dom) != 2:
        raise TypeError('dom must be a 2-element tuple')

    try:
        dom0 = float(dom[0])
        dom1 = float(dom[1])
    except TypeError:
        raise TypeError('The elems of dom must be numbers')

    if dom0 == dom1:
        raise ValueError('The elems of dom must be distinct')

    if dom0 < dom1:
        return dom0, dom1
    else:
        return dom1, dom0


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


''' 
not the better way to solve this problem of a flying projectile
since in this case we're often interest in the time interval where
the solution for z is positive, i. e. the projectile is still flying

with this algorithm, we do not ever check if z < 0, and, depending on
the supplied arguments, you may want to rerun the algorithm with a
different dom.

but it was intended in order to make the method general
'''

# using argparse to parse command line arguments

parser = argparse.ArgumentParser(
    description='plot the trajectory of a projectile in the atmosphere, considering Coriolis effect and air resistance',
    epilog='the images of the generated graphics for positions and velocities are saved in current dir'
)

parser.add_argument(
    'velocity', type=float, help='the initial velocity of the shot, in m/s')
parser.add_argument('angle', type=float,
                    help='the azimuthal angle of the shot')
parser.add_argument(
    'direction', type=float, help='the direction of the shot, as in a compass')
parser.add_argument('-l', '--latitude', type=float, nargs='?',
                    default=40.0, const=40.0, help='the latitude from where the shot is taken')
parser.add_argument('-m', '--mass', type=float, nargs='?',
                    default=1.0, const=1.0, help='the mass of the projectile')
parser.add_argument('-d', '--drag', type=float, nargs='?',
                    default=0, const=1.0, help='the drag coefficient')
parser.add_argument('-t', '--time', type=float, nargs='?',
                    default=120.0, const=120.0, help='time domain of the solution')
parser.add_argument('-s', '--scale', type=float, nargs='?',
                    default=100, const=100, help='0 - 100 scale between the plots with and without air resistance')

args = parser.parse_args()

velocity = abs(args.velocity)
angle = radians(abs(first_quad(args.angle)))
direction = radians(args.direction)
lat = radians(first_quad(args.latitude))
mass = abs(args.mass)
alpha = abs(args.drag)
t_max = abs(args.time)

if args.scale == 100:
    scale = 1.0
else:
    frac, whole = modf((args.scale % 100) / 100)
    scale = abs(frac)

# this strings will be the labels of some plots if drag coeff was supplied
plot_label = 'With Air Resistance'
plot_label_nodrag = 'Without Air Resistance'
drag = True     # this flag indicates that a drag coefficient was supplied
if not alpha:
    plot_label = plot_label_nodrag
    drag = False

# print the parameters of the shot
arg_info = (
    'Speed',
    'Angle of Shot',
    'Direction',
    'Latitude',
    'Mass',
    'Drag Coefficient'
)

arg = [
    velocity,
    abs(first_quad(args.angle)),
    args.direction,
    abs(degrees(lat)),
    mass,
    alpha
]

arg_unit = [
    ' m/s',
    '°',
    '°',
    '°',
    ' Kg',
    ''
]

if lat >= 0:
    arg_unit[3] += ' N'
else:
    arg_unit[3] += ' S'

for i in range(len(arg)):
    print(arg_info[i] + f' = {arg[i]}' + arg_unit[i])

print('\nComputing solution...')

omega = 7.292e-5    # angular velocity of the Earth, in Hz
g = 9.807           # Earth's gravitational acceleration near the surface


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
dom = (0, t_max)

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

if drag:
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
    'axes.labelsize': 15,
    'axes.formatter.limits': (-1, 4)
})

# creating figures and axes
vel_fig, vel_axs = plt.subplots(1, 3)
pos_fig, pos_axs = plt.subplots(1, 2)

fig_3d = plt.figure()

ax1 = fig_3d.add_subplot(1, 3, 1, projection='3d')
ax2 = fig_3d.add_subplot(1, 3, 2, projection='3d')
ax3 = fig_3d.add_subplot(1, 3, 3, projection='3d')
ax_3d = [ax1, ax2, ax3]

fig_iter, ax_iter = plt.subplots(subplot_kw={"projection": "3d"})
fig_iter.set_size_inches(6.4, 4.8)

# using tight_layout to prevent overlap in subplots
pos_fig.tight_layout(pad=4.0, w_pad=2.5)
vel_fig.tight_layout(pad=3.0, w_pad=2.5)
fig_3d.tight_layout(pad=4.0, w_pad=1.5)

# setting initial view for the 3d plots of fig_3d
ax1.view_init(azim=45, elev=15)
ax2.view_init(azim=-20, elev=20)
ax3.view_init(azim=135, elev=20)

# setting axes labels and grids and ploting
axes_labels = ('e', 'n', 'z')
max_index = int(scale * (len(t_nodrag) - len(t) + len(t)))

for i in range(3):
    vel_axs[i].grid(ls='--')
    vel_axs[i].set(
        xlabel=r'$t$',
        ylabel=r'$\dot {sub}$'.format(sub=axes_labels[i])
    )

    ax_3d[i].set(
        xlabel=r'$e$',
        ylabel=r'$n$'
    )

    ax_3d[i].plot(*pos, label=plot_label, color='red')
    vel_axs[i].plot(t, vel[i], label=plot_label, color='red')

    if drag:
        ax_3d[i].plot(
            pos_nodrag[0][:max_index],
            pos_nodrag[1][:max_index],
            pos_nodrag[2][:max_index],
            label=plot_label_nodrag,
            color='blue'
        )
        vel_axs[i].plot(t_nodrag, vel_nodrag[i],
                        label=plot_label_nodrag, color='blue')

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

ax_iter.set(
    xlabel=r'$e$',
    ylabel=r'$n$',
    zlabel=r'$z$'
)

pos_axs[0].plot(pos[0], pos[1], label=plot_label, color='red')
pos_axs[1].plot(t, pos[2], label=plot_label, color='red')
ax_iter.plot(*pos, label=plot_label, color='red')

if drag:
    pos_axs[0].plot(pos_nodrag[0][:max_index], pos_nodrag[1][:max_index],
                    label=plot_label_nodrag, color='blue')
    pos_axs[1].plot(t_nodrag, pos_nodrag[2],
                    label=plot_label_nodrag, color='blue')
    ax_iter.plot(
        pos_nodrag[0][:max_index],
        pos_nodrag[1][:max_index],
        pos_nodrag[2][:max_index],
        label=plot_label_nodrag,
        color='blue'
    )

# setting up legends and texts
text_str_1 = '\n'.join((
    r'$v_0 = %.2f$' % (velocity, ) + ' m/s',
    r'$\theta = %.2f$' % (abs(first_quad(args.angle)), ) + '°',
    r'$\phi = %.2f$' % (args.direction, ) + '°'
))

text_str_2 = '\n'.join((
    r'$\lambda = %.2f$' % (first_quad(args.latitude), ) + '°',
    r'$m = %.2f$' % (mass, ) + ' Kg',
    r'$\alpha = %.1e$' % (alpha, )
))

ax_iter.text2D(-0.2, 0.95, text_str_1,
               transform=ax_iter.transAxes, fontsize=14)
ax_iter.text2D(0.90, 0.95, text_str_2,
               transform=ax_iter.transAxes, fontsize=14)

ax1.text2D(0.05, 0.90, text_str_1, transform=ax1.transAxes, fontsize=12)
ax1.text2D(0.75, 0.90, text_str_2, transform=ax1.transAxes, fontsize=12)


ax_iter.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1))
ax2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05))

plt.close(pos_fig)
plt.close(vel_fig)
plt.close(fig_3d)

print('Saving images...')

pos_fig.savefig('position.png', dpi=300)
vel_fig.savefig('velocity.png', dpi=300)
fig_3d.savefig('trajectory.png', dpi=300)

print('Done!')

shot_range = sqrt(pos[0][-1] ** 2 + pos[1][-1] ** 2)
max_height = pos[2].max()
flight_time = t[-1]

if drag:
    shot_range_nodrag = sqrt(pos_nodrag[0][-1] ** 2 + pos_nodrag[1][-1] ** 2)
    terminal_velocity = abs(vel[2][-1])
    max_height_nodrag = pos_nodrag[2].max()
    flight_time_nodrag = t_nodrag[-1]

print('\nShot Distance:')
print(plot_label + f': {shot_range:.2f} m')

if drag:
    print(plot_label_nodrag + f': {shot_range_nodrag:.2f} m')
    print(f'\nTerminal Velocity: {terminal_velocity:.2f} m/s')

print('\nMax Height:')
print(plot_label + f': {max_height:.2f} m')

if drag:
    print(plot_label_nodrag + f': {max_height_nodrag:.2f} m')

print('\nFlight Time:')
print(plot_label + f': {flight_time:.2f} s')

if drag:
    print(plot_label_nodrag + f': {flight_time_nodrag:.2f} s')

plt.show()
