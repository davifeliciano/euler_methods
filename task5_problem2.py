import sys
import euler_methods as em
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
    a, b = em.check_dom(dom)

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
    alpha = abs(float(sys.argv[6]))
except ValueError:
    raise SystemExit(exit_str + '\nAll the arguments must be numbers')
except IndexError:
    alpha = 1.0
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
    print(arg_info[i] + f'= {float(sys.argv[i + 1])}' + arg_unit[i])

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

dom = (0, 120)

vel_funcs = (e_dot_func, n_dot_func, z_dot_func)
vel_init_values = (e_dot_init, n_dot_init, z_dot_init)

t, vel = euler_system(vel_funcs, vel_init_values, dom, 0.1)

pos_funcs = (
    lambda i, t, e, n, z: vel[0][i],
    lambda i, t, e, n, z: vel[1][i],
    lambda i, t, e, n, z: vel[2][i]
)

pos_init_values = (0, 0, 0)

t, pos = euler_system(pos_funcs, pos_init_values, dom, 0.1)

# pick only the subdomain where z is positive
for i in range(len(t)):
    if pos[2][i] < 0:
        t = t[:i]
        pos = [p[:i] for p in pos]
        vel = [v[:i] for v in vel]
        break

print('Ploting graphics')

''' setting up tex in mpl
    you must install texlive-latex-extra texlive-fonts-recommended
    cm-super dvipng ghostscript packages for this '''
plt.rcParams.update({
    'text.usetex': True,
    'figure.figsize': [10.8, 4.8]
})

# setting up figures
pos_fig, pos_axs = plt.subplots(1, 2)

fig_3d = plt.figure()
ax_3d = fig_3d.add_subplot(projection='3d')

# plotting
pos_axs[0].plot(pos[0], pos[1])
pos_axs[1].plot(t, pos[2])

traj = ax_3d.plot(*pos)

print('Saving images...')

pos_fig.savefig('position.png', dpi=300)
plt.close(pos_fig)
plt.show()
