import euler_methods as em
import numpy as np
from math import sin, cos, pi, radians, sqrt
import sys


def first_quad(deg):
    if abs(deg) == 90:
        return 90
    if deg >= 0:
        return deg % 90
    else:
        return - abs(deg) % 90


usg_str = f'Usage:\n{sys.argv[0]} [latitude] [velocity] [angle] [direction]\nOptional arguments: [mass] [alpha]'

try:
    lat = radians(first_quad(float(sys.argv[1])))
    velocity = float(sys.argv[2])
    angle = radians(abs(first_quad(float(sys.argv[3]))))
    direction = float(sys.argv[4])
except ValueError:
    raise SystemExit(usg_str + '\nAll the arguments must be numbers')
except:
    raise SystemExit(usg_str)

try:
    alpha = float(sys.argv[6])
except ValueError:
    raise SystemExit(usg_str + '\nAll the arguments must be numbers')
except IndexError:
    alpha = 1.0
    pass
except:
    raise SystemExit(usg_str)

try:
    mass = float(sys.argv[4])
except ValueError:
    raise SystemExit(usg_str + '\nAll the arguments must be numbers')
except IndexError:
    mass = 1.0
    pass
except:
    raise SystemExit(usg_str)

omega = 7.292e-5
g = 9.807

print(lat, angle, direction, mass, alpha)

e_dot_init = velocity * cos(angle) * sin(direction)
n_dot_init = velocity * cos(angle) * cos(direction)
z_dot_init = velocity * sin(angle)


def abs_vec(e, n, z): return sqrt(e ** 2 + n ** 2 + z ** 2)


def e_dot(e, n, z):
    return (- alpha * abs_vec(e, n, z) * e + 2 * omega * (n * sin(lat) - z * cos(lat))) / mass


def n_dot(e, n, z):
    return (- alpha * abs_vec(e, n, z) * n + 2 * omega * e * sin(lat)) / mass


def z_dot(e, n, z):
    return (- alpha * abs_vec(e, n, z) * z + 2 * omega * cos(lat)) / mass - g
