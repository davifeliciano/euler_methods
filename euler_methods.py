import numpy as np


def sec_method(func, dom=(0., 1.), err=10e-12, max_iter=50):
    '''
    finds the root of the equation func = 0 in the interval dom
    returns None if no root is found in max_iter iteractions
    '''
    if not isinstance(dom, tuple):
        raise TypeError('dom must be a 2-element tuple')

    if dom[0] < dom[1]:
        a = dom[0]
        b = dom[1]
    else:
        a = dom[1]
        b = dom[0]

    for i in range(max_iter):
        x = b - func(b) * (a - b) / (func(a) - func(b))
        if abs((x - b) / b) < abs(err):
            return x
        a = b
        b = x
    return None


def euler(func, init_value, dom=(0., 1.), h=0.1, method='explicit'):
    '''
    solves y' = func(x, y); y(dom[0]) = init_value in the given dom
    using the given type of method (explicit, implicit, or modified)
    '''
    if not isinstance(dom, tuple):
        raise TypeError('dom value must be a 2-element tuple')
    if dom[0] >= dom[1]:
        raise ValueError('Invalid domain')

    x = np.arange(dom[0], dom[1] + h, h)
    y = np.array([init_value])

    if not isinstance(method, str):
        raise TypeError('method value must be a string')

    if method == 'explicit':
        for i in range(len(x) - 1):
            y = np.append(y, y[i] + h * func(x[i], y[i]))
        return x, y
    elif method == 'implicit':
        for i in range(len(x) - 1):
            new_y = sec_method(
                lambda t: t - y[i] - h * func(x[i + 1], t),
                dom=(y[i], y[i] + h * func(x[i], y[i]))
            )
            y = np.append(y, new_y)
        return x, y
    elif method == 'modified':
        for i in range(len(x) - 1):
            f = func(x[i], y[i])
            new_y = sec_method(
                lambda t: t - y[i] - h *
                (f + func(x[i + 1], y[i] + h * f)) / 2,
                dom=(y[i], y[i] + h * f)
            )
            y = np.append(y, new_y)
        return x, y
    else:
        raise ValueError('Invalid method. Try explicit, implicit or modified')
