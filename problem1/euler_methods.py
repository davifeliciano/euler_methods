import numpy as np


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


def sec_method(func, dom=(0., 1.), err=10e-12, max_iter=50):
    '''
    finds the root of the equation func = 0 in the interval dom
    returns None if no root is found in max_iter iteractions
    '''
    a, b = check_dom(dom)

    for i in range(max_iter):
        x = b - func(b) * (a - b) / (func(a) - func(b))
        if abs((x - b) / b) < abs(err):
            return x
        a = b
        b = x
    return None


def euler(func, init_value, dom=(0., 1.), h=0.1, **kwargs):
    '''
    solves y' = func(x, y); y(dom[0]) = init_value in the given dom
    using the given type of method as a kwarg 
    (method='explicit', ='implicit', or ='modified')
    '''
    a, b = check_dom(dom)

    methods = ('explicit', 'implicit', 'modified')

    try:
        method = kwargs['method']
    except KeyError:
        method = 'modified'
        pass

    if not isinstance(method, str):
        raise TypeError('method kwarg value must be a string')
    if method not in methods:
        print('Invalid method in euler(). Using modified method instead')

    x = np.arange(a, b + h, h)
    y = np.array([init_value])

    if method == 'explicit':
        for i in range(len(x) - 1):
            y = np.append(y, y[i] + h * func(x[i], y[i]))
        return x, y

    if method == 'implicit':
        for i in range(len(x) - 1):
            new_y = sec_method(
                lambda t: t - y[i] - h * func(x[i + 1], t),
                dom=(y[i], y[i] + h * func(x[i], y[i]))
            )
            y = np.append(y, new_y)
        return x, y

    if method == 'modified':
        for i in range(len(x) - 1):
            f = func(x[i], y[i])
            new_y = sec_method(
                lambda t: t - y[i] - h *
                (f + func(x[i + 1], y[i] + h * f)) / 2,
                dom=(y[i], y[i] + h * f)
            )
            y = np.append(y, new_y)
        return x, y
