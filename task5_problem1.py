import euler_methods as em
import numpy as np
import matplotlib.pyplot as plt

lbd = 0.1
n_init = 1000
dom = (0., 50.)
h = (5, 1)
def func(x, y): return - lbd * y

''' setting up tex in mpl
    you must install texlive-latex-extra texlive-fonts-recommended
    cm-super dvipng ghostscript packages for this '''
plt.rcParams.update({
    'text.usetex': True,
    'figure.figsize': [8.0, 4.8]
})

fig, ax = plt.subplots(1, 2)

ax[0].set(ylabel=r'$N(t)$')
for i in range(2):
    #plt.title(r'$N_0 = 1000\qquad\lambda = 0.1\qquad\Delta t = 5$')
    t_exp, n_exp = em.euler(func, n_init, dom, h[i], method='explicit')
    t_imp, n_imp = em.euler(func, n_init, dom, h[i], method='implicit')
    t_mod, n_mod = em.euler(func, n_init, dom, h[i], method='modified')

    ax[i].set(
        title=r'$N_0 = 1000\qquad\lambda = 0.1\qquad\Delta t = {h}$'.format(
            h=h[i]),
        xlabel=r'$t$',
        box_aspect=1
    )
    ax[i].grid(ls='--')

    t = np.arange(dom[0], dom[1] + 1)
    ax[i].plot(t, n_init * np.exp(- lbd * t), label='Solução Exata')
    ax[i].plot(t_exp, n_exp, label='Método de Euler Explícito')
    ax[i].plot(t_imp, n_imp, label='Método de Euler Implícito')
    ax[i].plot(t_mod, n_mod, label='Método de Euler Modificado')

    ax[i].legend(loc='upper right')

plt.savefig('problem1.png', dpi=300)
plt.show()
