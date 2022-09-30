import numpy as np
import matplotlib.pyplot as plt

def pltfit(ax, x, y, s=None, e=None, label=''):
    fit=np.polyfit(np.log2(x[s:e]), np.log2(y[s:e]), 1)
    print(fit)
    y_fit = (2**fit[1])*x**(fit[0])
    ax.loglog(x, y, '.', label=label)
    ax.loglog(x, y_fit, color='black')

def plots1d():
    se = {'s':10, 'e':-5}

    fig = plt.figure(figsize=(5,4), dpi=100)
    # left bottom width height
    ax1 = fig.add_axes([0.1,0.5,0.8,0.4])
    ax2 = fig.add_axes([0.1,0.1,0.8,0.4])

    ####
    errL = np.fromfile('data/mon6-L-error.data')
    errR = np.fromfile('data/mon6-R-error.data')
    errtrap = np.fromfile('data/mon6-trap-error.data')
    errmid = np.fromfile('data/mon6-mid-error.data')
    errsimp = np.fromfile('data/mon6-simp-error.data')
    numpts = (2**(np.arange(errtrap.size)+1))

    pltfit(ax1,numpts, errL,  label='Left end', **se)
    pltfit(ax1,numpts, errR,  label='Right end', **se)
    pltfit(ax1,numpts, errtrap,  label='Trapezoidal', **se)
    pltfit(ax1,numpts, errmid,  label='Midpoint', **se)
    pltfit(ax1,numpts, errsimp,  label='Simpson', s=3, e=10)
    ax1.set_title('mon')
    ax1.legend()

    ####
    errL = np.fromfile('data/circ-L-error.data')
    errR = np.fromfile('data/circ-R-error.data')
    errtrap = np.fromfile('data/circ-trap-error.data')
    errmid = np.fromfile('data/circ-mid-error.data')
    errsimp = np.fromfile('data/circ-simp-error.data')
    numpts = (2**(np.arange(errtrap.size)+1))

    pltfit(ax2,numpts, errL,  label='Left end', **se)
    pltfit(ax2,numpts, errR,  label='Right end', **se)
    pltfit(ax2,numpts, errtrap,  label='Trapezoidal', **se)
    pltfit(ax2,numpts, errmid,  label='Midpoint', **se)
    pltfit(ax2,numpts, errsimp,  label='Simpson', **se)
    ax2.set_title('circ')
    ax2.legend()

    fig.savefig('1dplots.pdf')

def plots3d():
    fig = plt.figure(figsize=(5,4), dpi=100)
    # left bottom width height
    ax = fig.add_axes([0.1,0.1,0.8,0.8])

    errmid= np.fromfile('data/3-ball-mid-error.data')
    errmc= np.fromfile('data/3-ball-mc-error.data')
    numpts = (2**(np.arange(errmid.size)+1))**2

    pltfit(ax,numpts, errmid,  label='Midpoint')
    pltfit(ax,numpts, errmc,  label='Monte-Carlo')
    ax.set_title('n=3')
    ax.legend()
    fig.savefig('3dplots.pdf')

plots1d()
plots3d()
