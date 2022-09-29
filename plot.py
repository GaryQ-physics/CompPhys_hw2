import numpy as np
import matplotlib.pyplot as plt

def plots1d():
    def pltfit(ax, x, y, label=''):
        fit=np.polyfit(np.log2(x[10:-5]), np.log2(y[10:-5]), 1)
        print(fit)
        y_fit = (2**fit[1])*x**(fit[0])
        ax.loglog(x, y, '.', label=label)
        ax.loglog(x, y_fit, color='black')
    
    #fig = plt.figure(figsize=(8,8), dpi=300)
    fig = plt.figure(figsize=(5,4), dpi=100)
    # left bottom width height
    ax1 = fig.add_axes([0.1,0.5,0.8,0.4])
    ax2 = fig.add_axes([0.1,0.1,0.8,0.4])
    #ax = fig.add_axes([0.1,0.1,0.9,0.9])
    #ax1=ax2=ax
    
    ####
    errL = np.fromfile('data/mon6-L-error.data')
    errR = np.fromfile('data/mon6-R-error.data')
    errtrap = np.fromfile('data/mon6-trap-error.data')
    errmid = np.fromfile('data/mon6-mid-error.data')
    numpts = (2**(np.arange(errtrap.size)+1))
    
    pltfit(ax1,numpts, errL,  label='Left end')
    pltfit(ax1,numpts, errR,  label='Right end')
    pltfit(ax1,numpts, errtrap,  label='Trapezoidal')
    pltfit(ax1,numpts, errmid,  label='Midpoint')
    ax1.set_title('mon')
    ax1.legend()
    
    ####
    errL = np.fromfile('data/circ-L-error.data')
    errR = np.fromfile('data/circ-R-error.data')
    errtrap = np.fromfile('data/circ-trap-error.data')
    errmid = np.fromfile('data/circ-mid-error.data')
    numpts = (2**(np.arange(errtrap.size)+1))
    
    pltfit(ax2,numpts, errL,  label='Left end')
    pltfit(ax2,numpts, errR,  label='Right end')
    pltfit(ax2,numpts, errtrap,  label='Trapezoidal')
    pltfit(ax2,numpts, errmid,  label='Midpoint')
    ax2.set_title('circ')
    ax2.legend()

    fig.savefig('1dplots.pdf')

def plots3d():
    def pltfit(ax, x, y, label=''):
        fit=np.polyfit(np.log2(x[5:-2]), np.log2(y[5:-2]), 1)
        print(fit)
        y_fit = (2**fit[1])*x**(fit[0])
        ax.loglog(x, y, '.', label=label)
        ax.loglog(x, y_fit, color='black')

    fig = plt.figure(figsize=(5,4), dpi=100)
    # left bottom width height
    ax = fig.add_axes([0.1,0.1,0.8,0.8])

    err3d= np.fromfile('data/3-ball-error.data')
    numpts = (2**(np.arange(err3d.size)+1))**2
    print(numpts)
    pltfit(ax,numpts, err3d,  label='3D')

    fig.savefig('3dplots.pdf')

plots1d()
plots3d()
