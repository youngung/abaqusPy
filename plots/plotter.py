from MP.lib import checkX
if checkX.main()!=0:
    import matplotlib as mpl
    mpl.use('Agg')

import numpy as np
import matplotlib.pyplot as plt

def strstr(fn):
    """
    Argument
    --------
    fn
    """
    dat=np.loadtxt(fn).T

    fig=plt.figure(figsize=(3.5,3))
    ax=fig.add_subplot(111)
    ax.plot(dat[0],dat[1])
    ax.set_xlabel(r'$E_{11}$')
    ax.set_ylabel(r'$S_{11}$')
    fnFig='%s.pdf'%fn.split('.')[0]
    fig.savefig(fnFig,bbox_to_inches='tight')
    print 'Flow stress curve from %s saved to %s'%(fn,fnFig)
