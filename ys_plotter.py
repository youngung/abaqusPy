import matplotlib.pyplot as plt
import os
import glob
import numpy as np

def comb(fns):
    fig=plt.figure(figsize=(3.5,3))
    ax=fig.add_subplot(111)
    for i in xrange(len(fns)):
        dat=np.loadtxt(fns[i]).T
        label = fns[i].split('.')[0]
        ax.plot(dat[0],dat[1],label=label)

    ax.grid()
    ax.legend(fontsize=6,ncol=4,bbox_to_anchor=(1,1.1))
    ax.set_aspect('equal')
    fnPdf='all_ys.pdf'
    fig.savefig(fnPdf,bbox_inches='tight')

def main(fn='ys.txt'):
    """
    Arguments
    ---------

    fn='ys.txt'
    """
    fig=plt.figure(figsize=(3.5,3))
    ax=fig.add_subplot(111)
    dat=np.loadtxt(fn).T
    ax.plot(dat[0],dat[1])
    ax.grid()
    fnPdf='%s.pdf'%fn.split('.')[0]
    fig.savefig(fnPdf,bbox_inches='tight')
    return fnPdf

if __name__=='__main__':
    fns=glob.glob('hah_*.txt')
    comb(fns)
    # for i in xrange(len(fns)):
    #     fnPdf=main(fns[i])
    #     print '%s->%s'%(fns[i],fnPdf)
