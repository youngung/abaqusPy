"""
Plotting data generated by hah_test

Youngung Jeong
youngung.jeong@gmail.com
"""
import matplotlib.pyplot as plt
import os
import glob
import numpy as np

import vpscyld.lib_dat


def comb(fns):
    """
    Arguments
    ---------
    fns
    """

    ls=['-','-','-','-','--',':','--','-','--']

    fig=plt.figure(figsize=(7,3))
    ax1=fig.add_subplot(121)
    ax2=fig.add_subplot(122)
    for i in range(len(fns)):
        dat=np.loadtxt(fns[i]).T
        label = fns[i].split('.')[0]
        ax1.plot(dat[0],dat[1],label=label,ls=ls[i])
        ax2.plot(dat[4],dat[5],label=label,ls=ls[i])

    ax1.grid(); ax2.grid()
    ax1.legend(fontsize=6,ncol=3,bbox_to_anchor=(1,1.2))
    ax2.legend(fontsize=6,ncol=3,bbox_to_anchor=(1,1.2))
    ax1.set_aspect('equal'); ax2.set_aspect('equal')
    fnPdf='all_ys.pdf'
    print('Figure saved to <%s>'%fnPdf)

    vpscyld.lib_dat.pi_rad(ax2,150)
    fig.savefig(fnPdf,bbox_inches='tight')

def main(fn='ys.txt'):
    """
    Argument
    --------
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
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f',type=str,help='a list of filenames',default=None)
    args        = parser.parse_args()
    if type(args.f)==type(None):
        hah_fns = glob.glob('hah*.txt')
    else:
        hah_fns = args.f.split()

    comb(hah_fns)
    # for i in xrange(len(fns)):
    #     fnPdf=main(fns[i])
    #     print '%s->%s'%(fns[i],fnPdf)
