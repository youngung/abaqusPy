"""
"""
import os, glob
path_site_packages='/home/younguj/anaconda2/lib/python2.7/site-packages/'
os.sys.path.append(path_site_packages)

import abaquspy.plots.plotter
import abaquspy.mats.ifsteel as ifsteel
import numpy as np
## stress-strain curve

ref_dat=ifsteel.hard() ## ((yied stress, eps))
dat=np.array(ref_dat).T
stress=dat[0]; strain=dat[1]
ref_dat=np.array([strain,stress])

fnsRst=glob.glob('OneElement_??_*.txt')
fig=None
rd=None
line_style=['-','--','-.',':','|']
for i in xrange(len(fnsRst)):
    # if i==0: rd=ref_dat
    # else: rd=None
    fn = fnsRst[i]
    fig=abaquspy.plots.plotter.strstr(
        fn=fnsRst[i],ref_dat=rd,fig=fig,
        label=fn.split('OneElement_')[-1],
        ls=line_style[i])

fnFig='OneElementResult.pdf'
fig.tight_layout()
fig.savefig(fnFig,bbox_to_inches='tight')
print 'Flow stress curve saved to %s'%(fnFig)
