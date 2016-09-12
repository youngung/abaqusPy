"""
"""
import abaquspy.plots.plotter
import abaquspy.mats.ifsteel as ifsteel
import numpy as np
## stress-strain curve

ref_dat=ifsteel.hard() ## ((yied stress, eps))
dat=np.array(ref_dat).T
stress=dat[0]; strain=dat[1]
ref_dat=np.array([strain,stress])

abaquspy.plots.plotter.strstr(fn='strstr.txt',ref_dat=ref_dat)
