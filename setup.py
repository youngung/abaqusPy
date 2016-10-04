#from distutils.core import setup
from numpy.distutils.core import setup

from numpy.distutils.core import Extension
ext1 = Extension(name='yld2000',
                 sources=['umats/yld/yld2000_2d.f',
                          'umats/lib/algb.f',
                          'umats/lib/lib_write.f',
                          'umats/lib/lib.f',
                          'umats/lib/is.f',
                          'umats/lib/cnv.f'
                          ],
                 # f2py_options=['-DF2PY_REPORT_ON_ARRAY_COPY=1']

                 )

setup(name='abaqusPy',
      version='0.01',
      description='Collection of Abaqus python scripts',
      author='Youngung Jeong',
      author_email='youngung.jeong@gmail.com',

      packages=['abaquspy','abaquspy.lib','abaquspy.examples','abaquspy.examples.E8',
                'abaquspy.examples.one','abaquspy.sketches','abaquspy.plots',
                'abaquspy.mats','abaquspy.yld'],
      package_dir={
          'abaquspy':'main',
          'abaquspy.lib':'lib',
          'abaquspy.examples':'examples',
          'abaquspy.examples.E8':'examples/E8',
          'abaquspy.examples.one':'examples/one',
          'abaquspy.sketches':'sketches',
          'abaquspy.plots':'plots',
          'abaquspy.mats':'mats',
          'abaquspy.yld':'umats/yld'
          },

      ext_modules=[ext1]
)
