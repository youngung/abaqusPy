
from distutils.core import setup
setup(name='abaqusPy',
      version='0.01',
      description='Collection of Abaqus modules',
      author='Youngung Jeong',
      author_email='youngung.jeong@gmail.com',

      packages=['abaquspy.lib'],
      package_dir={'abaquspy.lib':'lib'}

      )

      # packages=['MP','MP.lib','MP.mat','MP.opt','MP.cal','MP.geom'],
      # package_dir={'MP':'src',
      #              'MP.lib':'src/lib',
      #              'MP.mat':'src/mat',
      #              'MP.opt':'src/opt',
      #              'MP.cal':'src/cal',
      #              'MP.geom':'src/geom'}
      # )
