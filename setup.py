from distutils.core import setup
setup(name='abaqusPy',
      version='0.01',
      description='Collection of Abaqus modules',
      author='Youngung Jeong',
      author_email='youngung.jeong@gmail.com',

      packages=['abaquspy','abaquspy.lib','abaquspy.examples','abaquspy.examples.E8',
                'abaquspy.examples.one','abaquspy.sketches','abaquspy.plots',
                'abaquspy.mats'],
      package_dir={
          'abaquspy':'main',
          'abaquspy.lib':'lib',
          'abaquspy.examples':'examples',
          'abaquspy.examples.E8':'examples/E8',
          'abaquspy.examples.one':'examples/one',
          'abaquspy.sketches':'sketches',
          'abaquspy.plots':'plots',
          'abaquspy.mats':'mats'
          }
)
