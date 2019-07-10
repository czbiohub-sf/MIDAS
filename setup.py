from setuptools import setup
from iggtools import __version__


setup(name='iggtools',
      version=__version__,
      summary='Integrated Gut Genome Tools',
      description='https://github.com/czbiohub/iggtools/wiki',
      url='http://github.com/czbiohub/iggtools',
      author='IdSeq Team @ Chan Zuckerberg Initiative',
      author_email='bdimitrov@chanzuckerberg.com',
      license='MIT',
      packages=['iggtools'],
      install_requires=[],
      dependency_links=[],
      entry_points={
        'console_scripts': [
          'iggtools = iggtools.__main__:main'
        ]
      },
      zip_safe=False
)
