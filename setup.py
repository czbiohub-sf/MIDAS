from setuptools import setup
from iggtools import version


setup(name='iggtools',
      version=version,
      summary='Integrated Gut Genome Tools',
      description='https://github.com/czbiohub/iggtools/wiki',
      url='http://github.com/czbiohub/iggtools',
      author='Microbiome Team At Pollard Lab and CZ Biohub',
      author_email='bdimitrov@chanzuckerberg.com',
      license='MIT',
      packages=['iggtools', 'iggtools/subcommands', 'iggtools/common', 'iggtools/params'],
      install_requires=[],
      dependency_links=[],
      entry_points={
        'console_scripts': [
          'iggtools = iggtools.__main__:main'
        ]
      },
      zip_safe=False
)
