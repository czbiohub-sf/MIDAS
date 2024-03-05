from setuptools import setup
from midas2 import version

setup(name='midas2',
      version=version,
      description='https://midas2.readthedocs.io/en/latest/',
      url='https://github.com/czbiohub/MIDAS2',
      author='Chunyu Zhao and Boris Dimitrov',
      author_email='chunyu.zhao@czbiohub.org',
      license='MIT',
      packages=['midas2', 'midas2/subcommands', 'midas2/common', 'midas2/params', 'midas2/models'],
      install_requires=[
        'biopython <= 1.83',
        'numpy >= 1.26.4',
        'pysam >= 0.22.0',
        'gffutils >= 0.12',
        'pandas >= 2.2.0',
        'pybedtools >= 0.9.1',
      ],
      dependency_links=[],
      entry_points={
        'console_scripts': [
          'midas2 = midas2.__main__:main'
        ]
      },
      zip_safe=False
)
