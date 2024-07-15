from setuptools import setup
#from midas import version

setup(name='midasv3',
      version='1.0.1',
      description='https://midas2.readthedocs.io/en/latest/',
      url='https://github.com/czbiohub/MIDAS',
      author='Chunyu Zhao and Boris Dimitrov',
      author_email='chunyu.zhao@czbiohub.org',
      license='MIT',
      packages=['midas', 'midas/subcommands', 'midas/common', 'midas/params', 'midas/models'],
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
          'midas = midas.__main__:main'
        ]
      },
      zip_safe=False
)
