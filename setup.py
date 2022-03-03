from setuptools import setup
from iggtools import version


setup(name='iggtools',
      version=version,
      description='https://github.com/czbiohub/iggtools/wiki',
      url='http://github.com/czbiohub/iggtools',
      author='Chunyu Zhao, Boris Dimitrov',
      author_email='chunyu.zhao@czbiohub.org',
      license='MIT',
      packages=['iggtools', 'iggtools/subcommands', 'iggtools/common', 'iggtools/params', 'iggtools/models'],
      install_requires=[],
      dependency_links=[],
      entry_points={
        'console_scripts': [
          'iggtools = iggtools.__main__:main'
        ]
      },
      zip_safe=False
)
