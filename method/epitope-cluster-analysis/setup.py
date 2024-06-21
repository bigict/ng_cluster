import os
from setuptools import setup

# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

with open(os.path.join(os.path.dirname(__file__), 'README.md')) as readme:
    README = readme.read()


def recursive_walk(startdir):
    file_list = []
    for dirname, subdirList, fileList in os.walk(startdir, followlinks=True):
        for fname in fileList:
            file_list.append(os.path.join(dirname, fname))
            for sdname in subdirList:
                sdfiles = recursive_walk(sdname)
                file_list.extend(sdfiles)
    return file_list

pkg_data_files = []
os.chdir('epitope_cluster_analysis')
for subdir in ['.',]:
    file_list = recursive_walk(subdir)
    pkg_data_files.extend(file_list)
os.chdir(os.pardir)


setup(
    name="epitope_cluster_analysis",
    version="2.0.3",
    author='Jason Yan',
    author_email='jyan@lji.org',
    install_requires=[
          'biopython>=1.67',
          'networkx>=2.2',
          'pandas>=1.1.4',
          'matplotlib>=3.1.1'
      ],
    packages=['epitope_cluster_analysis' ],
    package_data={
        'epitope_cluster_analysis': pkg_data_files
    },
    test_suite = 'tests',
    description='PyPA package for mhci smm and smmpmbec prediction method.',
    long_description=README,
    # Important only if the package will be widely distributed.  See more at:
    #    https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Intended Audience :: Developers',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
    ]
)
