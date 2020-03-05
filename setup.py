from setuptools import setup, find_packages

files = ["model/data/*"]

setup(
    name = 'mewpy',
    version = '0.0.2',
    package_dir = {'':'src'},
    packages = find_packages('src'),
    include_package_data=True,
    install_requires = [
        'inspyred',
        'reframed',
	    'jmetalpy'],
    author = 'Vítor Pereira',
    author_email = 'vitor.pereira@algoritmi.uminho.pt',
    description = 'mewpy - Metabolic Engineering in Python ',
    license = 'Apache License Version 2.0',
    keywords = 'strain optimization',
    url = '',
    package_data = {"":["*.xml","*.csv","*.txt"],'mewpy' : files },
    long_description = open('README.rst').read(),
    classifiers = [
        'Topic :: Utilities',
        'Programming Language :: Python :: 3.7',
    ],
)
