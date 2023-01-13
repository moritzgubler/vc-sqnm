from setuptools import setup

setup(
    name='sqnm',
    version='1.0.0',    
    description='Implementation of the SQNM optimization method',
    url='https://github.com/moritzgubler/vc-sqnm',
    author='Moritz Guber',
    author_email='moritz.gubler@e.email',
    license='BSD 2-clause',
    packages=['sqnm'],
    install_requires=[
                    'numpy',   
                      ]
)
