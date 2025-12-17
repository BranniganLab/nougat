from setuptools import setup

setup(
    name='nougat',
    version='1.0',
    description='Toolkit for analysis of membrane disruption by proteins and other inclusions',
    url='https://github.com/BranniganLab/nougat',
    author='Brannigan Lab',
    author_email='grace.brannigan@rutgers.edu',
    packages=['nougat'],
    install_requires=['numpy>=2.0.0', 'matplotlib>=3.9.1'],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5'
    ],
)
