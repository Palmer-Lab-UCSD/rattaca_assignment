from setuptools import setup, find_packages

setup(
    name='rattaca_assign',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
    ],
    entry_points={
        'console_scripts': [
            'rattaca_assign=rattaca_assign.cli:main',
        ],
    },
    author='Benjamin Johnson (Palmer Lab), bbjohnson@health.ucsd.edu',
    description='A tool for assigning HS West rats to RATTACA projects',
    python_requires='>=3.6',
)