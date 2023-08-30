from setuptools import setup

setup(
    name='plind',
    version='0.1.0',
    packages=['plind'],
    install_requires=[
        'requests',
        'importlib; python_version == "3.8"',
    ],
)