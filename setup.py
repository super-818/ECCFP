from setuptools import setup

setup(
    name='eccfp',
    version='1.0.0',
    packages=['ECCFP', 'ECCFP.cli'],
    url='',
    license='',
    author='Li Wang, Zhang Tangxuan',
    author_email='',
    description='',
entry_points={
        "console_scripts": [
            "eccfp=ECCFP.main:main",
        ],
    },
)
