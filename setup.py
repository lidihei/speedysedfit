from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

install_requires = [
    "numpy",
    "emcee",
    "astropy",
    "astroquery",
    "pyvo",
    "matplotlib",
    "pyaml",
    "corner",
    "tqdm",
    "scipy",
    "h5py",
    "pandas",
    "gaiadr3-zeropoint",
]

VERSION = '0.1.0' 
setup(
    name="speedysedfit",
    version=VERSION,
    author="Jiao Li",
    author_email="lijiao@bao.ac.cn",
    description="MC approach to fit photometric SEDs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lidihei/speedysedfit",
    packages=find_packages(),
    include_package_data=True,
    install_requires=install_requires,
    test_suite='pytest.collector',
    tests_require=['pytest'],
    entry_points={
        'console_scripts': ['speedysedfit=speedysedfit.main:main', 'speedysedfit-batch=speedysedfit.speedyfit_batch:main'],
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
    python_requires='>=3.7',
)
