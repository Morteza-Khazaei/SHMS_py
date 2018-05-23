import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SHMS_py",
    version="0.0.1",
	py_modules=['shms_py.DisPATCh'],
    author="Morteza Khazaei",
    author_email="Mortezakhazaei1370@gmail.com",
    description="Python library for smap satellite soil moisture data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="http://shms.ut.ac.ir",
	install_requires=['GDAL', 'numpy', 'pymodis',  'requests', 'osr',
					  'numexpr', 'gdalnumeric', 'datetime', 'osgeo',
					  'psycopg2', 'inspect', 'suds', 'xml', 'zipfile'],
	license='GNU GPL 2 or later',
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)