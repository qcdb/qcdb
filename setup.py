import sys
import setuptools
import versioneer

short_description = "QCDB is a module for quantum chemistry facilitating a uniform interface and interoperable post-processing procedures."

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except FileNotFoundError:
    long_description = short_description

if __name__ == "__main__":
    setuptools.setup(
        name='qcdb',
        description='Quantum Chemistry Common Driver and Databases.',
        author='The QCDB Development Team',
        author_email='psi4aiqc+qcdb.gmail.com',
        url="https://github.com/qcdb/qcdb",
        license='BSD-3C',
        version=versioneer.get_version(),
        cmdclass=versioneer.get_cmdclass(),
        packages=setuptools.find_packages(), #exclude=['*checkup*']),
        include_package_data=True,
        package_data={'': ['qcdb/basis/GENBAS', 'qcdb/basis/*.gbs']},
        setup_requires=[] + pytest_runner,
        python_requires=">=3.9",
        install_requires=["numpy >= 1.12.0", "qcelemental >= 0.22.1", "qcengine >= 0.19.1", "networkx>=2.4.0"],
        extras_require={
#            'docs': [
#                'numpydoc',
#                'sphinx',  # autodoc was broken in 1.3.1
#                'sphinxcontrib-napoleon',
#                'sphinx_rtd_theme',
#            ],
            'tests': [
                'pytest >= 4.0.0',
                'pytest-cov',
                # 'jsonschema',  # needed for speciality `pytest --validate`
            ],
            'align': [
                'networkx>=2.4.0',
            ],
#            'viz': [
#                'nglview',
#            ],
            'lint': [
                'autoflake',
                'black',
                'isort',
            ],
        },
        tests_require=[
            'pytest >= 4.0.0',
            'pytest-cov',
            # 'jsonschema',  # needed for speciality `pytest --validate`
        ],
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 3 :: Only',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.9',
        ],
        zip_safe=False,
        long_description=long_description,
        long_description_content_type="text/markdown")
