from setuptools import setup, find_packages

__version__ = '0.1.0'

if __name__ == '__main__':
    with open("README.md", "r", encoding="utf-8") as fh:
        long_description = fh.read()

    setup(
        name='variant-extractor',
        version=__version__,
        author='Rapsssito',
        author_email='contact@rodrigomartin.dev',
        description='Extractor of INDELs, SNVs and SVs from VCF files',
        long_description=long_description,
        long_description_content_type='text/markdown',
        url='https://github.com/Rapsssito/variant-extractor',
        project_urls={
            'Bug Tracker': 'https://github.com/Rapsssito/variant-extractor/issues',
        },
        classifiers=[
            'Programming Language:: Python : : 3'
            'License:: OSI Approved : : MIT License'
            'Operating System: : OS Independent'
        ],
        package_dir={'': 'src'},
        packages=find_packages(where="src"),
        python_requires='>= 3.6',
        install_requires=['pysam>=0.11.2.2']
    )
