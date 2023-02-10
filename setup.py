from setuptools import setup, find_packages
import re

if __name__ == '__main__':
    with open("README.md", "r", encoding="utf-8") as fh:
        long_description = fh.read()

    with open("requirements.txt", "r") as f:
        requirements = f.read().splitlines()

    with open("src/variant_extractor/__init__.py", "r") as fd:
        init_content = fd.read()
    version = re.search(
        r'^__version__\s*=\s*[\'\"]([^\'\"]*)[\'\"]', init_content, re.MULTILINE
    ).group(1)
    author = re.search(
        r'^__author__\s*=\s*[\'\"]([^\'\"]*)[\'\"]', init_content, re.MULTILINE
    ).group(1)

    setup(
        name='variant-extractor',
        version=version,
        author=author,
        author_email='contact@rodrigomartin.dev',
        description='Deterministic and standard extractor of indels, SNVs and structural variants (SVs) from VCF files',
        keywords='vcf genetics bioinformatics variant indel snv sv',
        long_description=long_description,
        long_description_content_type='text/markdown',
        url='https://github.com/EUCANCan/variant-extractor',
        project_urls={
            'Bug Tracker': 'https://github.com/EUCANCan/variant-extractor/issues',
        },
        classifiers=[
            'Programming Language :: Python :: 3 :: Only',
            'Development Status :: 5 - Production/Stable',
            'Operating System :: OS Independent'
        ],
        package_dir={'': 'src'},
        packages=find_packages(where="src"),
        python_requires='>= 3.6',
        install_requires=requirements,
        extras_require={
            "docs": ["sphinx", "sphinx-rtd-theme", "myst_parser", "docutils>=0.18.0"],
        }
    )
