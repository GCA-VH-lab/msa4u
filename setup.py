from setuptools import setup
import os

with open("README.md", "r") as fh:
    long_description = fh.read()

def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths

extra_files = package_files("msa4u/msa4u_data")
extra_files.append("docs/pypi.md")

setup(name="msa4u",
      version="0.4.0",
      description="A tool for short uORF annotation.",
      url="https://art-egorov.github.io/msa4u/",
      author="Artyom Egorov",
      author_email="artem.egorov@med.lu.se",
      license="WTFPL",
      packages=["msa4u"],
      package_data={"msa4u": extra_files},
      install_requires=["biopython", "configs", "argparse", "statistics", "pandas", "reportlab"],
      long_description=long_description,
      long_description_content_type="text/markdown",
      scripts=["bin/msa4u"],
      zip_safe=False)