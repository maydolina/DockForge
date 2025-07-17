from setuptools import setup, find_packages

setup(
    name="DockForge",
    version="1.0",
    description="Preparation package for molecular docking",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Mina Aybüke Yakıcı",
    author_email="myakici27@my.uaa.k12.tr",
    url="https://github.com/maydolina/DockForge",
    packages=find_packages(where="src"),
    python_requires=">=3.7",
    install_requires=["openbabel", "pdbfixer", "openmm"],
    keywords= ["molecule", "docking", "molecular docking", "bioinformatics"]
)
