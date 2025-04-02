from setuptools import setup, find_packages

setup(
    name="Magpie",
    version="0.1.0",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    url='https://github.com/jmmcd',
    author='James McDermott',
    author_email='jamesmichaelmcdermott@gmail.com',
    description='Magpie: multi-objective archive genetic programming (for symbolic regression) from Ireland',
    install_requires=[
        "numpy",
        "pandas",
        "scikit-learn",
        "matplotlib",
        "sympy"
    ],
    include_package_data=True,
    package_data={
        "Magpie": ["data/*", "grammars/*"]
    },
)