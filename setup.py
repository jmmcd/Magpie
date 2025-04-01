from setuptools import setup, find_packages

setup(
    name="GPAFS",
    version="0.1.0",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    url='https://github.com/jmmcd',
    author='James McDermott',
    author_email='jamesmichaelmcdermott@gmail.com',
    description='Grammatical Pareto Archive Function Synthesis',
    install_requires=[
        "numpy",
        "pandas",
        "scikit-learn",
        "matplotlib",
        "sympy"
    ],
    include_package_data=True,
    package_data={
        "GPAFS": ["data/*", "grammars/*"]
    },
)