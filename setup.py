from setuptools import setup

setup(
    name="nccrPipe",
    version="0.0.1",
    author="Anna Sintsova",
    author_email="ansintsova@ethz.ch",
    description=("Microbial Genomics Pipeline"),
    license="LICENSE",
    keywords="nccrPipe",
    url="https://github.com/MicrobiologyETHZ/NCCR_genomicsPipeline",
    install_requires=[
        'click'],
    packages=['workflow'],
    entry_points={
        'console_scripts': ['nccrPipe=workflow.main:main'],
    }
)

