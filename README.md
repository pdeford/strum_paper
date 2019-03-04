# DNA shape complements sequence-based representations of transcription factor binding sites

This GitHub repository is associated with the paper _DNA shape complements sequence-based representations of transcription factor binding sites_ by P. DeFord and J. Taylor, and reproduces all of the analysis and figures from that paper.

This analysis relies heavily on the **StruM package**. Source code, installation instructions, and documentation can be found [here](https://github.com/pdeford/StructuralMotifs).

In addition this analysis has the following dependencies:

* [Python 2.7](https://www.python.org/downloads/)
    - [NumPy](http://www.numpy.org/)
    - [Matplotlib](https://matplotlib.org/)
    - [SciPy](https://www.scipy.org/)
    - [Biopython](https://biopython.org/)
    - [scikit-learn](https://scikit-learn.org/stable/)
* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* [SAMtools](http://www.htslib.org/)
* [MEME suite](http://meme-suite.org)

These can all be installed via [`conda`](https://conda.io/docs/). Below is an example of how to set up an appropriate environment via Conda to run this analysis.


## The Working Environment

```
WORKING_DIRECTORY="~/scratch/working_draft"

cd $WORKING_DIRECTORY
git clone https://github.com/pdeford/strum_paper.git
cd strum_paper
mkdir src
cd src
git clone https://github.com/pdeford/StructuralMotifs.git

conda create -n strum_paper python=2.7.15
source activate strum_paper
conda install -c bioconda \
    bedtools=2.27.1 biopython=1.68 meme=4.12.0 \
    regex=2016.06.24 samtools=1.9
conda install -c default \
    matplotlib=2.2.3 numpy=1.15.4 scikit-learn=0.20.1 \
    scipy=1.1.0 cython=0.29.2 python-dateutil=2.7.5 \
    requests=2.20.1 libiconv=1.15
conda install -c conda-forge \
    rdflib=4.2.2
cd StructuralMotifs
python setup.py install
cd ../..
```

## Running the Code

Once your environment is initialize appropriately, all of the code can be produced using the command:

```
./do_all.sh $n_processes
```

where `$n_processes` is the number of processors that you have available to devote to the analysis.
