# DESeq2 Analysis using R and Python for GSE344886 treatment vs control

You can do the same by copying the renv.lock file to your project folder and make sure follow these steps for reproducibility
```r
# Install the renv package if you haven't already
install.packages("renv")

# Restore the renv.lock file to use in your project
renv::restore()
```

In case you are using python, create a virtual environment and install these libraries in them.
```python
pip install pandas numpy matplotlib seaborn scikit-learn pydeseq2
```