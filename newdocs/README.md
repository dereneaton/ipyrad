
## README


### Writing the docs

Tips on using REST styling here:
https://docs.typo3.org/m/typo3/docs-how-to-document/master/en-us/WritingReST/MenuHierarchy.html

Try to make the style of notebooks match .rst files in terms of the use of headers, etc. 

### Building the docs (https://nbsphinx.readthedocs.io/en/0.4.3/)

1. Install nbsphinx
2. Edit your conf.py and add 'nbsphinx' to extensions.
3. Edit index.rst to add .rst and .ipynb files to the toctree.
4. Run Sphinx!

```bash
cd newdocs/
sphinx-build . html
```

### Using readthedocs

We do not need to be able to RUN the notebooks again, just to use nbconvert to make them look pretty. Working continuous integration into the building of docs is too messy (and slow).

requirements.txt
	- sphinx>=1.4
	- ipykernel
	- nbsphinx
	- sphinx_rtd_theme