

## Convert jupyter notebooks to rtf documentation
jupyter nbconvert test_rad.ipynb --to rst ../docs/examples_rad.rst
jupyter nbconvert test_gbs.ipynb --to rst ../docs/examples_gbs.rst

## Create Quick-guides
jupyter nbconvert API_Quick_guide.ipynb --to rst ../docs/Quickguide_API.rst


