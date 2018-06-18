{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {
    "collapsed": true
   },
   "source": [
    "## Cookbook: *PCA* analyses\n",
    "\n",
    "As part of the `ipyrad.analysis` toolkit we've created convenience functions for easily performing exploratory principal component analysis (PCA) on your data. PCA is a very standard dimension-reduction technique that is often used to get a general sense of how samples are related to one another. PCA has the advantage over STRUCTURE type analyeses in that it is very fast. Similar to STRUCTURE, PCA can be used to produce simple and intuitive plots that can be used to guide downstream analysis. There are three very nice papers that talk about the application and interpretation of PCA in the context of population genetics:\n",
    "\n",
    "[Reich et al (2008) Principal component analysis of genetic data](https://www.nature.com/articles/ng0508-491)\n",
    "\n",
    "[Novembre & Stephens (2008) Interpreting principal component analyses of spatial population genetic variation](https://www.nature.com/articles/ng.139)\n",
    "\n",
    "[McVean (2009) A genealogical interpretation of principal components analysis](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000686)\n",
    "\n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### A note on Jupyter/IPython\n",
    "This is a Jupyter notebook, a reproducible and executable document. The code in this notebook is Python (2.7), and should be executed either in a jupyter-notebook, like this one, or in an IPython terminal. Execute each cell in order to reproduce our entire analysis. The example data set used in this analysis is from the [empirical example ipyrad tutorial](http://ipyrad.readthedocs.io/pedicularis_.html)."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Required software\n",
    "You can easily install the required software for this notebook with a locally installed `conda` environment. Just run the commented code below in a terminal. If you are working on an HPC cluster you **do not need** administrator privileges to install the software in this way, since it is only installed locally."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {},
   "outputs": [],
   "source": [
    "## conda install ipyrad -c ipyrad\n",
    "## conda install -c conda-forge scikit-allel"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Import Python libraries\n",
    "The call to `%matplotlib inline` here enables support for plotting directly inside the notebook."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {},
   "outputs": [],
   "source": [
    "%matplotlib inline\n",
    "import ipyrad\n",
    "import ipyrad.analysis as ipa      ## ipyrad analysis toolkit"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Quick guide (tldr;)\n",
    "The following cell shows the quickest way to results. Further explanation of all of the features and options is provided further below. "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "loading Assembly: rad\n",
      "from saved path: /tmp/ipyrad-test/rad.json\n",
      "Using default cmap: Spectral\n"
     ]
    },
    {
     "data": {
      "text/plain": [
       "<matplotlib.axes._subplots.AxesSubplot at 0x7fb6fdf82050>"
      ]
     },
     "execution_count": 2,
     "metadata": {},
     "output_type": "execute_result"
    },
    {
     "data": {
      "image/png": "iVBORw0KGgoAAAANSUhEUgAAAagAAAFgCAYAAADuCe0ZAAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAAALEgAACxIB0t1+/AAAADl0RVh0U29mdHdhcmUAbWF0cGxvdGxpYiB2ZXJzaW9uIDIuMi4yLCBodHRwOi8vbWF0cGxvdGxpYi5vcmcvhp/UCwAAIABJREFUeJzt3Xt4VNW9//HPNyEkIYKIREFAEy9EqWiVHPV4afGIilZFbfFHqxZbvICKYkVEawWqtgiitlixVERULPZYKz4t2oK3qv2pDd5AEaSaCggabxUhhIR8zx+zBwMGhkAme03m/XoeHmavvWfz3Y7kw9p7zVrm7gIAIDQ5cRcAAEBjCCgAQJAIKABAkAgoAECQCCgAQJAIKABAkAgoAECQCCgAQJAIKABAkNrEXUBz6Ny5s5eUlMRdBgA0yfz58z929+K46whVqwiokpISVVRUxF0GADSJmf077hpCxi0+AECQCCgAQJAIKABAkFrFMygAaC3mz5+/W5s2be6WdKBafyeiXtLCurq68/v06fPR5jsJKAAISJs2be7u0qXLAcXFxZ/l5OS06gX76uvrraqqqteqVavulnTa5vtbezoDQKY5sLi4+IvWHk6SlJOT48XFxf9Rorf49f0tXM8mzKyjmT1sZm+b2SIz+28z62Rmc83snej3XeKsEQBaWE42hFNSdK2NZlHcPahfSXrC3feXdLCkRZJGS3rS3feT9GS0DQDIMrEFlJl1kPQtSdMkyd3Xu/vnkgZImhEdNkPS6fFUCADhW7zo7bZDB3y/dPDh/XsOHfD90sWL3m4bd01JEyZMKO7Zs2ev/fffv1efPn3K5s+fX9CU98c5SGJvSVWSppvZwZLmS7pc0u7uvlKS3H2lme3W2JvN7EJJF0rSnnvu2ezFVVZWasq0yVpX+6UK8nbSsCHDxXRKAEKyeNHbba8+8eye/Ze1z8+3XNX4Sl396tlFN/915pKyA/ZfH3d9559//iejRo2qkqSZM2fuPGLEiB7PPffcO9v6/jhv8bWRdKikKe5+iKQ1asLtPHef6u7l7l5eXNy8U1lVVlZqzC2j1Hdwmc647Ej1HVymMbeMUmVlZbP+OQCwI24bPa5bMpwkKd9y1X9Z+/zbRo/rtiPnXbx4cdvS0tJvnHnmmSU9e/bs1b9//71Xr16dM3v27PYHHHBAr549e/YaOHBgSXV1tUlSt27deg8bNqxb7969D+jdu/cBCxcuzJekTp061SfP+eWXX+aaWZPqiDOglkta7u4vRdsPKxFYH5pZV0mKfv/a2Ph0mzJtsgaN6KfCdvmSpMJ2+Ro0op+mTJvc0qUAwBZVr/osLxlOSfmWq+oPP8vb0XNXVlYWDB06tGrJkiVvtW/fvv6GG27Y/aKLLip96KGH/rVkyZK36urqNHHixI29gw4dOmxYsGDBoosuuuij4cOH90i2//KXvyzu0aPHgWPGjOn+m9/85v2m1BBbQLn7KknLzKwsajpO0luSHpM0OGobLGl2S9e2rvbLjeGUVNguX+tq17R0KQCwRYVddqmt8Q2btNX4BhXuvkvtjp67S5cu60844YQ1knTuued+8uyzz7bv3r17zUEHHVQjSeedd94nzz//fPvk8YMHD/5Uki644IJPX3311Z2S7ddcc03VsmXLFo4dO3b5mDFjujalhrhH8Q2XNNPM3pD0TUm/kDRe0vFm9o6k46PtFlWQt5Oq19Zs0la9tkYbarJm5CeADHDF+DErnuixuiYZUjW+QU/0WF1zxfgxK3b03E2+HZfzVZyY2dd+WF5wwQWfzp07t2OTztmkCpqZu78WPUc6yN1Pd/fP3P0Tdz/O3feLfv+0pesaNmS4Hpj4+MaQql5bo+kT/qiVnyznORSAYJQdsP/6m/86c8lrp3X99JnD8la/dlrXT5trgMTKlSvbzps3r0iSHnzwwU59+/b9YsWKFW2Tz5fuu+++XY855pjVyePvu+++TpI0bdq0XQ455JA1krRgwYKNt6Ieeuihnffaa69N/+WfAlMdNaKkpETt23bWg3f8WTk5OcrJMZ15/gnauVN7TZk2WTffMCnuEgFAUiKk7pr9+/ea+7x77733unvuuWfXiy++eK/S0tKa3/3ud8uOPPLINQMHDtxnw4YNOvjgg9eOHDmyKnl8TU2NHXTQQfvX19fbrFmz3pWkW2+9dbfnnnuuQ5s2bXznnXeuu/fee5tUJwG1BflFORp02Xe/1s5zKADZICcnRw8++OAmgxoGDBiwesCAAW81dvzIkSOrJk2atLJh2/Tp05ftUA078ubWbEvPoQryimKqCACyCwG1BcOGDNes2+dt8hxq1u3zNGzI8JgrA4D0KisrW//OO++8ua3Hr1ixYkHXrl3rmrsObvFtQUlJicaNnBDNJrFGBXlFGjdyArNJAEALIaC2oqSkhAERABATbvEBAIJEQAEAgsQtPgDIYIsXv932njvGd1Pt6jzlta/98aWjV5SVxT+TuSSNHTt29/vvv79zbm6u77rrrnUzZsyo7Nmz5zbXRg8KADLU4sVvt73z58N6XndyUadfDurR/rqTizrd+fNhPRcvDmNNqD59+qx97bXXFi1ZsuSt008//bMrrriie1PeT0ABQIa6547x3W48p1d+UWFi8vKiwjzdeE6v/HvuGB/Echunnnrq6vbt29dL0tFHH/3lypUrmxScBBQAZKra1XnJcEoqKsyT1q8OZrmNpN/+9rfF/fr1+09TaiCgACBT5bWvXVO96coaa6prpbbtg1luQ5LuvPPOTq+//nq7cePGrWpKDQQUAGSoH186esV1D7xVkwypNdW1uu6Bt2p+fOnoYJbbePTRR9vfcsstXefMmbO0sLCwSWsWEVAAkKHKyvZff/H1U5bcOGfNp9f8ftnqG+es+fTi66csaY5RfM2x3MYLL7xQOHz48L1mz569tFu3bk2eColh5gCQwcrK9l9/8+SmLWOxLZpjuY2rrrqqx9q1a3MHDhy4jyTtscce65966qml21oDAQUA+JrmWG7jH//4x5IdqmFH3gwAQLrQgwIAbGJ7lttIRx30oAAAQSKgAABBIqAAAEEioAAAQSKgACCDLV68uO3wK4eWDr3ivJ7Drxxaunjx4iBmMpekxx9/fKdevXod0KZNmz7Tp0/fpanvJ6AAIEMtXry47c9vHd3z5GGHdhp01bHtTx52aKef3zq6Zyghtffee6+fPn165amnnvrJ9ryfgAKADHXH1Nu6nXPVSfmF7fIlSYXt8nXOVSfl3zH1tiCW2ygrK1t/+OGHVzecp68pCCgAyFC19evykuGUVNguX+vrq4NbbmN7EFDNoLKyUlf/7EpdPvoiXf2zK1VZWRl3SQCyQF5OQW312ppN2qrX1qhtTmFQy21sLwJqB1VWVmrMLaPUd3CZzrjsSPUdXKYxt4wipACk3aUXXrHigYmP1yRDqnptjR6Y+HjNpRdeEcxyGzuCgNpBU6ZN1qAR/dTwHvCgEf00ZdrkmCsD0NqVlZWtv/4n45fMmfLKp7+f+NTqOVNe+fT6n4xfUlZWFsRyGzuKufh20LraL9XYPeB1tc3y+QDAVpWVla2fPOmuIJfbePbZZ9udddZZ+37xxRe5Tz75ZMebbrppj6VLl27zHH8E1A4qyNtJ1WtrNgmp6rU1KsgrirEqANgxzbHcxre//e21H3744RvbXcP2vhEJw4YM16zb56nhPeBZt8/TsCHDY64MADIbPagdVFJSonEjJ2jKtMlaV7tGBXlFGjdygkpKSuIuDQC2SyjLbRBQzaCkpEQ33zAp7jIAtA719fX1lpOT0ywj4UJXX19vkuob28ctPgAIy8Kqqqqdox/crVp9fb1VVVXtLGlhY/vpQQFAQOrq6s5ftWrV3atWrTpQrb8TUS9pYV1d3fmN7SSgACAgffr0+UjSaXHXEYLWns4AgAxFQAEAgkRAAQCCFHtAmVmumb1qZn+OtkvN7CUze8fMHjKzIBbeAgC0rNgDStLlkhY12L5Z0m3uvp+kzyQNiaUqAECsYg0oM+su6TuS7o62TdL/SHo4OmSGpNPjqQ4AEKe4e1C3Sxqlr75FvKukz929LtpeLqnRpYvN7EIzqzCziqqqqsYOAQBksNgCysxOkfSRu89v2NzIoY1O9+HuU9293N3Li4uLGzsEAJDB4vyi7lGSTjOzkyUVSOqgRI+qo5m1iXpR3SV9EGONAICYxNaDcvdr3L27u5dIGiTpKXc/W9LTkr4XHTZY0uyYSgQAxCjuZ1CNuVrST8xsqRLPpKbFXA8AIAZBzMXn7s9IeiZ6/a6kw+KspzlVVlZGa0V9qYK8nTRsyHDWigKAbRBiD6rVqKys1JhbRqnv4DKdcdmR6ju4TGNuGaXKysq4SwOA4BFQaTRl2mQNGtFPhe3yJUmF7fI1aEQ/TZk2OebKACB8BFQarav9cmM4JRW2y9e62jUxVQQAmYOASqOCvJ1UvbZmk7bqtTUqyCuKqSIAyBwEVBoNGzJcs26ftzGkqtfWaNbt8zRsyPCYKwOA8AUxiq+1Kikp0biRE6JRfGtUkFekcSMnMIoPALaBuTc6k1BGKS8v94qKirjLAIAmMbP57l4edx2h4hYfACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIsQWUmfUws6fNbJGZvWlml0ftncxsrpm9E/2+S1w1AgDiE2cPqk7Sle5+gKQjJF1iZr0kjZb0pLvvJ+nJaBsAkGViCyh3X+nur0SvV0taJKmbpAGSZkSHzZB0ejwVAgDiFMQzKDMrkXSIpJck7e7uK6VEiEnabQvvudDMKsysoqqqqqVKBQC0kNgDysx2kvRHSSPc/YttfZ+7T3X3cncvLy4uTl+BAIBYxBpQZpanRDjNdPdHouYPzaxrtL+rpI/iqg8AEJ84R/GZpGmSFrn7rQ12PSZpcPR6sKTZLV0bACB+bWL8s4+SdK6kBWb2WtR2raTxkv5gZkMkvS9pYEz1AQBiFFtAufvzkmwLu49ryVoAAOGJfZAEAACNIaAAAEEioAAAQSKgAABBIqAAAEEioAAAQWpSQJlZkZnlpqsYAACSthpQZpZjZj8ws7+Y2UeS3pa0Mlq/aaKZ7dcyZQIAsk2qHtTTkvaRdI2kLu7ew913k3SMpBcljTezc9JcIwAgC6WaSaKfu9du3ujunyoxyesfowlfAQBoVlsNqM3DycwKJJ0jqVDSg+7+SWMBBgDAjmrqKL5fScqVtE7So81fDgAACakGSTxoZvs0aOokaaak30vaJZ2FAQCyW6pnUNdJutHMPpB0g6RblFivqUDS2PSWBgDIZqmeQb0r6QdmdrSkhyT9RdLx7r6hJYoDAGSvVLf4djGzSyT1knSWpP9I+quZndISxQEAsleqQRKPSqpR4pbe/e5+n6RTJfUxs8fSXRwAIHulega1q6QHlRhW/kNJcvdqSePMrGuaawMAZLFUAXW9pLmSNkga3XCHu69MV1EAAKQaJPGIpEdaqBYAADbaakCZWY6kwZK+K6mHpDpJ70i6y92fSXt1AICsleoW3zRJ/5Y0XtL3JH0h6TlJ15lZb3efnOb6AABZKlVA9XH3H0WvnzezF939ejP7u6TXJBFQAIC0SDXMvDY51ZGZHSppvSS5e40kT3NtAIAslqoHdZWkp81snaQ8SYMkycyKJf05zbUBALJYqlF8T5nZXpJ2dfePG7RXSRqV7uIAANkr5XIbnvDx5u1m1iU9JQEA0PT1oBqa1mxVAACwme0OKHf/TnMWAgBAQ00OKDPrlI5CAABoKNVyG9c1eN3LzJZImm9mlWZ2eNqrAwBkrVQ9qDMbvJ4o6XJ3L1Vibajb0lYVACDrNeUW3x7u/rgkufvLSizBAQBAWqT6ou7e0cKEJqm7mbVz97XRvrz0lgYAyGapAmrAZts5kmRmu0uakpaKAABQ6pkknt1C+4eSfpOWigAA0A58D8rMpjZnIQAANJRqwcItfefJJJ3c/OUAAJCQ6hlUlRILFlqDNo+2d0tXUQAApAqodyUd5+7vb77DzJalpyQAAFI/g7pd0i5b2DehmWsBAGCjVKP4tjhSz91Z7h0AkDap5uI7OsX+DmZ2YPOWtPHc/c1ssZktNbPR6fgzAADhSvUM6rtmNkHSE5LmKzFookDSvpKOlbSXpCubuygzy1Xie1bHS1ou6Z9m9pi7v9XcfxYAIEypbvFdYWa7SPqepIGSukqqlrRI0m/d/fk01XWYpKXu/q4kmdksJWa1IKAAIEuk6kHJ3T+T9LvoV0vpJqnhKMHlkjZZ3sPMLpR0oSTtueeeLVcZAKBF7MiS7+lkjbT5JhvuU9293N3Li4uLW6gsAEBLSdmDislyST0abHeX9EFMtQBAoyrfe0+//tkvtGbFxyrq1lmX3XCtSkpL4y6r1Qg1oP4paT8zK5W0QtIgST+ItyQA2aJh8GzYua1yPUf6Yt3GEJKkX464Vgv/9oI6rcvRt7SH2uszXfXiDzRx7oOEVDMxd9/6AWYdJBW7+782az/I3d9IW2FmJyvxReFcSfe4+01bOra8vNwrKirSVQqALFL53nu66vgf6Lh/5Wu11muelusM7a18y1WNb9Bf9/pSNfX1Om1Zh41tj+pdHafuaq+2WnL2frr1gW17ZG9m8929PM2XlLFSfQ/qLElvS/qjmb1pZv/VYPe96SzM3ee4e09332dr4QQAzenXP/uFjvtXvvItVy9o1cZwkqR8y9WJ/95Jucs+3aTtdO2tF7RK+ZarNR98Emf5rUqqQRLXSurj7t+U9CNJ95vZmdG+xgYyAEBGW7Pi443h4/KNr5PyLVc5m/34y7dcuVw1vkFFe+zaYrW2dqmeQeW6+0pJcveXzexYSX82s+7abFQdALQGRd06q8Y/U77lymSq8Q2bhFSNb1D9Zj/+km1P7lOjidEzKuy4VD2o1Wa2T3IjCqu+Snxp9htprAsAYnHZDdfqyX1qVOMbdJS66E96VzW+QZI2PoPa0KPTJm2zCiu1+4DDGSDRzFL1oIZps1t57r7azPpLOittVQFATEpKSzVx7oOJUXwffKIeHXroVc+RVteoaI9ddVvUQ0ruL9pjV917w3SCKQ22OorPzPaVtLu7v7BZ+zGSPth8ZF9cGMUHIBMxim/rtmU9qNWNtFdH+wAASItUAVXS2Hed3L1CUklaKgIAQKkDqmAr+wqbsxAAABpKFVD/NLMLNm80syFKrA8FAEBapBrFN0LSn8zsbH0VSOWS2ko6I52FAQCyW6oFCz+UdGT0Bd3k0u5/cfen0l4ZACCrbTWgzKxA0lAllnhfIGmau9e1RGEAgOyW6hnUDCVu6S2QdJKkW9JeEQAASv0Mqpe795YkM5sm6eX0lwQAQOoeVG3yBbf2AAAtKVUP6mAz+yJ6bZIKo22T5O7eIa3VAQCyVqpRfLlb2w8AQLqkusUHAEAsCCgAQJAIKABAkAgoAECQCCgAQJAIKABAkAgoAECQCCgAQJAIKABAkAgoAECQCCgAQJAIKABAkAgoAECQCCgAQJAIKABAkAgoAECQCCgAQJAIKABAkAgoAECQCCgAQJAIKABAkAgoAECQCCgAQJAIKABAkGIJKDObaGZvm9kbZvYnM+vYYN81ZrbUzBab2Ylx1AcAiF9cPai5kg5094MkLZF0jSSZWS9JgyR9Q1J/SXeaWW5MNQIAYhRLQLn739y9Ltp8UVL36PUASbPcvcbd35O0VNJhcdQIAIhXCM+gfizp8eh1N0nLGuxbHrV9jZldaGYVZlZRVVWV5hIBAC2tTbpObGbzJHVpZNdP3X12dMxPJdVJmpl8WyPHe2Pnd/epkqZKUnl5eaPHAAAyV9oCyt37bW2/mQ2WdIqk49w9GTDLJfVocFh3SR+kp0IAQMjiGsXXX9LVkk5z97UNdj0maZCZ5ZtZqaT9JL0cR40AgHilrQeVwh2S8iXNNTNJetHdh7r7m2b2B0lvKXHr7xJ33xBTjQCAGMUSUO6+71b23STpphYsBwAQoBBG8QEA8DUEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSLEGlJmNNDM3s87RtpnZr81sqZm9YWaHxlkfACA+sQWUmfWQdLyk9xs0nyRpv+jXhZKmxFAaACAAcfagbpM0SpI3aBsg6T5PeFFSRzPrGkt1AIBYxRJQZnaapBXu/vpmu7pJWtZge3nU1tg5LjSzCjOrqKqqSlOlAIC4tEnXic1snqQujez6qaRrJZ3Q2NsaafNG2uTuUyVNlaTy8vJGjwEAZK60BZS792us3cx6SyqV9LqZSVJ3Sa+Y2WFK9Jh6NDi8u6QP0lVjYyor39O9Uyapft3nyinoqPOGXamSktKWLAEAoBhu8bn7Anffzd1L3L1EiVA61N1XSXpM0g+j0XxHSPqPu69sqdoqK9/T5DFDdWVf19gzOuvKvq7JY4aqsvK9lioBABAJ7XtQcyS9K2mppN9Jurgl//B7p0zS2EH7qqgwT5JUVJinsYP21b1TJrVkGQAApfEW37aKelHJ1y7pkrhqqV/3uYoKO2/SVlSYp/p1H8dUEQBkr9B6ULHKKeioNdW1m7Stqa5VTkHHmCoCgOxFQDVw3rArNXbW0o0htaa6VmNnLdV5w66MuTIAyD6x3+ILSUlJqYaPu0uTpkxS/bqPlVPQUcPH3cUoPgCIgSUe+2S28vJyr6ioiLsMAGgSM5vv7uVx1xEqbvEBAIJEQAEAgkRAAQCCREABAIJEQAEAgpT1w8yZHBYAwpTVPSgmhwWAcGV1QDE5LACEK6sDKjE5bN4mbYnJYT+PqSIAQFJWBxSTwwJAuLI6oJgcFgDCldWj+JgcFgDCldU9KIaYA0C4sjagGGIOAGHL2oBiiDkAhC0rA6qy8j29++bLDDEHgIBlXUAlb+3tuYsYYg4AAcu6gEre2hvynQM19t4XGWIOAIHKumHmidkjOquoME+XnnGwJv1hvurrXQuXr9et0x5hFB8ABCLrAio5e0RRYZ726tJB1w8+QmuqazXpGSOcACAgWXeLj9kjACAzZF0PitkjACAzmLvHXcMOKy8v94qKirjLAIAmMbP57l4edx2hyrpbfACAzEBAAQCCREABAIJEQAEAgkRAAQCCREABAIJEQAEAgkRAAQCCREABAILUKmaSMLMqSf/ewu7Okj5uwXLSjesJW2u6ntZ0LVKY17OXuxfHXUSoWkVAbY2ZVbSmqUS4nrC1putpTdcitb7ryQbc4gMABImAAgAEKRsCamrcBTQzridsrel6WtO1SK3velq9Vv8MCgCQmbKhBwUAyEAEFAAgSK06oMxspJm5mXWOts3Mfm1mS83sDTM7NO4at4WZTTSzt6Oa/2RmHRvsuya6nsVmdmKcdTaFmfWPal5qZqPjrqepzKyHmT1tZovM7E0zuzxq72Rmc83snej3XeKutSnMLNfMXjWzP0fbpWb2UnQ9D5lZ27hr3FZm1tHMHo7+7iwys//O9M8n27TagDKzHpKOl/R+g+aTJO0X/bpQ0pQYStsecyUd6O4HSVoi6RpJMrNekgZJ+oak/pLuNLPc2KrcRlGNv1Hi8+gl6fvRtWSSOklXuvsBko6QdEl0DaMlPenu+0l6MtrOJJdLWtRg+2ZJt0XX85mkIbFUtX1+JekJd99f0sFKXFemfz5ZpdUGlKTbJI2S1HAUyABJ93nCi5I6mlnXWKprAnf/m7vXRZsvSuoevR4gaZa717j7e5KWSjosjhqb6DBJS939XXdfL2mWEteSMdx9pbu/Er1ercQPv25KXMeM6LAZkk6Pp8KmM7Pukr4j6e5o2yT9j6SHo0My5nrMrIOkb0maJknuvt7dP1cGfz7ZqFUGlJmdJmmFu7++2a5ukpY12F4etWWSH0t6PHqdqdeTqXU3ysxKJB0i6SVJu7v7SikRYpJ2i6+yJrtdiX/U1Ufbu0r6vME/jjLpc9pbUpWk6dEty7vNrEiZ/flknTZxF7C9zGyepC6N7PqppGslndDY2xppC2Kc/daux91nR8f8VIlbSzOTb2vk+CCuJ4VMrftrzGwnSX+UNMLdv0h0OjKPmZ0i6SN3n29mfZPNjRyaKZ9TG0mHShru7i+Z2a/E7byMk7EB5e79Gms3s96SSiW9Hv2w6C7pFTM7TIl/AfZocHh3SR+kudRtsqXrSTKzwZJOkXScf/XltWCvJ4VMrXsTZpanRDjNdPdHouYPzayru6+Mbh9/FF+FTXKUpNPM7GRJBZI6KNGj6mhmbaJeVCZ9TsslLXf3l6Lth5UIqEz9fLJSq7vF5+4L3H03dy9x9xIl/kc91N1XSXpM0g+j0XxHSPpPsrsfMjPrL+lqSae5+9oGux6TNMjM8s2sVInBHy/HUWMT/VPSftEIsbZKDPR4LOaamiR6PjNN0iJ3v7XBrsckDY5eD5Y0u6Vr2x7ufo27d4/+zgyS9JS7ny3paUnfiw7LpOtZJWmZmZVFTcdJeksZ+vlkq4ztQW2nOZJOVmIwwVpJP4q3nG12h6R8SXOjXuGL7j7U3d80sz8o8RevTtIl7r4hxjq3ibvXmdmlkv4qKVfSPe7+ZsxlNdVRks6VtMDMXovarpU0XtIfzGyIEiNIB8ZUX3O5WtIsM7tR0quKBh1kiOGSZkb/CHpXib/vOWpdn0+rxlRHAIAgtbpbfACA1oGAAgAEiYACAASJgAIABImAAgAEiYBCxjCzDWa+vDMRAAADNklEQVT2mpktNLP/NbN2UXsXM5tlZv8ys7fMbI6Z9Yz2PWFmnydn597KuW83s29Fr2dGM60vNLN7oi/kyszOjmaUf8PM/mFmB2/hXNPM7PXouIej2SZkZsOjc85JzgpuZkeb2a0N3ltsZk80x38vINMRUMgk1e7+TXc/UNJ6SUOjL8z+SdIz7r6Pu/dS4vtIu0fvmajE95W2yMw6STrC3f8eNc2UtL+k3pIKJZ0ftb8n6dvRrPI3aMtLiF/h7gdHx70v6dKo/XxJBynxfaITo9p/Fp1LkuTuVZJWmtlRqf9zAK0bAYVM9ZykfSUdK6nW3e9K7nD319z9uej1k5JWpzjX9yRt7LW4+5xoxntXYmaO7lH7P9z9s+iwhrPKb8Ldv5A2zjZRqE3nr8uT1E5SrRLBOafBOZMelXR2ipqBVo+AQsYxszZKrCW1QNKBkubv4CmPauwc0a29c9UgvBoYoq9mlW+sxumSVinRE5scNd+iRLAVS3pBial27mzk7RWSjtn28oHWiYBCJimMphWqUOLWWXNNu9NViaUZNnenpL8ne2NJZnasEgF19ZZO6O4/krSHEutE/b+o7X53P8Tdz5H0E0m/lnRS9JzqNjNL/n38KHovkNUIKGSS5DOob7r78Gixwzcl9dnR8yoxg/dGZjZGiZ7OTzZrP0iJBf0GuPsnWztpNC/iQ5K+u9k59pD0X9EyKtcpEWA1SkxoqqiW6u29GKC1IKCQ6Z6SlG9mFyQbzOy/zOzbTTjHIiWeZyXff76kEyV9393rG7TvKekRSee6+5LGThTNlL9v8rWkUyW9vdlhNygxOEL66hlVvRLPpiSpp6SFTagfaJUIKGS0aCDDGZKOj4aZvylprKJ1i8zsOUn/K+k4M1tuZic2cpq/SOrbYPsuJUYB/v9oWPv1Ufv1Sqwye2fUXpF8QzR0fA8lFvmbYWYLlHhG1lXSzxscd0hU96tR07TouEP11bOuY6OagKzGbOaAJDN7XtIp7v55ALX8XYlbiJuP7gOyCgEFSDKzw5V4xvVGzHUUSzrK3R+Nsw4gBAQUACBIPIMCAASJgAIABImAAgAEiYACAASJgAIABOn/ADOE6hYbyUctAAAAAElFTkSuQmCC\n",
      "text/plain": [
       "<Figure size 432x360 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "## Load your assembly\n",
    "data = ipyrad.load_json(\"/tmp/ipyrad-test/rad.json\")\n",
    "## Create they pca object\n",
    "pca = ipa.pca(data)\n",
    "## Bam!\n",
    "pca.plot()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Full guide\n",
    "In the most common use you'll want to plot the first two PCs, then inspect the output, remove any obvious outliers, and then redo the pca. It's often desirable to import a vcf file directly rather than to use the ipyrad assembly, so here we'll demonstrate this."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "metadata": {},
   "outputs": [],
   "source": [
    "## Path to the input vcf, in this case it's just the vcf from our ipyrad pedicularis assembly\n",
    "vcffile = \"/home/isaac/ipyrad/test-data/pedicularis/ped_outfiles/ped.vcf\""
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Here we can just load the vcf file directly into the pca analysis module. Then ask for the samples in `samples_vcforder`, which is the order in which they are written in the vcf."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "[u'29154_superba_SRR1754715' u'30556_thamno_SRR1754720'\n",
      " u'30686_cyathophylla_SRR1754730' u'32082_przewalskii_SRR1754729'\n",
      " u'33413_thamno_SRR1754728' u'33588_przewalskii_SRR1754727'\n",
      " u'35236_rex_SRR1754731' u'35855_rex_SRR1754726' u'38362_rex_SRR1754725'\n",
      " u'39618_rex_SRR1754723' u'40578_rex_SRR1754724'\n",
      " u'41478_cyathophylloides_SRR1754722' u'41954_cyathophylloides_SRR1754721']\n"
     ]
    }
   ],
   "source": [
    "pca = ipa.pca(vcffile)\n",
    "print(pca.samples_vcforder)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Now construct the default plot, which shows all samples and PCs 1 and 2. By default all samples are assigned to one population, so everything will be the same color."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Using default cmap: Spectral\n"
     ]
    },
    {
     "data": {
      "text/plain": [
       "<matplotlib.axes._subplots.AxesSubplot at 0x7fe0beb3a650>"
      ]
     },
     "execution_count": 5,
     "metadata": {},
     "output_type": "execute_result"
    },
    {
     "data": {
      "image/png": "iVBORw0KGgoAAAANSUhEUgAAAagAAAFgCAYAAADuCe0ZAAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAAALEgAACxIB0t1+/AAAADl0RVh0U29mdHdhcmUAbWF0cGxvdGxpYiB2ZXJzaW9uIDIuMi4yLCBodHRwOi8vbWF0cGxvdGxpYi5vcmcvhp/UCwAAIABJREFUeJzt3Xt8VeWd7/HvLwlsICjXKBAYs6WJEPECBMS2p0dKqeJYseq0Cra246WdGQeFSlVaTqelzauW8dLM6ehYrbUYvFSrUqV4KLVVzwzaUKtFkZuJEi4SKAJNmQjhN3/slbjBkE1IdtYD+bxfr/3K2s9aWevHAvLNetazn2XuLgAAQpMTdwEAALSEgAIABImAAgAEiYACAASJgAIABImAAgAEiYACAASJgAIABImAAgAEKS/uArJp4MCBXlRUFHcZAAKyYsWKbe5eEHcdyOyYDqiioiJVVVXFXQaAgJjZ23HXgMNDFx8AIEgEFAAgSAQUACBIx/Q9KAA4HCtWrDghLy/vXkmjxC/unWm/pJX79u27euzYsVsPXklAAejy8vLy7h00aNDIgoKCHTk5OTwkr5Ps37/f6urqSrds2XKvpAsPXs9vCgAgjSooKNhFOHWunJwcLygo2KnUleuH13dyPQAQohzCKR7ReW8xiwgoAECQCCgAaKPVq97s/tWplyevPOu8kq9OvTy5etWb3eOu6VhEQEmqqa7WrCuu0VcmflazrrhGNdXVcZcEIFCrV73Z/aZzp5ecuWhz/3N+v/e4Mxdt7n/TudNLOiqkfvazn/U1s7GvvPJKD0lavXp19+Li4lMl6emnnz5u4sSJH+mI47RXel3Z0uUDqqa6WrMnT1NJ5VpN+N0ulVSu1ezJ0wgpAC264+ZvF5634bhEwnIlSQnL1XkbjkvccfO3Czti/w8//HD/MWPG/GXBggX9O2J/R7MuH1AVc8s1aX1C6f/YJq1PqGJuecyVAQjRni07ujX9vGiSsFzteXdHt/bue+fOnTlVVVW977///ponnniiX1u//5lnnuk9YsSI0hEjRpSOHDmydMeOHTk7d+7MOfvss0tKS0tHlpSUlD744IN9pdQVUDKZPPXzn//8ScXFxadeeOGFySeffPK4MWPGjDjppJNGPffcc70kadasWUMuuuii5IQJE0pOOumkUbfddtvAg4+7b98+feUrXxk6atSokSUlJaXz588fKElvv/12t7KyslNGjBhRWlxcfOqSJUt6t+XP0+U/B1W/cZta+sdWv2l7TBUBCFnPQf32NvjmA35uNHijep54wt727ruysrLvOeecs/P0009v6Nu3b+OLL77Yq6CgYN/hfv9tt902qKKi4u1Pf/rT9Tt37szp1avXfkl65pln1vXv33//5s2b884666wR06ZNe0+SNmzY0OORRx55a+zYsW+ffvrpIysrKwdUVVW9uXDhwr7f+973Bk+cOHG9JK1atarnihUrVu3evTt39OjRpZdccsnO9OPeeeedA/v06dO4cuXKVXv27LFx48aN+MxnPrProYce6jdp0qSdt95665Z9+/Zp9+7dbboo6vJXUPmFA9XgjQe0NXij8ocMiKkiACGb+f1vbVwybHdD08+NBm/UkmG7G2Z+/1sb27vvRx99tP/ll1++Q5IuueSSP7e1m2/ChAl/ufHGG4d997vfPWHbtm253bp10/79++2GG24YWlJSUjpx4sSSrVu3dq+trc2TpMLCwobx48fvyc3NVUlJyZ5PfvKTu3JycjRmzJi/1tbWJpr2O2XKlPd69+7tgwcP3nf22WfveuGFF/LTj/vrX//6+EcffXTAiBEjSkePHj1yx44deW+88UaPCRMm1D/00EMDZ82aNeTll1/u2a9fv/1t+fN0+SuoGfPmaPbyac3dfA3eqGXDGzR/3py4SwMQoFNGjnj/1mcr19xx87cL97y7o1vPE0/Ye+v37954ysgR77dnv1u2bMldvnz58WvWrOl53XXXqbGx0czMZ86c+aEpgA6lvLx8y0UXXbTzqaee6vPRj3505JIlS9a88MIL+du3b8/705/+tCqRSHhhYeFpe/bsyZGk7t27N3/2KycnRz169HBJys3NVWNjozWtM7MDjnPwe3e322677Z1LLrlk18E1Pf/886sff/zxPl/60peSM2bMePe666477O6pLh9QRcmk5i9dqIq55arftF35QwZo/rw5Kkom4y4NQKBOGTni/bufeqhDR1ItWLCg38UXX7x94cKFzc+rGjdu3Ck1NTWHPTrw9ddfT4wfP37P+PHj97z00kv5K1eu7LFz587cgQMH7k0kEv7LX/7yuE2bNrV5tOGvfvWrvt/73vc279q1K2f58uXH3XHHHRsbGhqaU2ry5Mk777rrroILLrhgdyKR8Ndeey1RVFS0d8uWLXnJZPL9r33ta9vq6+tz/vCHP/SSREC1RVEyqdsf/HHcZQDown7+858P+PrXv745vW3q1Kk7ysvLBx/uPn7wgx+c8J//+Z/H5+TkeElJyZ5LL71053vvvZc7ZcqUj4waNWrkqaee+tdkMvnfba1t9OjR9ZMmTSretGlT9xtvvHFzUVHR3tWrVzcH3cyZM7fV1NQkTjvttJHubv3799+7ePHi9c8+++xxFRUVg/Ly8rxXr16NlZWVbQp1cz92Z/coKytznqgLIJ2ZrXD3svS2V199teaMM87YFldNIZs1a9aQ3r17N37nO995N1vHePXVVweeccYZRQe3d/lBEgCAMNHFBwBHmR/+8IcD7rrrrhPT28aNG/eXBQsWvNPRx7r99ts3dfQ+DxcBBQDS/v3799vRMqP59ddfv/36668/Jj6suX//flPqwYUfQhcfAEgr6+rq+kQ/LNFJogcW9pG0sqX1XEEB6PL27dt39ZYtW+7dsmULj3zvXM2PfG9pJQEFoMsbO3bsVrXwyHHEi98UAABBIqAAAEEioAAAQSKgAABBIqAAAEEioAAAQYo1oMysr5k9ZmZvmtkqMzvbzPqb2VIzWxt97Rdta2ZWYWbrzOw1MxsTZ+0AgOyK+wrqh5KWuPsISWdIWiXpZknL3L1Y0rLovSRNkVQcva6VdFe2iqqprtasK67RVyZ+VrOuuEY11R362BcAwGGILaDM7HhJn5B0nyS5+/vu/p6kqZIeiDZ7QNJF0fJUST/zlOWS+prZYT8n5XDVVFdr9uRpKqlcqwm/26WSyrWaPXkaIQUAnSzOK6iTJdVJut/MXjGze80sX9KJ7r5ZkqKvJ0TbF0rakPb9tVHbAczsWjOrMrOqurq6NhdVMbe8+fHvkpSwXE1an1DF3PI27wsAcOTiDKg8SWMk3eXuoyXV64PuvJa0NInjh2Yedvd73L3M3csKCgraXFT9xm3N4dQkYbmq33RMTBwMAEeNOAOqVlKtu78UvX9MqcB6t6nrLvq6NW37YWnfP1RShz+nJL9woBq88YC2Bm9U/pABHX0oAEArYgsod98iaYOZnRI1TZL0hqRFkq6M2q6U9FS0vEjSF6PRfBMk7WzqCuxIM+bN0bLhDc0h1eCNWja8QTPmzenoQwEAWhH3bOb/LKnSzLpLekvSl5UKzUfN7CpJ70j6u2jbxZLOl7RO0l+jbTtcUTKp+UsXqmJuueo3bVf+kAGaP2+OipLJbBwOAHAI5n5UPEDyiJSVlXlVVVXcZQAIiJmtcPeyuOtAZnF/DgoAgBYRUACAIBFQAIAgEVAAgCARUACAIBFQAIAgEVAAgCARUACAIBFQAIAgEVAAgCARUACAIBFQAIAgEVAAgCARUACAIBFQAIAgEVAAgCARUACAIBFQAIAgEVAAgCARUACAIBFQAIAgEVAAgCARUACAIBFQAIAgEVAAgCARUACAIBFQAIAgEVAAgCARUACAIBFQAIAgEVAAgCARUACAIBFQAIAgEVAAgCARUACAIBFQAIAgEVAAgCARUACAIBFQAIAgEVAAgCARUACAIBFQAIAgEVAAgCARUACAIBFQAIAgEVAAgCARUACAIBFQAIAgxR5QZpZrZq+Y2dPR+6SZvWRma83sETPrHrUnovfrovVFcdYNAMiu2ANK0vWSVqW9v1XSHe5eLGmHpKui9qsk7XD3j0i6I9oOAHCMijWgzGyopL+VdG/03iR9UtJj0SYPSLooWp4avVe0flK0PQDgGBT3FdSdkr4uaX/0foCk99x9X/S+VlJhtFwoaYMkRet3RtsfwMyuNbMqM6uqq6vLZu0AgCyKLaDM7AJJW919RXpzC5v6Yaz7oMH9Hncvc/eygoKCDqgUABCHvBiP/TFJF5rZ+ZJ6SDpeqSuqvmaWF10lDZW0Kdq+VtIwSbVmliepj6Q/d37ZAIDOENsVlLvf4u5D3b1I0mWSfuPu0yU9J+nSaLMrJT0VLS+K3ita/xt3/9AVFADg2BD3PaiW3CRplpmtU+oe031R+32SBkTtsyTdHFN9AIBOEGcXXzN3/62k30bLb0ka38I2/y3p7zq1MABAbEK8ggIAgIACAISJgAIABImAAgAEiYACAASJgAIABImAAgAEiYACAASJgAIABImAAgAEiYACAASJgAIABImAAgAEiYACAASJgAIABImAAgAEiYACAASJgAIABImAAgAEiYACAASJgAIABImAAgAEiYACAASJgAIABImAAgAEqU0BZWb5ZpabrWIAAGjSakCZWY6ZTTOzZ8xsq6Q3JW02s9fNbL6ZFXdOmQCAribTFdRzkoZLukXSIHcf5u4nSPpfkpZL+r6ZXZHlGgEAXVBehvWfcve9Bze6+58lPS7pcTPrlpXKAABdWqsBdXA4mVkPSVdI6ilpobtvbynAAABor7aO4vuhpFxJ/y3pyY4vBwCAlEyDJBaa2fC0pv6SKiU9JKlfNgsDAHRtme5BfVPSd81sk6R5kv5V0iJJPST9S3ZLAwB0ZZnuQb0laZqZfVzSI5KekTTZ3Rs7ozgAQNeVqYuvn5n9k6RSSZ+TtFPSs2Z2QWcUBwDoujINknhSUoNSXXoL3P1nkj4jaayZLcp2cQCArivTPagBkhYqNaz8i5Lk7nskfdvMBme5NgBAF5YpoL4laamkRkk3p69w983ZKgoAgEyDJB5XasYIAAA6VaZBEteZ2cBoebiZPW9m75nZS2Z2WueUCADoijINkvgHd98WLVdIusPd+0q6SdLdWa0MANClZQqo9C7AE9z9CUly999KOi5bRQEAkCmgHjOzn5rZyZKeMLMbzOxvzOzLkt7phPoAAF1UpkES3zCzLyk1995wSQlJ1yr1+ajpWa8OANBlZRpmLnf/qaSfZr0SAADStPVxG83MbFBHFgIAQLojDihJ93VYFQAAHOSIA8rd/7YjCwEAIF2bA8rM+mejEAAA0mWaSeKbaculZrZG0gozqzGzs9pzYDMbZmbPmdkqM3vdzK6P2vub2VIzWxt97Re1m5lVmNk6M3vNzMa05/gAgLBluoK6OG15vqTr3T2p1LOh7mjnsfdJ+pq7j5Q0QdI/mVmpUpPSLnP3YknL9MEktVMkFUevayXd1c7jd4ia6mrNuuIafWXiZzXrimtUU10dd0kAcExoSxffEHf/lSS5+8tKPYLjiLn7Znf/Q7S8W9IqSYWSpkp6INrsAUkXRctTJf3MU5ZL6hv3Iz9qqqs1e/I0lVSu1YTf7VJJ5VrNnjyNkAKADpApoE42s0Vm9ktJQ82sV9q6bh1VhJkVSRot6SVJJzY9yiP6ekK0WaGkDWnfVhu1Hbyva82sysyq6urqOqrEFlXMLdek9QklLFeSlLBcTVqfUMXc8qweFwC6gkwf1J160PscSTKzE9VBXWxm1lupR3rc4O67zOyQm7bQ5h9qcL9H0j2SVFZW9qH1Hal+47bmcGqSsFzVb9qezcMCQJeQaaqj3x2i/V1JP2rvwc2sm1LhVOnuv4ia3zWzwe6+OerC2xq110oalvbtQyVtam8N7ZFfOFANvuOAkGrwRuUPGRBjVQBwbGjPTBL3tOfAlrpUuk/SKne/PW3VIklXRstXSnoqrf2L0Wi+CZJ2xv1U3xnz5mjZ8AY1eKOkVDgtG96gGfPmxFkWABwTWr2CauUzTybp/HYe+2OSviDpT2b2x6htjqTvS3rUzK5Sasb0v4vWLY6OuU7SXyV9uZ3Hb7eiZFLzly5Uxdxy1W/arvwhAzR/3hwVJZNxlwYARz1zP/RtGjNrlPS2Drz/49H7Qnfvnt3y2qesrMyrqqriLgNAQMxshbuXxV0HMss0SOItSZPc/UPPfjKzDS1sDwBAh8h0D+pOSf0Ose4HHVwLAADNMo3iO+RIPXf/t44vBwCAlExz8X08w/rjzWxUx5YEAEDme1CXmNkPJC2RtEJSnaQekj4iaaKkkyR9LasVAgC6pExdfDOj2cQvVWq492BJe5SaN+8/3P3F7JcIAOiKMl1Byd13SPpx9AIAoFO055HvAABkDQEFAAgSAQUACFLGgIqGkg9vof307JQEAEDmz0F9TtKbkh43s9fNbFza6p9mszAAQNeW6QpqjqSx7n6mUrOHLzCzi6N1h3yyIAAA7ZVpmHlu2uPXXzaziZKeNrOhauFptgAAdJRMV1C70+8/RWF1jlKPgj81i3UBALq4TFdQ/6CDuvLcfbeZnSfpc1mrCgDQ5WW6gqqXdGIL7RMkLe/4cgAASDmc50HtbqF9T7QOAICsyBRQRe7+2sGN7l4lqSgrFQEAoMwB1aOVdT07shAAANJlCqjfm9k1Bzea2VVKPR8KAICsyDSK7wZJT5jZdH0QSGWSukv6bDYLA4BMaqqrVTG3XPUbtym/cKBmzJujomQy7rLQQTI9sPBdSR+NPqDb9Gj3Z9z9N1mvDABaUVNdrdmTp2nS+oQSlqsG36HZy6dp/tKFhNQxItNcfD3M7AZJl0h6X9JdhBOAEFTMLW8OJ0lKWK4mrU+oYm55zJWho2S6B/WAUl16f5I0RdK/Zr0iADgM9Ru3NYdTk4Tlqn7T9pgqQkfLdA+q1N1PkyQzu0/Sy9kv6dhCHzmQHfmFA9XgOw4IqQZvVP6QATFWhY6U6Qpqb9OCu+/Lci3HnKY+8pLKtZrwu10qqVyr2ZOnqaa6Ou7SgKPejHlztGx4gxq8UVIqnJYNb9CMeXNirgwdJVNAnWFmu6LXbkmnNy2b2a7OKPBoRh85kD1FyaTmL12oNdOLtfycPlozvZgBEseYTKP4cltbj9bRRw5kV1Eyqdsf/HHcZSBLMj7yHUcu1UfeeEAbfeQAcHgIqCyijxwAjhwBlUX0kQPAkcs0zBztRB85ABwZrqAAAEHiCgpAq/iwOeJCQAE4JCZkRZzo4gNwSHzYHHEioAAcEh82R5wIKACHxIfNEScCCsAh8WFzxImAAnBIfNgccWIUH4BW8WFzxIUrKABAkAgoAECQCCgAQJAIKABAkAgoAECQCCgAQJCOuoAys/PMbLWZrTOzm+OuBwCQHUdVQJlZrqQfSZoiqVTS5WZWGm9VAIBsOKoCStJ4Sevc/S13f1/Sw5KmxlwTACALjraAKpS0Ie19bdTWzMyuNbMqM6uqq6vr1OIAAB3naAsoa6HND3jjfo+7l7l7WUFBQSeVBQDoaEdbQNVKGpb2fqikTTHVAgDIoqMtoH4vqdjMkmbWXdJlkhbFXBMAIAuOqtnM3X2fmV0n6VlJuZJ+4u6vx1wWACALjqqAkiR3Xyxpcdx1AACy62jr4gMAdBEEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSHlxFwCpprpaFXPLVb9xm/ILB2rGvDkqSibjLgsAYkVAxaymulqzJ0/TpPUJJSxXDb5Ds5dP0/ylCwkpAF0aXXwxq5hb3hxOkpSwXE1an1DF3PKYKwOAeBFQMavfuK05nJokLFf1m7bHVBEAhIGAill+4UA1eOMBbQ3eqPwhA2KqCADCQEDFbMa8OVo2vKE5pBq8UcuGN2jGvDkxVwYA8SKgYlaUTGr+0oVaM71Yy8/pozXTixkgAQBiFF8QipJJ3f7gj+MuAwCCwhUUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSAQUACBIBBQAIEgEFAAgSLEElJnNN7M3zew1M3vCzPqmrbvFzNaZ2WozOzet/byobZ2Z3RxH3QCAzhPXFdRSSaPc/XRJayTdIklmVirpMkmnSjpP0r+bWa6Z5Ur6kaQpkkolXR5tCwA4RsUSUO7+/9x9X/R2uaSh0fJUSQ+7e4O7V0taJ2l89Frn7m+5+/uSHo62BQAco0K4B/X3kn4VLRdK2pC2rjZqO1T7h5jZtWZWZWZVdXV1WSgXANAZsvbAQjP7taRBLaz6hrs/FW3zDUn7JFU2fVsL27taDlJv6bjufo+keySprKysxW0AAOHLWkC5+6daW29mV0q6QNIkd28KklpJw9I2GyppU7R8qHYAwDEorlF850m6SdKF7v7XtFWLJF1mZgkzS0oqlvSypN9LKjazpJl1V2ogxaLOrhsA0HmydgWVwf+VlJC01Mwkabm7f9XdXzezRyW9oVTX3z+5e6Mkmdl1kp6VlCvpJ+7+ejylAwA6g33Qu3bsKSsr86qqqrjLABAQM1vh7mVx14HMQhjFBwDAhxBQAIAgxXUPqsupqa5Wxdxy1W/cpvzCgZoxb46Kksm4ywKAYBFQnaCmulqzJ0/TpPUJJSxXDb5Ds5dP0/ylCwkpADgEuvg6QcXc8uZwkqSE5WrS+oQq5pbHXBkAhIuA6gTb19c2h1OThOWqftP2mCoCgPARUFn24vPP65WqFWpIfZyrWYM3Kn/IgJiqAoDwcQ8qi2qqq/WP512mQft66KdapQLvqU9oiI5Tdz3ae4N+Mu/+uEsEgGARUFn03Zk36/g9rot1cjQ4olGVWq2eytOgUaUMkACAVtDFl0Vv/tcrukKnHDA4YrpOUUK5GjT8b2KuDgDCxhVUFiUsr8XBEe9ag/5t3pyYqgKAowNXUFk05LTiFgdHnPLJ8XTvAUAGBFSW1FRXa/vqDXpM65tDqsEb9Uzhe5p3z50xVwcA4aOLL0sq5pbrwg3Ha7d6aInekbtrv1wnlp3F1RMAHAYCKkvqN25TwnKVUE9N1QeBtHzX3hirAoCjB118WZJfOJAP5wJAO3AFlSUXX/sFzXzqCg36S45ylaNxOkF/HG6az+g9ADgsBFQW1FRX64d/f5Ou+ssHH9B9tPcGfecn93P/CQAOE118WdDS7OWf+8sw/eKeBTFXBgBHDwIqC5oGSKRj9nIAaBsCKgsYIAEA7UdAZcGMeXO0bHjDAR/QXTa8QTMYIAEAh42AyoKiZFLzly7UmunFWn5OH62ZXszj3QGgjRjFlyVFyaRuf/DHcZcBAEctrqAAAEEioAAAQSKgAABBIqAAAEEioAAAQSKgAABBIqAAAEEioAAAQSKgAABBMnePu4asMbM6SW9nafcDJW3L0r7bi9qOXMj1UduRS6/vJHcviLMYHJ5jOqCyycyq3L0s7jpaQm1HLuT6qO3IhV4fWkYXHwAgSAQUACBIBNSRuyfuAlpBbUcu5Pqo7ciFXh9awD0oAECQuIICAASJgAIABImAysDM5pvZm2b2mpk9YWZ909bdYmbrzGy1mZ2b1n5e1LbOzG7u5HpjO3Z0/GFm9pyZrTKz183s+qi9v5ktNbO10dd+UbuZWUVU72tmNqYTasw1s1fM7OnofdLMXopqe8TMukftiej9umh9UZbr6mtmj0X/3laZ2dmBnbeZ0d/pSjN7yMx6xHXuzOwnZrbVzFamtbX5XJnZldH2a83syo6sER3A3Xm18pL0aUl50fKtkm6NlkslvSopISkpab2k3Oi1XtLJkrpH25R2Uq2xHTuthsGSxkTLx0laE52rH0i6OWq/Oe08ni/pV5JM0gRJL3VCjbMkLZT0dPT+UUmXRct3S/qHaPkfJd0dLV8m6ZEs1/WApKuj5e6S+oZy3iQVSqqW1DPtnH0prnMn6ROSxkhamdbWpnMlqb+kt6Kv/aLlftn+98erDX/PcRdwNL0kfVZSZbR8i6Rb0tY9K+ns6PVsWvsB22W5vtiO3UpNT0maLGm1pMFR22BJq6Pl/5B0edr2zdtlqZ6hkpZJ+qSkp6MfWtv0wS8hzeew6e80Ws6LtrMs1XV8FAB2UHso561Q0oboh3ledO7OjfPcSSo6KKDadK4kXS7pP9LaD9iOV/wvuvja5u+V+k1M+uA/bJPaqO1Q7Z0hzmN/SNStM1rSS5JOdPfNkhR9PSHarLNrvlPS1yXtj94PkPSeu+9r4fjNtUXrd0bbZ8PJkuok3R91P95rZvkK5Ly5+0ZJ/yrpHUmblToXKxTGuWvS1nMV1P8XfBgBJcnMfh31qx/8mpq2zTck7ZNU2dTUwq68lfbOEOexD2BmvSU9LukGd9/V2qYttGWlZjO7QNJWd19xmMfvzPOZp1SX1V3uPlpSvVLdVIfSqX/X0f2cqUp1Zw+RlC9pSis1BPNvUWH+X8VhyIu7gBC4+6daWx/dPL1A0iSP+gKU+m1rWNpmQyVtipYP1Z5trdXUacysm1LhVOnuv4ia3zWzwe6+2cwGS9oatXdmzR+TdKGZnS+ph1LdandK6mtmedFv+unHb6qt1szyJPWR9Ocs1VYrqdbdX4reP6ZUQIVw3iTpU5Kq3b1OkszsF5I+qjDOXZO2nqtaSecc1P7bLNeINuAKKgMzO0/STZIudPe/pq1aJOmyaLRSUlKxpJcl/V5ScTS6qbtSN4gXdVK5cR5bUmrElKT7JK1y99vTVi2S1DRK6kql7k01tX8xGmk1QdLOpm6ajubut7j7UHcvUurc/Mbdp0t6TtKlh6itqeZLo+2z8hu2u2+RtMHMTomaJkl6QwGct8g7kiaYWa/o77ipvtjPXZq2nqtnJX3azPpFV4ifjtoQirhvgoX+krROqX7qP0avu9PWfUOpUXOrJU1Jaz9fqdFr6yV9o5Prje2jaMDIAAADvUlEQVTY0fE/rlQ3yWtp5+x8pe4/LJO0NvraP9reJP0oqvdPkso6qc5z9MEovpOV+uVinaSfS0pE7T2i9+ui9SdnuaYzJVVF5+5JpUaWBXPeJH1b0puSVkpaoNQI1ljOnaSHlLoXtlepK6GrjuRcKXVfeV30+nJn/3/h1fqLqY4AAEGiiw8AECQCCgAQJAIKABAkAgoAECQCCgAQJAIKnc7MGs3sj9FsHT83s15R+yAze9jM1pvZG2a22MxKonVLzOw9i2Ygb2Xfd5rZJ6LlSkvN7L4ymv26W9Q+wsz+y8wazOzGVvY1ycz+ENX6opl9JGr/52ifi9Nm7/64md2e9r0FZrakvecK6MoIKMRhj7uf6e6jJL0v6avRhz+fkPRbdx/u7qWS5kg6Mfqe+ZK+0NpOzay/pAnu/nzUVClphKTTJPWUdHXU/mdJM5SaW641d0ma7u5nKjX7+Tej9qslnS7pFUnnRrXPlTSv6Rs9NePCZjP7WIZjADgEAgpxe0HSRyRNlLTX3e9uWuHuf3T3F6LlZZJ2Z9jXpZKar1rcfbFHlPqw6NCofau7/16pD3m2xpWaDklKTdWTPpVQN0m9on18QdJid99x0Pc/KWl6hmMAOATm4kNsojnapigVKqOUmh27PT6m1Bx2Bx+nm1Ihcn0b93e1pMVmtkfSLqWeJSSlrryWS3pd0v9XKojOa+H7qyR9t43HBBDhCgpx6Glmf1TqB/g7Ss3d1xEGK/XIioP9u6Tnm67G2mCmpPPdfaik+yXdLknuvsDdR7v7FUo9/LBC0hRLPQ33DjNr+n+1VamZvwEcAQIKcWi6B3Wmu/+zu7+v1NXI2PbuV6k54JqZ2bckFSgVJIfNzAokneEfzC7+iFKzd6dvM0TSOHd/Sqn7U5+X1KDURKqKatnTxj8DgAgBhVD8RlLCzK5pajCzcWb2v9uwj1VK3c9q+v6rlXrq6+Xuvv+Q39WyHZL6NI0iVOqpwKsO2maeUoMjpNQgDFfqQYi9orYSpSZWBXAECCgEIRrI8FlJk6Nh5q9L+hdFAxPM7AWlZseeZGa1ZnZuC7t5Rgc+3+dupUYB/lc0VPz/RPsaZGa1Sl1VfTPa3/HRusVmNsRTzze6RtLjZvaqUvewZjft2MxGR3W/EjXdp9RM2WP0wUCNiVFNAI4As5njmGJmL0q6wN3fC6CW5yVNbWF0H4DDQEDhmGJmZyl1j+u1mOsokPQxd38yzjqAoxkBBQAIEvegAABBIqAAAEEioAAAQSKgAABBIqAAAEH6HzC1kd2ohtOsAAAAAElFTkSuQmCC\n",
      "text/plain": [
       "<Figure size 432x360 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "pca.plot()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Population assignment for sample colors\n",
    "In the tl;dr example the assembly of our simulated data had included a `pop_assign_file` so the pca() was smart enough to find this and color samples accordingly. In some cases you might not have used a pops file, so it's also possible to specify population assignments in a dictionary. The format of the dictionary should have populations as keys and lists of samples as values. Sample names need to be identical to the names in the vcf file, which we can verify with the `samples_vcforder` property of the pca object."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "metadata": {},
   "outputs": [],
   "source": [
    "pops_dict = {\n",
    "    \"superba\":[\"29154_superba_SRR1754715\"],\n",
    "    \"thamno\":[\"30556_thamno_SRR1754720\", \"33413_thamno_SRR1754728\"],\n",
    "    \"cyathophylla\":[\"30686_cyathophylla_SRR1754730\"],\n",
    "    \"przewalskii\":[\"32082_przewalskii_SRR1754729\", \"33588_przewalskii_SRR1754727\"],\n",
    "    \"rex\":[\"35236_rex_SRR1754731\", \"35855_rex_SRR1754726\", \"38362_rex_SRR1754725\",\\\n",
    "            \"39618_rex_SRR1754723\", \"40578_rex_SRR1754724\"],\n",
    "    \"cyathophylloides\":[\"41478_cyathophylloides_SRR1754722\", \"41954_cyathophylloides_SRR1754721\"]\n",
    "}"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Using default cmap: Spectral\n"
     ]
    },
    {
     "data": {
      "text/plain": [
       "<matplotlib.axes._subplots.AxesSubplot at 0x7fe092fbbe50>"
      ]
     },
     "execution_count": 7,
     "metadata": {},
     "output_type": "execute_result"
    },
    {
     "data": {
      "image/png": "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q68OrFu3LmDJkiWtfvfGBoIBRUSNitFoLH5r/pKspW+szHxr/pIse/Xe69y5c9H7778fYDAYTOfOnXObOnXq6fLLt2/f7pucnOw7e/bsdraOEtnZ2e7vvPPOkZ9++snXYDCYunTpErlw4cI2tvdce+21+aGhoUUAMGDAgLxTp065Dxo0KA8Ann766VPr1q0L6NGjR3hmZqaXt7d3Ga6yatWqVgaDITI8PNy0b98+r0cfffTyeH3NmjUr+/rrr/cvXLiw7erVq1vY4xjUN063QUS6odfpNjIyMjxuvfXWsH379qU4s47GqrLpNngGRUREusSAIiKqhtFoLObZU/1jQBERkS4xoIiISJcYUEREpEsMKCIi0iUOdUREjUp6RobH7DffDD5vNru38PQsmf7UU0fDHTySub35+Pj0LCgo+NnZdTgbz6CIqNFIz8jwiH1huiG/X99W3rfe4p/fr2+r2BemG9IbyJxQZWVlqGjEiKaKAUVEjcbsN98MDrrrLk83Ly8AgJuXF4Luustz9ptv1mm6jYsXL7oMGDCgq9FoNIWFhUW+++67LYODg7sdP37cDbBOBti7d28jAEyePLndnXfeGXr99dcbOnXqFDV//vzWtu288MILbaOioiIMBoNp0qRJ7QDrl4A7d+4cOXr06I6RkZGmAwcOeADAX//61/YmkymiT58+hmPHjrkBwPz581tHRUVFGI1G09ChQ7vk5eU16s/wRv3LEVHTct5sdreFk42blxcumM11mm7js88+axYYGFiSkZGRum/fvpS77777YlXrp6WleW/fvn3f7t270+fOndsuOzvb/bPPPmu2f/9+rz179qSlpaWl/vLLLz5bt271A6xTeYwbN+5sWlpaqsFgKC4sLHTp1atXQWpqatqNN96Yl5CQ0A4ARo0adS45OTktIyMj1Wg0Fi5YsKB1VXU0dE4PKBFxFZGfRWSz9jpURP4jIvtE5EMR8dDaPbXX+7XlIc6sm4j0p4WnZ4ml6MrZJSxFRWju6Vmn6TZ69epVuGvXrmYTJkwI/uqrr/wCAgKqvA43fPjw835+fiooKMjSp0+fi7t27fL96quvmn333XfNTCaTSTtT8kpPT/cCgKCgoOKBAwfm297v4uKCRx55JBcAHn744bP//e9//QAgKSnJOzo62mgwGEyffvppQEpKilfFFTQOTg8oAE8BSCv3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\n",
      "text/plain": [
       "<Figure size 432x360 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "pca = ipa.pca(vcffile, pops_dict)\n",
    "pca.plot()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "This is just much nicer looking now, and it's also much more straightforward to interpret."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Removing \"bad\" samples and replotting.\n",
    "In PC analysis, it's common for \"bad\" samples to dominate several of the first PCs, and thus \"pop out\" in a degenerate looking way. Bad samples of this kind can often be attributed to poor sequence quality or sample misidentifcation. Samples with lots of missing data tend to pop way out on their own, causing distortion in the signal in the PCs. Normally it's best to evaluate the quality of the sample, and if it can be seen to be of poor quality, to remove it and replot the PCA. The Pedicularis dataset is actually very nice, and clean, but for the sake of demonstration lets imagine the cyathophylloides samples are \"bad samples\".\n",
    "\n",
    "We can see that the cyathophylloides samples have particularly high values of PC2, so we can target them for removal in this way."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 13,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "[False False False False False False False False False False False  True\n",
      "  True]\n",
      "[u'41478_cyathophylloides_SRR1754722' u'41954_cyathophylloides_SRR1754721']\n"
     ]
    }
   ],
   "source": [
    "## pca.pcs is a property of the pca object that is populated after the plot() function is called. It contains\n",
    "## the first 10 PCs for each sample. We construct a 'mask' based on the value of PC2, which here is the '1' in\n",
    "## the first line of code (numpy arrays are 0-indexed and it's typical for PCs to be 1-indexed)\n",
    "mask = pca.pcs.values[:, 1] > 500\n",
    "print(mask)\n",
    "\n",
    "## You can see here that the mask is a list of booleans that is the same length as the number of samples.\n",
    "## We can use this list to print out the names of just the samples of interest\n",
    "print(pca.samples_vcforder[mask])"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 14,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "[u'29154_superba_SRR1754715' u'30556_thamno_SRR1754720'\n",
      " u'30686_cyathophylla_SRR1754730' u'32082_przewalskii_SRR1754729'\n",
      " u'33413_thamno_SRR1754728' u'33588_przewalskii_SRR1754727'\n",
      " u'35236_rex_SRR1754731' u'35855_rex_SRR1754726' u'38362_rex_SRR1754725'\n",
      " u'39618_rex_SRR1754723' u'40578_rex_SRR1754724']\n"
     ]
    }
   ],
   "source": [
    "## We can then use this list of \"bad\" samples in a call to pca.remove_samples\n",
    "## and then replot the new pca\n",
    "pca.remove_samples(pca.samples_vcforder[mask])\n",
    "\n",
    "## Lets prove that they're gone now\n",
    "print(pca.samples_vcforder)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 15,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Using default cmap: Spectral\n"
     ]
    },
    {
     "data": {
      "text/plain": [
       "<matplotlib.axes._subplots.AxesSubplot at 0x7fe0f8c25410>"
      ]
     },
     "execution_count": 15,
     "metadata": {},
     "output_type": "execute_result"
    },
    {
     "data": {
      "image/png": "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+y2v48ePuwLAuHHjTt1///29J06ceG4g29LSUpeePXtWWq1W2bBhQyd7u6+vb/WZM2ca/bwdMmTI2bVr13YCgL179xry8vI8+vfvXw4A3333XYfjx4+7nj17VjZv3hwwbNiwsw0dy9vbWw0bNqxoxowZPSdPntzgYLunT5929ff3r/bz86v56aefPH/++WefxmptDQwoInIqpoiIigVr3st+deOWrAVr3stuSjgBQGxsbHliYmLe0KFDI0wmk/nhhx/uAQDx8fGnzpw54xYfH19g3zYpKenooEGDIocOHWoMDw8vt7dPmDChYNGiRUGRkZHnOknU5cknnzxRXV0tRqPRfM899/RZunRpjpeXl9LqOHvPPfeERUdHR918882F1157bWljtU+aNKkAAO64444zDW135513FlVVVYnRaDQ//fTT3QYMGFDS+DvT8jjdBhHpRluabmPFihUdN27cGPDZZ5+12NxTdosWLQpMTk72Wb169eGL2e/ZZ5/tWlRU5Praa6/pek6q+qbb0O0XtIiI9CouLq7H9u3b/b/44ov9jq6lPiNHjuxz6NAhw44dO7IcXculYkAREV2kVatWHQFwpLXON23atFMAGuyFd6GtW7cebKFyWg2fQRERkS4xoIiISJcYUEREpEsMKCIi0iUGFBE5FYsl02PWY5PDZk250zjrsclhFkvTRjI/efKk68svv9wFsI24cN111/VtnkqpMQwoInIaFkumx5I5U42zb/Dp9PdxPfxm3+DTacmcqcamhNSpU6dcly9ffllz1km/DwOKiJzGu4tfDpk70Wzw8bINXu7j5Y65E82Gdxe/fMnTbSQmJnY/cuSIISIiwpyUlNS9pKTEdcyYMb3DwsKibrnllrCamhoAwMyZM4Ojo6Mjw8PDo+69995e9vZBgwaZ4uPje8TGxpp69+4dtWPHDu9Ro0b16dWrV/S0adO6AbaJB3v37h114bQZAPDDDz94DRgwIMJoNJpHjhzZJz8/3+GjjLcWBhQROY/KYnd7ONn5eLkDFcWXPN3GggULcnv06GHNzMxMf/nll3MzMjK83njjjSMHDhxIO3z4sGHr1q2+APCXv/zlRGpqasb+/fvTysrKXDZs2OBvP4aHh0dNcnKy5f7778+/++67+7799tuHMzMz095///3O9jH+6po2AwAmT54c9tJLL+VmZWWlR0VFlc2aNavbpf4ubQ0Dioich7tfZUnZ+TNrlJRVAh5+lzzdxoX69etX0qdPn0pXV1dERUWVHjx40AMAtmzZ4te/f/8Io9Fo/uGHH/xSU1PPzaB7++23nwaAAQMGlPXt27esV69elV5eXqpHjx7WX375xQOoe9qMU6dOuRYXF7veeOONZwHgz3/+86ndu3f7NtfvoncMKCJyGg88mvTr7LXpVntIlZRVYvbadOsDjyZd8nQbFzIYDLWnt0BVVZWUlpZKYmJir08++eRgVlZW+sSJE0+Wl5ef+3y1T63h4uJy3v4uLi6wT49R37QZ7RkDioichskUUfHws29mzd1cUvDU+iPFczeXFDz87JtZJtOlj2ju7+9fXVJS0uBnZWlpqQsABAUFVRUVFbl8/vnnHS/1fLUFBgZWd+jQodo+Jfzy5csDBw8e3OA0G86EY/ERkVMxmSIq5r2+stlGGA8KCqqOiYk5Gx4eHmUwGGq6dOnym9uFnTt3rp4wYUK+2WyO6t69e0VzTlexYsWK7KlTp/aaNm2aS8+ePa3r16/Paa5j6x2n2yAi3WhL021Q86lvug3e4iMiIl1iQBERkS4xoIiISJcYUEREpEsMKCIi0iUGFBER6RIDioicSlaWxWPWrMfCkpIeMs6a9VhYVpalSdNtNIc1a9YEpKSkeNqXBw0aZNq5c6d3U49rsVg8wsPDoy5mn9rnDgkJ6ZeXl9fg92EHDhwYUVf7nXfeGbpixYpm+UJyffhFXSJyGllZFo8lS+YYX3xxgsHHxwslJWV45pk5Pg8//GyW0Wi65NEkmuqzzz4LqKqqKoqJiSl3VA2X6qeffsp01Ll5BUVETmP58sUh9nACAB8fL7z44gTD8uWLL3m6DQBYvHhxoNFoNJtMJvPIkSP7hISE9LNarQIABQUFLvblBQsWdI6Ojo40mUzm0aNH9ykuLnbZunWrz7Zt2wJmz57dPSIiwpyWlmYAgPXr13fs169fZGhoaLR9KKPS0lK56667Qo1GozkyMtL8+eef+wHAokWLAkeMGNFn6NCh4aGhodGJiYnB9tqqq6tx4TQdaWlpBrPZHGnfZt++fYaoqKhINOD555/vGh4eHhUeHh41Z86cc/NfeXt7DwSAmpoaTJo0qWefPn2ihg8f3vfkyZPnLnB27drlfeWVV5qioqIihwwZEn7o0CF3AJg7d+5lffr0iTIajeabbrqp98W+7wwoInIaIhXu9nCy8/HxgkjFJU+3kZyc7PnKK68E79ixI8tisaSvXbs2Z/DgwcUffPCBPwC8++67nW644YZCg8GgJkyYUJiampphsVjSTSZT2aJFizqPHDmy5Prrrz89d+7c3MzMzPSoqCgrAFRVVcm+ffsy5s2bd2TOnDndAGDevHmXAUBWVlb6e++990tCQkJoaWmpAMDevXt9Pvzww19SU1PTNm3a1Ml+m66uaTqioqKsfn5+1T/88IMXACxdurTz+PHjT9X3O+7atcv7vffeC0xJSclITk7OWL16dZfvv//+vDdyzZo1AQcOHDBYLJa0lStXHvrxxx99AcBqtcq0adN6bty48WBaWlpGXFzcyZkzZ4YAwKJFi4JSU1PTs7Ky0leuXHnoYt97BhQROQ2lPCpLSsrOayspKYNSHpc83cZXX33V4eabby4MDg6uAoCuXbtWJyQk5K9cuTIQANauXds5ISHhJACkpKR4xcTEmIxGo/njjz8OTEtL86zvuHfffXchAFx99dUlubm5HgDwww8/+E6aNOkUAAwcOLC8W7duFfv27fMEgCFDhpwJCgqq9vX1VTfeeGPht99+6wvUPU0HAEyePPnk22+/3bmqqgobN27sGB8fX29Affvtt7433HDD6Q4dOtT4+/vX3HjjjYXbt2/3q73Njh07/P70pz8VuLm5ITQ0tHLw4MHFALB3717D/v37vf74xz8aIyIizPPnzw8+evSoOwCYTKay22+/PWzJkiWd3N3dL3pcPQYUETmN+PhHf33mmXVWe0jZnkGts8bHP3rJ020opSAi5324jho1qiQ3N9fwr3/9y7e6ulquvPLKcgBISEgIW7x48eGsrKz0WbNmHbVarfV+xtqn4HBzc0N1dbXYz1UfEalzub5pOuLi4gq3b9/uv2HDhoB+/fqVBgUFVTf0O/4eF9ag7St9+/Yty8zMTM/MzEzPyspK//777/cDwPbt2/c/8sgj+SkpKT4DBgwwV1Ze3N8JDCgichpGo6ni4YefzZozZ0tBUtKG4jlzthQ0tYPEmDFjzmzatKmTfebb48ePuwLAuHHjTt1///29J06ceG4g29LSUpeePXtWWq1W2bBhQyd7u6+vb/WZM2ca/bwdMmTI2bVr13YCbFcmeXl5Hv379y8HgO+++67D8ePHXc+ePSubN28OGDZsWIPTbnh7e6thw4YVzZgxo+fkyZMbHGz3j3/849nNmzcHFBcXu5w5c8Zl8+bNHa+77rri2tsMGzas+MMPP+xUVVWFQ4cOue/evdsPAPr3719eUFDgtm3bNh/AdssvOTnZs7q6GgcPHvS4+eabi5csWZJbXFzsWlRUdFHT1Tu0F5+IBAB4B0A0AAXgAQAWAO8DCAWQA+BPSqlCsUX3awBuAFAKYLJS6kcHlE1EOmY0mirmzXu92abbiI2NLU9MTMwbOnRohIuLi4qOji79+OOPc+Lj40/NmzcvJD4+vsC+bVJS0tFBgwZFhoSEVERGRpaePXvWFQAmTJhQMHXq1NC33nqr60cffXSwvnM9+eSTJ+67775eRqPR7OrqiqVLl+Z4eXkprY6z99xzT1hOTo7nnXfeeeraa68ttVga7kI/adKkgi1btnS84447zjS03ZAhQ0rHjx9/6oorrogEgPvuuy//mmuuOe9e6X333Xf666+/7mAymaLCwsLKBw0aVAzYrgQ3bNhwcNq0aT2Li4tdq6urZerUqcf79etnHT9+fFhxcbGrUkoeeuih4507d673Kq4uDp1uQ0RWAdillHpHRDwAeAN4GkCBUuplEUkC0FEpNUtEbgDwGGwB9QcAryml/tDQ8TndBlHb0pam21ixYkXHjRs3Bnz22WfNFob1WbRoUWBycrLP6tWrD1/Mfs8++2zXoqIi19dee+1oS9XWHOqbbsNhV1Ai0gHAtQAmA4BSqgJAhYjcCmC4ttkqAN8CmAXgVgCrlS1Rd4tIgIgEK6XyWrl0Imrn4uLiemzfvt3/iy++2O/oWuozcuTIPocOHTLs2LEjy9G1XCpH3uLrDSAfwAoRGQAgBcDjALraQ0cplSci9v74IQCO1No/V2s7L6BEJAFAAgD07NmzRX8BImqfVq1adQTnfx61qGnTpp0CUG8vvLps3bq13luJbYUjO0m4AbgCwJtKqYEASgAkNbD9b7uP2J5bnd+g1DKlVKxSKrZLly7NUykREbU6RwZULoBcpdR/tOWPYAus4yISDADavydqbd+j1v7dAej6vioREV06hwWUUuoYgCMiYtKaRgBIB7AJQJzWFgdgo/Z6E4BJYnMVgCI+fyIicl6OHiz2MQDrtB58vwC4H7bQ/EBE4gEcBnC3tu1m2HrwHYCtm/n9rV8uERG1Fod+UVcptUd7XtRfKXWbUqpQKXVKKTVCKRWu/VugbauUUo8opfoopfoppdh/nIh+w2KxeDyWOCVsyvTJxscSp4Q19l2htqaxKTbqW79z507vyZMn9wCAdevW+T/99NNBLVlnc3D0FRQRUbOxWCwec15NMk78y1iDl7cBZaVWzJmf5PPsjJezTKaWnW6jqqoKbm76/Ui99tprS6+99tpSAJgwYUIRgCIHl9QoDnVERE5j8bKFIfZwAgAvbwMm/mWsYfGyhU2absNisXiEhYVF3XHHHaFGo9E8ZsyY3sXFxS4hISH9Zs6cGRwTE2N65513OkVERJjtP66urjFZWVkeR48edRs9enSf6OjoyOjo6Mh///vfPgBgNBrNJ0+edK2pqUFAQMDlixcvDgSA2267Leyzzz7zs1gsHjExMSaz2RxpNpsjt27d6nNhXcnJyZ79+vWLjIiIMBuNRvO+ffsMtdenp6d7REZGmnfs2OH9xRdf+F133XV9AdsXfydNmqT77+EwoIjIaVTWlLvbw8nOy9uAipqyS55uwy4nJ8dzypQp+VlZWel+fn418+fP7wIAnp6eNSkpKZYpU6YU2AdMjYuLyx89enSh0WiseOihh3rMmDHjeGpqasann356cMqUKaGAbeiibdu2+aakpHh2797d+t133/kCwE8//eRz3XXXlXTr1q1q165dWenp6Rnvv//+L9OnT/9NoLz++utdHn744eOZmZnpe/fuzQgLCzt3lfjzzz8b7rzzzr7Lly/PHjZsWGlTf39H0O/1KBHRRXJ38awsK7WidkiVlVrh4eJ1ydNt2AUFBVWMGjWqBADuu+++U4sWLboMACZNmlRYe7t///vfPqtXr+6ye/fuTAD4/vvvO+zfv//c3Epnz551LSwsdBk6dOjZHTt2+Obk5Hg8+OCDJ1asWNElOzvb3d/fv8rf37/m1KlTrvHx8b3S09O9XFxccOjQofOTF8DgwYNLXnnlleDc3FyPcePGFfbr188KAAUFBW633XZb3w8//PBgbGxsm5vF145XUETkNB5NmP7r2vlbrGWlVgC2cFo7f4v10YTplzzdhl190134+fnV2NsOHTrk/tBDD4W+//77B/39/WsA21QWycnJGfarqxMnTuzt2LFjzciRI4t3797t9/333/uOGjWqODAwsGrt2rUdr7rqqrMA8OKLL3a97LLLKjMyMtL37duX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\n",
      "text/plain": [
       "<Figure size 432x360 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "## and do the plot\n",
    "pca.plot()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Inspecting PCs directly\n",
    "At any time after calling plot() you can inspect the PCs for all the samples using the `pca.pcs` property. The PC values are saved internally in a convenient pandas dataframe format."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 17,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>PC1</th>\n",
       "      <th>PC2</th>\n",
       "      <th>PC3</th>\n",
       "      <th>PC4</th>\n",
       "      <th>PC5</th>\n",
       "      <th>PC6</th>\n",
       "      <th>PC7</th>\n",
       "      <th>PC8</th>\n",
       "      <th>PC9</th>\n",
       "      <th>PC10</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>29154_superba_SRR1754715</th>\n",
       "      <td>-143.458</td>\n",
       "      <td>344.601</td>\n",
       "      <td>-9.146</td>\n",
       "      <td>654.063</td>\n",
       "      <td>-71.953</td>\n",
       "      <td>-7.616</td>\n",
       "      <td>-19.466</td>\n",
       "      <td>44.390</td>\n",
       "      <td>-52.568</td>\n",
       "      <td>-8.116</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>30556_thamno_SRR1754720</th>\n",
       "      <td>-194.318</td>\n",
       "      <td>-181.059</td>\n",
       "      <td>-348.673</td>\n",
       "      <td>-94.304</td>\n",
       "      <td>-212.550</td>\n",
       "      <td>-492.266</td>\n",
       "      <td>-199.647</td>\n",
       "      <td>54.872</td>\n",
       "      <td>-71.137</td>\n",
       "      <td>-5.081</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>30686_cyathophylla_SRR1754730</th>\n",
       "      <td>-171.720</td>\n",
       "      <td>783.009</td>\n",
       "      <td>21.897</td>\n",
       "      <td>-354.809</td>\n",
       "      <td>23.015</td>\n",
       "      <td>-0.905</td>\n",
       "      <td>4.389</td>\n",
       "      <td>15.448</td>\n",
       "      <td>-19.187</td>\n",
       "      <td>-3.718</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>32082_przewalskii_SRR1754729</th>\n",
       "      <td>693.254</td>\n",
       "      <td>-18.583</td>\n",
       "      <td>-4.085</td>\n",
       "      <td>35.981</td>\n",
       "      <td>527.664</td>\n",
       "      <td>-210.055</td>\n",
       "      <td>-10.588</td>\n",
       "      <td>19.116</td>\n",
       "      <td>-22.978</td>\n",
       "      <td>-3.683</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>33413_thamno_SRR1754728</th>\n",
       "      <td>-126.793</td>\n",
       "      <td>-59.102</td>\n",
       "      <td>-29.833</td>\n",
       "      <td>24.647</td>\n",
       "      <td>4.006</td>\n",
       "      <td>-17.379</td>\n",
       "      <td>8.998</td>\n",
       "      <td>-339.049</td>\n",
       "      <td>438.306</td>\n",
       "      <td>-32.892</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>33588_przewalskii_SRR1754727</th>\n",
       "      <td>881.139</td>\n",
       "      <td>-8.878</td>\n",
       "      <td>5.835</td>\n",
       "      <td>-53.687</td>\n",
       "      <td>-434.127</td>\n",
       "      <td>170.774</td>\n",
       "      <td>6.425</td>\n",
       "      <td>3.491</td>\n",
       "      <td>-3.660</td>\n",
       "      <td>-1.877</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>35236_rex_SRR1754731</th>\n",
       "      <td>-187.931</td>\n",
       "      <td>-165.702</td>\n",
       "      <td>-163.637</td>\n",
       "      <td>-47.395</td>\n",
       "      <td>148.425</td>\n",
       "      <td>430.936</td>\n",
       "      <td>-459.261</td>\n",
       "      <td>35.808</td>\n",
       "      <td>-54.179</td>\n",
       "      <td>-5.964</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>35855_rex_SRR1754726</th>\n",
       "      <td>-184.338</td>\n",
       "      <td>-161.701</td>\n",
       "      <td>-164.247</td>\n",
       "      <td>-36.742</td>\n",
       "      <td>41.453</td>\n",
       "      <td>125.039</td>\n",
       "      <td>357.653</td>\n",
       "      <td>-286.551</td>\n",
       "      <td>-318.039</td>\n",
       "      <td>-8.572</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>38362_rex_SRR1754725</th>\n",
       "      <td>-201.661</td>\n",
       "      <td>-205.271</td>\n",
       "      <td>502.125</td>\n",
       "      <td>-54.539</td>\n",
       "      <td>-41.762</td>\n",
       "      <td>-76.632</td>\n",
       "      <td>-30.824</td>\n",
       "      <td>58.575</td>\n",
       "      <td>-66.826</td>\n",
       "      <td>-260.359</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>39618_rex_SRR1754723</th>\n",
       "      <td>-175.793</td>\n",
       "      <td>-160.807</td>\n",
       "      <td>368.111</td>\n",
       "      <td>-31.844</td>\n",
       "      <td>-28.502</td>\n",
       "      <td>-56.008</td>\n",
       "      <td>-17.545</td>\n",
       "      <td>16.067</td>\n",
       "      <td>-16.588</td>\n",
       "      <td>337.616</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>40578_rex_SRR1754724</th>\n",
       "      <td>-188.110</td>\n",
       "      <td>-166.450</td>\n",
       "      <td>-178.318</td>\n",
       "      <td>-41.402</td>\n",
       "      <td>44.339</td>\n",
       "      <td>134.100</td>\n",
       "      <td>359.870</td>\n",
       "      <td>377.916</td>\n",
       "      <td>186.759</td>\n",
       "      <td>-7.397</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "</div>"
      ],
      "text/plain": [
       "                                   PC1      PC2      PC3      PC4      PC5  \\\n",
       "29154_superba_SRR1754715      -143.458  344.601   -9.146  654.063  -71.953   \n",
       "30556_thamno_SRR1754720       -194.318 -181.059 -348.673  -94.304 -212.550   \n",
       "30686_cyathophylla_SRR1754730 -171.720  783.009   21.897 -354.809   23.015   \n",
       "32082_przewalskii_SRR1754729   693.254  -18.583   -4.085   35.981  527.664   \n",
       "33413_thamno_SRR1754728       -126.793  -59.102  -29.833   24.647    4.006   \n",
       "33588_przewalskii_SRR1754727   881.139   -8.878    5.835  -53.687 -434.127   \n",
       "35236_rex_SRR1754731          -187.931 -165.702 -163.637  -47.395  148.425   \n",
       "35855_rex_SRR1754726          -184.338 -161.701 -164.247  -36.742   41.453   \n",
       "38362_rex_SRR1754725          -201.661 -205.271  502.125  -54.539  -41.762   \n",
       "39618_rex_SRR1754723          -175.793 -160.807  368.111  -31.844  -28.502   \n",
       "40578_rex_SRR1754724          -188.110 -166.450 -178.318  -41.402   44.339   \n",
       "\n",
       "                                   PC6      PC7      PC8      PC9     PC10  \n",
       "29154_superba_SRR1754715        -7.616  -19.466   44.390  -52.568   -8.116  \n",
       "30556_thamno_SRR1754720       -492.266 -199.647   54.872  -71.137   -5.081  \n",
       "30686_cyathophylla_SRR1754730   -0.905    4.389   15.448  -19.187   -3.718  \n",
       "32082_przewalskii_SRR1754729  -210.055  -10.588   19.116  -22.978   -3.683  \n",
       "33413_thamno_SRR1754728        -17.379    8.998 -339.049  438.306  -32.892  \n",
       "33588_przewalskii_SRR1754727   170.774    6.425    3.491   -3.660   -1.877  \n",
       "35236_rex_SRR1754731           430.936 -459.261   35.808  -54.179   -5.964  \n",
       "35855_rex_SRR1754726           125.039  357.653 -286.551 -318.039   -8.572  \n",
       "38362_rex_SRR1754725           -76.632  -30.824   58.575  -66.826 -260.359  \n",
       "39618_rex_SRR1754723           -56.008  -17.545   16.067  -16.588  337.616  \n",
       "40578_rex_SRR1754724           134.100  359.870  377.916  186.759   -7.397  "
      ]
     },
     "execution_count": 17,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "pca.pcs"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Looking at PCs other than 1 & 2\n",
    "PCs 1 and 2 by definition explain the most variation in the data, but sometimes PCs further down the chain can also be useful and informative. The plot function makes it simple to ask for PCs directly."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 11,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Using default cmap: Spectral\n"
     ]
    },
    {
     "data": {
      "text/plain": [
       "<matplotlib.axes._subplots.AxesSubplot at 0x7fa3d05fd190>"
      ]
     },
     "execution_count": 11,
     "metadata": {},
     "output_type": "execute_result"
    },
    {
     "data": {
      "image/png": "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0tJ83NzcsHPnzjOTF0C/fv0KXnzxxeA9e/Z4DR8+PK9Hjx52AMjNzfW46aabur///vs7YmNjG90qvk48gzIZv5A2KEM5B0oQ1cJD8Y/uXTHnC3tRoR2AI5xWzPnC/lD8o7VebsOpquUuAgICyp1tO3fu9Lz//vut77777o7AwMBywLGURVJSUrrz7OrQoUNbWrVqVT5o0KD8zZs3B/zwww/+11xzTX5QUFDpihUrWl1yySUnAeD5559vd8EFF5Skp6enbd26Na2kpOQPn9fjxo3LXbNmzXYfH5/yoUOH2tauXRtg1FQWHBx8yrlmVGPFgDKZCTOm4oB/+envoJw4UIKoZuHh4aeemfhC1roFP+WumvPv/HULfsqtqwES+/fv93IuKfHOO++0vvTSS89Y7sJut8stt9zSdcaMGXt79uxpd7b379//hHOlXABwrnLbvXv3kry8PI/s7GzvqKioU/369Tv52muvtb/88stPAsDx48fdg4ODS9zd3TF//vygsrI/TgRufMdkf+qppw5dc801x3755RcfAPD09NQvv/xyx6pVq4Jef/311n94YiPBgDIZa2goXv58Bd7z333GQImvu9kxYcZUF1dHZH7h4eGnXp37evbCl5dmvTr39ey6Gr3XtWvX4rfeeivIZrNF5eXleUyaNOmMqWo2bNjgl5KS4jdz5swOzoESOTk5nosWLdr9008/+dlstqhu3bpFz5s37/QMAhdeeGFBaGhoMQBceeWV+YcOHfK8+uqr8wHgkUceObRq1aqgXr16RWRlZXn7+PiU4yzLly9vbbPZoiMiIqK2bdvmff/995/uZmnRokX5V199tX3evHntVqxY0bIu3oOGxuU2TConOxuJT89Cwb6j8OsQhAkzpsIaGurqsojqlVmX28jMzPS67rrrwrZt25bqyjqaKtMtt0HVs4aGckAEETVr7OIjIqpBeHj4KZ49NTwGFBERmRIDioiITIkBRUREpsSAIiIiU+IoPiJqUjIyM71mvvJKyDG73bOlxVLy1N//vjeinmcyr2u+vr69CwsLf3Z1Ha7GMygiajIyMjO9xj79lK1gQP/WPtddG1AwoH/rsU8/ZctoJGtClZeXo7IZI5orBhQRNRkzX3klJPjmmy0e3t4AAA9vbwTffLNl5iuvnNdyGydOnHC78soru4eHh0eFhYVFv/HGG61CQkJ6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\n",
      "text/plain": [
       "<Figure size 432x360 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "## Lets reload the full dataset so we have all the samples\n",
    "pca = ipa.pca(vcffile, pops_dict)\n",
    "pca.plot(pcs=[3,4])"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 12,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Using default cmap: Spectral\n",
      "Using default cmap: Spectral\n"
     ]
    },
    {
     "data": {
      "text/plain": [
       "<matplotlib.axes._subplots.AxesSubplot at 0x7fa3d0a04290>"
      ]
     },
     "execution_count": 12,
     "metadata": {},
     "output_type": "execute_result"
    },
    {
     "data": {
      "image/png": "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PvnkEy8AWL16tffdd99dpNFo1KhRo4r27t2bnZubmxUWFlaxePFi36FDh5YNGTKkePbs2UdycnKyoqKijABgMplkz5492Rb54kREdsgaORuwft6eO3fu4VmzZnUBgLlz53YCgLy8vKyPPvroYGJioq68vFwAYPfu3e6ffvrpwb1792Zu2rTJu27I36FDh7STJk06tX///kwvL6/qtWvXdoiKijJ6enpW//jjj64AsHLlSt9HHnnk7LW+54cffmjIzMzM3rVrV9bKlSs7nzhxwrEpv1tT2fpO1iIAW5RS4QBiAGQDmAng30qpUAD/Ni8DwHAAoea/RAArmiPAggIDkpNnICnpaSQnz0BBgaE5DktE1GTqXJFz3cm6jruTE9T5IufG9vn111+3GzFiRJG/v78JADp37lydmJh4OiUlxQcA1q1b55uYmHgGANLT011jY2PD9Hp95GeffeaTmZmpvVq/Dz30UFFjYyIiag2skbMB6+ftW265pezIkSMuAPDjjz96jBkz5iwA9O7du7JLly4X9uzZowWA/v37n/fz86v28PBQ99xzT9GOHTs8ACAgIMB4yy23VJj3KTcYDBoAiI+PP/Puu+/6mkwmbNy4sUNCQsI1i6y5c+d2DgsLi4yNjY04ceKE87Vibw42K7JEpB2A2wCsAgCl1AWlVDGAPwJ437zZ+wDuM3/+I4C1qtZPANqLiD+sqKDAgCVLkjBt2kAkJ9+HadMGYsmSJBZaRGQXxKtDVZnJdElbmckEadehqrF9KqUgIqp+25133ll25MgRzRdffOFRXV0tN998cyUAJCYmBi1duvRQXl5e1owZM44ZjcarnnO0Wq262joiorbAGjkbsH7ednJyQnV1tdQd62pE5IrLLi4uF3dydHRUJpNJAGDs2LFF27dv90pNTW3fs2fPcj8/v+qr9b1582bPnTt3eqalpeXk5uZmRUREVFRUVNj0ZpItD94DwGkAa0TkFxF5T0TcAXRWSh0HAPO/nczbBwA4XG//I+a2S4hIooikiUja6dOnmxRgSsoKJCc/DHd3VwCAu7srkpMfRkpKs9xEIyJqksQXXj4691Cxse6kXWYyYe6hYmPiCy8fbWyfw4YNO79p0ybvumEYJ0+edASAkSNHnn3sscd6jB49+kzdtuXl5Q7du3evMhqNkpqa6l3X7uHhUX3+/Hlbj6QgImpR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\n",
      "text/plain": [
       "<Figure size 864x360 with 2 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "import matplotlib.pyplot as plt\n",
    "\n",
    "fig = plt.figure(figsize=(12, 5))\n",
    "ax1 = fig.add_subplot(1, 2, 1)\n",
    "ax2 = fig.add_subplot(1, 2, 2)\n",
    "\n",
    "pca.plot(ax=ax1, pcs=[1, 2])\n",
    "pca.plot(ax=ax2, pcs=[3, 4])"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "It's nice to see PCs 1-4 here, but it's kind of stupid to plot the legend twice, so we can just turn off the legend on the first plot."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 13,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Using default cmap: Spectral\n",
      "Using default cmap: Spectral\n"
     ]
    },
    {
     "data": {
      "text/plain": [
       "<matplotlib.axes._subplots.AxesSubplot at 0x7fa3d0a8db10>"
      ]
     },
     "execution_count": 13,
     "metadata": {},
     "output_type": "execute_result"
    },
    {
     "data": {
      "image/png": "iVBORw0KGgoAAAANSUhEUgAAA1kAAAFACAYAAABOR7ZJAAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAAALEgAACxIB0t1+/AAAADl0RVh0U29mdHdhcmUAbWF0cGxvdGxpYiB2ZXJzaW9uIDIuMi4yLCBodHRwOi8vbWF0cGxvdGxpYi5vcmcvhp/UCwAAIABJREFUeJzs3XlcVPX+P/DXm3XYRAEXxGVQGFYllfxmatpVS1tu2+1maYrRJc20UlNbSbOupGapaVoUbkWLLWpmV7umlT/vvVC5sLowKq6oqDADIwOf3x8MhoagMnBYXs/Hg4fnfM7nnPOauY/H3N7z+czniFIKREREREREZB8OWgcgIiIiIiJqSlhkERERERER2RGLLCIiIiIiIjtikUVERERERGRHLLKIiIiIiIjsiEUWERERERGRHbHIIiIiIiIisiMWWURERERERHbEIouIiIiIiMiOnLQOUJf8/PyUXq/XOgYR0XVJTU09pZRqrXWO+sTPbSJqrJrjZzZdWZMusvR6PVJSUrSOQUR0XUTkoNYZ6hs/t4mosWqOn9l0ZZwuSEREREREZEcssoiI6KqJiKOI/CYi6237gSLyHxHZKyKfioiLrd3Vtr/PdlyvZW4iIqL6xCKLiIiuxdMAMirtJwCYr5QKBpAPINbWHgsgXykVBGC+rR8REVGz0KR/k0VERPYjIh0A3AngdQCTREQA/AXAI7YuywG8CmAJgHts2wDwBYBFIiJKKVWfmYmItJSamtrGycnpAwCR4OBGU1MGYI/Van28V69eJy8/yCKLiIiu1tsApgLwsu37AjirlLLa9nMBBNi2AwAcBgCllFVEztn6n7r8oiISByAOADp16lRn4YmI6puTk9MH7dq1C2vdunW+g4MDv2RqQsrKyiQvLy/8+PHjHwD46+XHWVETEVGNROQuACeVUqmVm6voqq7i2KWNSi1TSkUrpaJbt+bqx0TUpES2bt36PAuspsfBwUG1bt36HMpHKf+EI1lERHQ1+gL4q4jcAUAHoAXKR7ZaioiTbTSrA4Cjtv65ADoCyBURJwDeAM7Uf2wiIk05sMBqumz/21Y5aMUiq5KDB41ISlqCsrJiODjoEBMzDp0767WORUSkOaXU8wCeBwARGQhgilJqhIh8DuBvAJIBjAbwje2Utbb9/2c7/m/+HgswGnOQtGQeyorPwkHXEjHjJkOvD9Q6FhER2RmLLJuDB41YuDAeM2Y8BA8PN5hMRYiPj8eECTNYaBERXdk0AMkiMgvAbwASbe2JAFaKyD6Uj2AN1yhfg2E05mBh/Fi8OjwIHm5+MBWV4NX4sZgw4z0WWkRETQx/k2WTlLTkYoEFAB4ebpgx4yEkJS3ROBkRUcOilPpRKXWXbfuAUqq3UipIKfWgUspiay+27QfZjh/QNrX2kpbMsxVYzgAADzdnvDo8CElL5mmcjIgaiqyMTJex9zwcOPr/hhrG3vNwYFZGpos9r19WVobS0lJ7XpKugEWWTVlZ8cUCq4KHhxvKyoo1SkRERE1JWfHZiwVWBQ83Z5QVn9UoERE1JFkZmS7Tbh9huGHtMZ+B/yvxumHtMZ9pt48w1LbQysrKcunSpUvEyJEjO0VERIQvXrzY94YbbggNDw8PGzZsWJdz5845nD592lGv10fu3LnTFQDuvvvuwHnz5vnZ55U1TyyybBwcdDCZii5pM5mK4OCg0ygRERE1JQ66ljAVlVzSZioqgYOupUaJiKghmT99RsDQw16uruIIAHAVRww97OU6f/qMgBpOrZHRaNSNGTPm9L///e/s5cuX+23bti07PT09o2fPnubXXnutra+vb+n8+fMPjR49OnDZsmWtzp496zR58uQ/PXKDrh6LLJuYmHGIj//0YqFV/pusTxETM07jZERE1BTEjJuMV5P3XSy0TEUleDV5H2LGTdY4GRE1BEXH850rCqwKruKIohP5zlc45ar5+/tfGDRokOnHH3/02L9/v653796hoaGh4cnJyb6HDh1yAYD77rvvfFhYWNHUqVM7JyUlGWt7z+aOC1/YdO6sx4QJMzB37h+rC3LRCyIishe9PhATZryHeUvmoaz4FBx0LbnoBRFd5NauVYlFHUPlQsuiSuHWtk1JNaddFXd39zIAUEqhX79+59etW5dzeZ/S0lJkZ2frXF1dy06dOuXUtWvXWt+3OWORVUnnznrExydoHYOIiJoovT4QryYs0joGETVAz86OPzLttxEeFVMGLaoUGzsWWBJmv3fEXvcYOHCgafLkyZ327NnjGhkZaSkoKHDIyclx7t69u2XmzJltDQZD8euvv34kNjZWn5qamunq6trsH71xvVhkERERERFpLCQs9ELC96uz50+fEVB0It/ZrW2bkoTZ7x0JCQu9YK97tG/f3rp06VLj8OHDu1y4cEEAID4+/ggArFy50i81NTWjVatWZV988UXB9OnT/efPn3+0+ivSlUhTfjZkdHS0SklJ0ToGEdF1EZFUpVS01jnqEz+3iaixquoze+fOncaoqCguINGE7dy50y8qKkp/eTsXviAiIiIiIrIjFllERERERER2xCKLiIiIiIjIjlhkERERERER2RGLLCIiIiIiIjtikUVERERERGRHLLKIiIiIiOiKVq5c2TI1NVVXsd+7d++Qbdu2udf2ullZWS7BwcER13JO5XsHBAR0O3bsWIN87m+DDEVERERE1NxkZWa6LHt9ZoA6l+8s3q1K4l585UhIqP0eRny9vv7665ZWq/Vcr169irXO0lhwJIuIiIiISGNZmZku82NHGp6ynvSZ7q28nrKe9JkfO9KQlZnpUttrL1q0yNdgMISHhISEDxkypGtAQEA3i8UiAHDmzBmHiv158+b5RUZGhoWEhITffvvtXQsKChw2bdrksXnz5pYvvfRSh9DQ0PC0tDRXAPjkk09adevWLUyv10du3LjREwDMZrP87W9/0xsMhvCwsLDwdevWeQHAggULfAcNGtS1f//+wXq9PnLy5Mn+FdlKS0sxfPjwzkFBQRF9+/YNLiwslLS0NNfw8PCwij67d+92jYiICEM1Bg8e3DUiIiIsKCgoYu7cuX61fc9qS9MiS0RaisgXIpIpIhki0kdEfERkk4jstf3bytZXRGSBiOwTkV0i0lPL7ERERERE9rLs9ZkB0zq1dPVwKp9o5uHkhGmdWroue31mQG2um5KSops7d67/1q1bs7OystJXrVpl7NOnT8Fnn33mDQAffvihzx133JHv6uqqRowYkb9nz56MrKys9JCQkKIFCxb4DRkyxDR48OCzs2bNys3MzEyPiIiwAIDVapXdu3dnJCQkHJ45c2Z7AEhISGgDANnZ2ekff/zxgbi4OL3ZbBYA2LVrl8fnn39+YM+ePWlr1671qZjyd+jQId3EiRNP7tu3L83b27t0xYoVrSIiIixeXl6l27dvdwOApUuX+j3yyCOnq3udq1evNqalpWX8/vvv6UuXLm17/Phxx9q8b7Wl9UjWOwA2KqVCAUQByAAwHcAPSqlgAD/Y9gFgGIBg218cgCX1EfDgQSNmzJiG+PinMWPGNBw8aKyP2xIRERFRM6LO5TtXFFgV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8LDw8PCwTZs2eVyeKyUlRdetW7ew0NDQcIPBEL57927XysfT09NdwsLCwrdu3eq+fv16r1tvvTUIKH+w8ahRozrV5j2pDyyyiIiIiIg0tmjZ/ICRzw1zdXMvrzXc3F0x8rlhrouWzQ+o7bWNRqNu7NixednZ2eleXl5lc+bMaQ0AOp2uLDU1NWvs2LFnMjMz0zMzM9NHjx6dd/vtt+cbDIYLTzzxRMdJkyad2LNnT8ZXX321f+zYsXoAiI6OLty8ebNnamqqrkOHDpaff/7ZEwB+++03j1tvvdXUvn17608//ZSdnp6e8emnnx549tnd86lZAAAgAElEQVRn/1QULVy4sPWTTz55IjMzM33Xrl0ZgYGBF0fsdu7c6frAAw8EJSYm5gwYMMBc29evhYY7LkhERERE1EyUlBU7VxRYFdzcXXGhrMi5ttdu167dhdtuu80EAI8++ujpBQsWtAGAUaNG5Vfu969//ctjxYoVrXfs2JEJAL/88kuLvXv3ulUcLywsdMzPz3fo379/4datWz2NRqPL448/fvKjjz5qnZOT4+zt7W319vYuO336tGNsbGzn9PR0NwcHBxw8ePDSFwagT58+prlz5/rn5ua6DB8+PL9bt24WADhz5ozTvffeG/T555/vj46OLq7ta9cKR7KIiIiIiDTm7KArKTJbLmkrMlvg4uBWUttri0iV+15eXmUVbQcPHnR+4okn9J9++ul+b2/vMgBQSiElJSWjYpTr5MmTu1q1alU2ZMiQgh07dnj98ssvnrfddluBr6+vddWqVa1uuummQgB4/fXX27Zp06YkIyMjfffu3eklJSV/qjnGjh175ptvvtnn5uZWNmzYMMPatWu9bJlK/f39L/z444+etX3dWmKRRURERESksafinj2yas53lopCq8hswao531meinv2SG2vfezYMZfNmzd7AMDHH3/sc/PNNxdWPm6xWOT+++/v8tprrx3p3r37xUqvX79+5xMSEtpU7G/fvt0NAIKCgkry8/OdcnJydOHh4Rf69OlT+O6777a75ZZbCgHg3Llzjv7+/iWOjo5YvHixb2nppc8EBC7+5sry0ksvnbztttvO/v77724A4OzsrDZu3Lj/k08+8X3vvfd8avvatcIii4iIiIhIYyEhIRdemTQ7e8OSX898MuffBRuW/HrmlUmzs+2xumCXLl2KP/zwQ1+DwRCen5/vNGXKlLzKxzdv3uyxZ88ej1mzZrWvWPzCaDQ6L1u27PCvv/7qYTAYwrt27RqxaNGi1hXn3HDDDabAwMBiABg4cGDByZMnnQcPHlwAAM8888zJTz75xDcqKio0Oztb5+bmVobLrFy50sdgMESEhoaG7927V/fEE0+crjjWokWLsu+//37fokWL2q5ataplbV+/FkQppXWGOhMdHa1SUlK0jkFEdF1EJFUpFa11jvrEz20iaqyq+szeuXOnMSoq6pRWmYDy1QXvuuuu4L1796ZpmaOp2rlzp19UVJT+8nbNR7JExFFEfhOR9bb9QBH5j4jsFZFPRcTF1u5q299nO67XMjcREREREVFVNC+yADwNIKPSfgKA+UqpYAD5AGJt7bEA8pVSQQDm2/oREREREdEVhISEXOAoVv3TtMgSkQ4A7gTwgW1fAPwFwBe2LssB3Gvbvse2D9vxQXL5UilEREREREQa03ok620AUwFU/BjOF8BZpZTVtp8LoOIBbAEADgOA7fg5W38iIiIiIqIGQ7MiS0TuAnBSKZVaubmKruoqjlW+bpyIpIhISl5eXhWnEBERERER1R0tR7L6AviriBgBJKN8muDbAFqKiJOtTwcAR23buQA6AoDtuDeAM5dfVCm1TCkVrZSKbt269eWHiYiIiIiI6pRmRZZS6nmlVAellB7AcAD/VkqNALAFwN9s3UYD+Ma2vda2D9vxf6umvP48EREREVEj4+7u3kPrDA2BU81d6t00AMkiMgvAbwASbe2JAFaKyD6Uj2AN1ygfEREREZHdZWZlucx6552AsxaLc0tX15KXnn76SKgdHkZcH8rKysDxjz9ovfAFAEAp9aNS6i7b9gGlVG+lVJBS6kGllMXWXmzbD7IdP6BtaiIiIiIi+8jMynKJefklg6l/Px+3u+70MvXv5xPz8kuGzKwsl9pc9/z58w4DBw4MCgkJCQ8ODo54//33WwUEBHQ7duyYEwBs27bNvXfv3iEAMGnSpPb33ntv4E033WTo3Llz5Lx58/wqrvPyyy+3jYyMDDMYDOHPPvtse6D8QcddunSJGDlyZKeIiIjw/fv3uwDAP/7xjw7h4eFhffr0MRw9etQJAObNm+cXGRkZFhISEn777bd3LSgoaBB1SF1p0i+OiIiIiKgxmPXOOwH+993n6qTTAQCcdDr433ef66x33gmo4dRqffnlly3atWtXkpWVlb537960+++//3x1/TMyMtw2b968d8eOHZlz5sxpbzQanb/88ssW+/bt0+3atSsjIyMj/ffff3f/7rvvPAHAaDTqxowZczojIyPdYDBcKCoqcujZs6c5PT09o2/fvgXTp09vDwAjRozI37NnT0ZWVlZ6SEhI0YIFC/yqy9HYscgiIiIiItLYWYvFuaLAquCk0+GcxeJcm+v27Nmz6Keffmoxbty4gI0bN3r6+vqWVtd/2LBhZz09PZW/v7+1T58+53/66SePjRs3tti2bVuL8PDwcNuIlS4zM1MHAP7+/hcGDRpkqjjfwcEBjz/++BkAeOyxx07/97//9QSA1NRUt169eoUYDIbwNWvW+KalpemqTtA0NMTfZBERERERNSstXV1LTMXFqFxoWYuL4e3qWlKb63bv3t3y66+/pq9Zs8b7xRdfDNi8efN5R0dHVVZW/pjaoqKiSwZdRC59apKIQCmFZ5555thzzz13qvKxrKwsF3d39zJUo+J6cXFxgV988cW+Pn36FC1YsMB369atXrV5XQ0dR7KIiIiIiDT20tNPHzn21VcWa3ExgPIC69hXX1leevrpI7W5rtFodPby8ip78sknzzzzzDMnfv/9d/cOHTpc+OWXX9wB4LPPPmtVuf93333X0mw2y/Hjxx137Njh1a9fP9OwYcPOr1y50u/cuXMOAJCTk+N85MiRKgdrysrK8NFHH7UCgKSkJN/evXsXAIDZbHbo1KlTicVikeTkZJ/avKbGgCNZREREREQaCw0JuZD02qzsWe+8E3DOYnH2dnUtSXptVq1XF0xNTXV7/vnnOzg4OMDJyUktXrz4oNlsdhg7dqw+ISGhpFevXqbK/Xv06GEaNGhQ8NGjR12mTJlyTK/Xl+j1+pK0tDTdjTfeGAoA7u7uZatXr85xcnL603KCbm5uZWlpaW4RERHtvLy8Sr/88ssDADB9+vSjvXv3DgsICLgQFhZmLiwsdKzN62ropCkvtRgdHa1SUlK0jkFEdF1EJFUpFa11jvrEz20iaqyq+szeuXOnMSoq6tSVzmloJk2a1N7T07N05syZJ7TO0ljs3LnTLyoqSn95O6cLEhERERER2RGnCxIREREREd56662jWmdoKjiSRURERERUN8rKysqk5m7UGNn+t61ydUUWWUREREREdWNPXl6e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\n",
      "text/plain": [
       "<Figure size 864x360 with 2 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "fig = plt.figure(figsize=(12, 5))\n",
    "ax1 = fig.add_subplot(1, 2, 1)\n",
    "ax2 = fig.add_subplot(1, 2, 2)\n",
    "\n",
    "pca.plot(ax=ax1, pcs=[1, 2], legend=False)\n",
    "pca.plot(ax=ax2, pcs=[3, 4])"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Controlling colors\n",
    "You might notice the default color scheme is unobtrusive, but perhaps not to your liking. There are two ways of modifying the color scheme, one simple and one more complicated, but which gives extremely fine grained control over colors.\n",
    "\n",
    "Colors for the more complicated method can be specified according to [python color conventions](https://matplotlib.org/users/colors.html). I find [this visual page of python color names useful](https://matplotlib.org/2.0.0/examples/color/named_colors.html)."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 14,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "<matplotlib.axes._subplots.AxesSubplot at 0x7fa3d099ac50>"
      ]
     },
     "execution_count": 14,
     "metadata": {},
     "output_type": "execute_result"
    },
    {
     "data": {
      "image/png": "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HhYWFnPfffeFGo1G05AhQ7qXlJQYgoOD4yZNmhQUHx8f+e6777aPiooy2R4uLi7xubm57gcOHHAdPHhweGxsbHRsbGz0d9995wMARqPRdOzYMReLxYK2bdtePX/+/AAAuOeee8JWr17tl5OT4x4fHx9pMpmiTSZT9IYNG3wurys1NdUzLi4uOioqymQ0Gk27du3yqLo8MzPTPTo62rR582bvdevW+d166609AOsXf0eNGtW1McekKTCgiKjFEMshN5/L5qn18QbEcrjB023YFBQUeD7++ONHc3NzM/38/CyvvvpqRwDw9PS0pKWl5Tz++OPFtgFTR48efXTw4MEnjEZj+WOPPdZlwoQJh9PT07O++OKLvMcffzwUsA5dtHHjRt+0tDTPkJCQsh9//NEXAH799VefW2+9tbRz584VW7Zsyc3MzMz66KOP9o4fP/53gfLGG290HDt27OHs7OzMnTt3ZoWFhV04S9yxY4fH/fff32Px4sX5N998s7mxv78z6Pd8lIionpQh8Hyp2RpKNqVmQBk6NXi6DZvAwMDyQYMGlQLAI488cnzevHlXAcCoUaNOVF3vu+++8/nggw86bt26NRsAfvrpJ//du3dfmFvpzJkzLidOnDDceOONZzZv3uxbUFDg/uijjx5ZsmRJx/z8fLc2bdpUtGnTxnL8+HGXpKSkbpmZmV4GgwGFhYWXnB0BQL9+/Upnz54dVFRU5D58+PATcXFxZQBQXFzses899/T45JNP8hISEs419nd3Fp5BEVGLkZQ8c/+Uf4aUlWrnC9Z7UCFlSckzGzzdhk1N0134+flZbG2FhYVujz32WOhHH32U16ZNGwtgncoiNTU1y3Z2deTIkZ3t2rWzDBw4sGTr1q1+P/30k++gQYNKAgICKpYvX97uuuuuOwMAr7zySqerrrrqfFZWVuauXbsyz58//7vP68cff7x4zZo1e7y8vCxDhw41rl271k+rqTIoKKjcNmdUc8WAIqIWwxgZXT52wne50xcMK07557Ul0xcMK7ZXB4mDBw+626aU+PDDD9tff/31l0x3UVZWJvfdd1/3GTNm7O/Zs2eZrb1///6nbTPlAoBtltsePXqcP3HihGt+fr6nyWQq79ev35k333wz8KabbjoDAKdOnXIJCgo67+LiggULFgRUVv5+IHDtHlPZCy+8cGTQoEEnt2/f7gUA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iEiXGFBERKRLDCgiItIlBhQREekSA4qIiHSJAUVERLrEgCIiIl1iQBERkS4xoIiISJcYUEREpEsMKCIi0iUGFBER6RIDioiIdIkBRUREusSAIiIiXWJAERGRLjGgiIhIlxhQRESkSwwoIiLSJQYUERHpUr0CSkR8RMTFUcUQERHZ1BpQImIQkYdF5CsROQIgG8BBEckQkVdFJKJpyiQiotamrjOoTQDCATwHIFAp1UUpdRWAGwFsBTBTREY6uEYiImqFXOtYfptS6vzljUqpYgCfAfhMRNwcUhkREbVqtQbU5eEkIp4ARgLwAvChUup4dQFGRETUWPXtxfc6ABcA5wCstn85REREVnV1kvhQRMKrNLUHsALASgDtHFkYERG1bnXdg3oBwMsicgDADACzAawF4Ang744tjYiIWrO67kHtBfCwiPQH8BGArwAMVEpVNkVxRETUetV1ia+diDwJwATgjwBOAfhWRO5siuKIiKj1qquTxGoAZbBe0lumlPoAwF0A4kVkraOLIyKi1quue1ABAD6EtVv5KABQSp0F8JKIBDm4NiIiasXqCqhpADYAqASQUnWBUuqgo4oiIiKqq5PEZ7COGEFERNSk6uok8ZSIdNCeh4vIDyJyUkT+KyJxTVMiERG1RnV1knhCKXVMez4PwFylVFsAkwG87dDKiIioVasroKpeArxKKfUFACilvgfg56iiiIiI6gqoT0VkqYh0B/CFiDwjIl1FZAyAfU1QHxERtVJ1dZKYIiKJsI69Fw7AA0AyrN+PGuHw6oiIqNWqq5s5lFJLASx1eCVERERV1He6jQtEJNCehRAREVXV4IACsNhuVRAREV2mwQGllLrDnoUQERFVVe+AEpH2jiiEiIioqrpGknihynOTiOQCSBORAhG5tjE7FpEuIrJJRLJEJENE/qq1txeRDSKyW/vZTmsXEZknIntEZKeI9GnM/omISN/qOoO6r8rzVwH8VSkVBuvcUHMbue8KABOVUtEArgPwpIiYYB2U9t9KqQgA/8bFQWqHAojQHskA3mrk/hslP78QI8e8hFvvmIaRY15Cfn6hM8shImpx6uxmXkVnpdTXAKCU+kVEvBqzY2009IPa8xIRyQIQDOBuALdoq70P4HtYh1a6G8AHSikFYKuItBWRIGeMqp6fX4iB976BvHMvAQYfIK8UW++dhg1fPI2wsG5N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\n",
      "text/plain": [
       "<Figure size 432x360 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "## Here's the simple way, just pass in a matplotlib cmap, or even better, the name of a cmap\n",
    "pca.plot(cmap=\"jet\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 22,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "<matplotlib.axes._subplots.AxesSubplot at 0x7fa3d0646b50>"
      ]
     },
     "execution_count": 22,
     "metadata": {},
     "output_type": "execute_result"
    },
    {
     "data": {
      "image/png": "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du7c2ezll18O6tKliyklJcUDADZu3Ni8a9euEcHBwVG2oYzMZrM8+OCDwUaj0RQREWH6/PPP/QAgISGhZf/+/UP79OkTFhwcHDV16tQAW20lJSW4fpqOlJQUD5PJFGFb58CBAx6RkZERqMQrr7zSNiwsLDIsLCxyzpw5pfNfeXt79wAAi8WC0aNHdwgNDY28/fbbO585c6b0BGfv3r3eN910U3hkZGRE7969w3JyctwA4NVXX20TGhoaaTQaTffcc0+n6h53BhQRNRhyQtxs4WTjAx/ISanxdBuJiYmeixYtCti9e3dmRkZG6vr167N79eqV/+GHHzYFgPfee6/FXXfddc7Dw0ONHDnyXHJyclpGRkZqeHj45YSEhFYDBgwouPPOO8+/+uqruenp6amRkZGFAFBcXCwHDhxIW7BgwW9z5sxpBwALFixoAwCZmZmp77///q/jx48PNpvNAgD79+/3+eijj35NTk5O2bp1awvbZbrypumIjIws9PPzK/nhhx+8AGD58uWtHnvssbMV/Yx79+71fv/991smJSWlJSYmpq1du7b1999/71V2nXXr1jU7dOiQR0ZGRsrq1atzfvzxR18AKCwslMmTJ3fYsmXL4ZSUlLQxY8acmTZtWiAAJCQk+CcnJ6dmZmamrl69Oqe6x54BRUQNhvJXRQW4dhi8AhRAtVU1nm7j66+/bjJ06NBzAQEBxQDQtm3bkvHjx59evXp1SwBYv359q/Hjx58BgKSkJK/o6Ohwo9Fo+uSTT1qmpKR4VrTdhx566BwA3HrrrQW5ubnuAPDDDz/4jh49+iwA9OjR40q7du2uHjhwwBMAevfufdHf37/E19dX3X333ee+/fZbX6D8aToAIDY29syKFStaFRcXY8uWLc3j4uIqDKhvv/3W96677jrfpEkTS9OmTS133333uV27dvmVXWf37t1+Dz/8cJ6rqyuCg4OLevXqlQ8A+/fv9zh48KDXn/70J2OXLl1MCxcuDDh27JgbAISHh1++//77Q5YtW9bCzc2t2uPqMaCIqMGImx93dEbQjEJbSBWgADOCZhTGzY+r8XQbSimIyDUfrgMHDizIzc31+OKLL3xLSkrkpptuugIA48ePD1m6dOmRzMzM1OnTpx8rLCys8DPWNgWHq6srSkpKxLaviohIua8rmqZjzJgx53bt2tV006ZNzbp27Wr29/e/dqKs637GP+L6GrT3SufOnS+np6enpqenp2ZmZqZ+//33BwFg165dB//85z+fTkpK8unevbupqKh6vycwoIiowTBGGK8+/c3TmXPunZMXf3N8/px75+TVtoPE4MGDL27durWFbebbkydPugDAiBEjzo4dO7bTqFGjSgeyNZvNhg4dOhQVFhbKpk2bWtjafX19Sy5evFjl523v3r0vrV+/vgVgPTM5fvy4e7du3a4AwHfffdfk5MmTLpcuXZLt27c369evX6XTbnh7e6t+/fpdmDJlSofY2NhKB9v905/+dGn79u3N8vPzDRcvXjRs3769+R133JFfdp1+/frlf/TRRy2Ki4uRk5Pjtm/fPj8A6Nat25W8vDzXnTt3+gDWS36JiYmeJSUlOHz4sPvQoUPzly1blpufn+9y4cKFak1Xz158RNSgGCOMVxdsWWC36TZiYmKuTJ069XifPn26GAwGFRUVZf7kk0+y4+Lizi5YsCAwLi4uz7ZufHz8sZ49e0YEBgZejYiIMF+6dMkFAEaOHJk3ceLE4Lfffrvtxx9/fLiifb3wwgunHn/88Y5Go9Hk4uKC5cuXZ3t5eSmtjkuPPPJISHZ2tucDDzxwtm/fvuaMjIxKu9CPHj0678svv2w+fPjwi5Wt17t3b/Njjz129sYbb4wAgMcff/z0bbfddrnsOo8//vj5f/zjH03Cw8MjQ0JCrvTs2TMfsJ4Jbtq06fDkyZM75Ofnu5SUlMjEiRNPdu3atfCxxx4Lyc/Pd1FKyVNPPXWyVatWFZ7FlYfTbRCRbtSn6TZWrVrVfMuWLc02b95stzCsSEJCQsvExESftWvXHqnO+/7617+2vXDhgssbb7xxzFG12UNF023wDIqIqJrGjBnTfteuXU23bdt20Nm1VGTAgAGhOTk5Hrt37850di01xYAiIqqmNWvW/Abgt7ra3+TJk88CqLAXXnl27NhR4aXE+oKdJIiISJcYUEREpEsMKCIi0iUGFBER6RIDiogalMy0dPfpw8eHxPceZZw+fHxIZlq606fbsKeqptioaPmePXu8Y2Nj2wPAhg0bmr700kv+jqzTHtiLj4gajMy0dPdlQ583zjt/u4ePeKAgvRAzhj7v8/TnCzONEV0cOt1GcXExXF31+5Hat29fc9++fc0AMHLkyAsALji5pCrxDIqIGoyVM14LtIUTAPiIB+adv91j5YzXajXdRkZGhntISEjk8OHDg41Go2nw4MGd8vPzDYGBgV2nTZsWEB0dHf7uu++26NKli8n2cHFxic7MzHQ/duyY66BBg0KjoqIioqKiIr755hsfADAajaYzZ864WCwWNGvW7IalS5e2BID77rsvZPPmzX4ZGRnu0dHR4SaTKcJkMkXs2LHD5/q6EhMTPbt27RrRpUsXk9FoNB04cMCj7PLU1FT3iIgI0+7du723bdvmd8cdd3QGrF/8HT16dIfaHJO6wIAiogZDTpndbOFk4yMekFPmGk+3YZOdne05YcKE05mZmal+fn6WhQsXtgYAT09PS1JSUsaECRPybAOmjhkz5vSgQYPOGY3Gq0899VT7KVOmnExOTk777LPPDk+YMCEYsA5dtHPnTt+kpCTPoKCgwu+++84XAH766SefO+64o6Bdu3bFe/fuzUxNTU374IMPfn3uued+FyhLlixp/fTTT59MT09P3b9/f1pISEjpWeIvv/zi8cADD3ReuXJlVr9+/cy1/fmdQb/no0RE1aTaeBcVpBeibEgVqEKoNt41nm7Dxt/f/+rAgQMLAODxxx8/m5CQ0AYARo8efa7set98843P2rVrW+/bty8dAL7//vsmBw8eLJ1b6dKlSy7nzp0z9OnT59Lu3bt9s7Oz3Z988slTq1atap2VleXWtGnT4qZNm1rOnj3rEhcX1zE1NdXLYDAgJyfn2uQF0KtXr4JFixYF5Obmuo8YMeJc165dCwEgLy/P9b777uv80UcfHY6Jial3s/ja8AyKiBqMuHlTjs5o9m1hgSoEYA2nGc2+LYybN6XG023YVDTdhZ+fn8XWlpOT4/bUU08Ff/DBB4ebNm1qAaxTWSQmJqbZzq5OnTq1v3nz5pYBAwbk79u3z+/777/3HThwYH7Lli2L169f3/yWW265BADz5s1r26ZNm6K0tLTUAwcOpBYVFf3u83rChAl5W7ZsOeTl5WUZMmSIcevWrX5aTSUBAQFXbXNG1VcMKCJqMIwRXa4+/fnCzDl9M/Liw/flz+mbkWevDhLHjx93t00p8f7777e49dZbr5nuorCwUIYPH95p7ty5R7t161Zoa+/du/dF20y5AGCb5bZz585F586dc83KyvI0mUxXe/XqdenNN9/079u37yUAuHDhgktAQECRi4sLli1b1rKk5PcDgWv3mApf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\n",
      "text/plain": [
       "<Figure size 432x360 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "## Here's the harder way that gives you uber control. Pass in a dictionary mapping populations to colors.\n",
    "my_colors = {\n",
    "    \"rex\":\"aliceblue\",\n",
    "    \"thamno\":\"crimson\",\n",
    "    \"przewalskii\":\"deeppink\",\n",
    "    \"cyathophylloides\":\"fuchsia\",\n",
    "    \"cyathophylla\":\"goldenrod\",\n",
    "    \"superba\":\"black\"\n",
    "}\n",
    "pca.plot(cdict=my_colors)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Dealing with missing data\n",
    "RAD-seq datasets are often characterized by moderate to high levels of missing data. While there may be many thousands or tens of thousands of loci recovered overall, the number of loci that are recovered in all sequenced samples is often quite small. The distribution of depth of coverage per locus is a complicated function of the size of the genome of the focal organism, the restriction enzyme(s) used, the size selection tolerances, and the sequencing effort. \n",
    "\n",
    "Both model-based (STRUCTURE and the like) and model-free (PCA/sNMF/etc) genetic \"clustering\" methods are sensitive to missing data. Light to moderate missing data that is distributed randomly among samples is often not enough to seriously impact the results. These are, after all, only exploratory methods. However, if missing data is biased in some way then it can distort the number of inferred populations and/or the relationships among these. For example, if several unrelated samples recover relatively few loci, for whatever reason (mistakes during library prep, failed sequencing, etc), clustering methods may erroniously identify this as true \"similarity\" with respect to the rest of the samples, and create spurious clusters.\n",
    "\n",
    "In the end, all these methods must do something with sites that are uncalled in some samples. Some methods adopt a strategy of silently asigning missing sites the \"Reference\" base. Others, assign missing sites the average base. \n",
    "\n",
    "There are several ways of dealing with this:\n",
    "\n",
    " * One method is to simply __eliminate all loci with missing data__. This can be ok for SNP chip type data, where missingness is very sparse. For RAD-Seq type data, eliminating data for all missing loci often results in a drastic reduction in the size of the final data matrix. Assemblies with thousands of loci can be pared down to only tens or hundreds of loci.\n",
    " * Another method is to __impute missing data__. This is rarely done for RAD-Seq type data, comparatively speaking. Or at least it is rarely done intentionally. \n",
    " * A third method is to __downsample using a hypergeometric projection__. This is the strategy adopted by dadi in the construction of the SFS (which abhors missing data). It's a little complicated though, so we'll only look at the first two strategies."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Inspect the amount of missing data under various conditions\n",
    "The pca module has various functions for inspecting missing data. The simples is the `get_missing_per_sample()` function, which does exactly what it says. It displays the number of ungenotyped snps per sample in the final data matrix. Here you can see that since we are using simulated data the amount of missing data is very low, but in real data these numbers will be considerable. "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "1A_0    2\n",
       "1B_0    2\n",
       "1C_0    1\n",
       "1D_0    4\n",
       "2E_0    0\n",
       "2F_0    0\n",
       "2G_0    0\n",
       "2H_0    1\n",
       "3I_0    2\n",
       "3J_0    2\n",
       "3K_0    1\n",
       "3L_0    2\n",
       "dtype: int64"
      ]
     },
     "execution_count": 4,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "pca.get_missing_per_sample()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "This is useful, but it doesn't give us a clear direction for how to go about dealing with the missingness. One way to reduce missing data is to reduce the tolerance for samples ungenotyped at a snp. The other way to reduce missing data is to remove samples with very poor sequencing. To this end, the `.missingness()` function will show a table of number of retained snps for various of these conditions."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>Full</th>\n",
       "      <th>2E_0</th>\n",
       "      <th>2F_0</th>\n",
       "      <th>2G_0</th>\n",
       "      <th>1C_0</th>\n",
       "      <th>2H_0</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0</th>\n",
       "      <td>2547</td>\n",
       "      <td>2452</td>\n",
       "      <td>2313</td>\n",
       "      <td>2093</td>\n",
       "      <td>1958</td>\n",
       "      <td>1640</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>1</th>\n",
       "      <td>2553</td>\n",
       "      <td>2458</td>\n",
       "      <td>2319</td>\n",
       "      <td>2098</td>\n",
       "      <td>1963</td>\n",
       "      <td>1643</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>3</th>\n",
       "      <td>2554</td>\n",
       "      <td>2459</td>\n",
       "      <td>2320</td>\n",
       "      <td>2099</td>\n",
       "      <td>1963</td>\n",
       "      <td>1643</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>8</th>\n",
       "      <td>2555</td>\n",
       "      <td>2460</td>\n",
       "      <td>2321</td>\n",
       "      <td>2099</td>\n",
       "      <td>1963</td>\n",
       "      <td>1643</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "</div>"
      ],
      "text/plain": [
       "   Full  2E_0  2F_0  2G_0  1C_0  2H_0\n",
       "0  2547  2452  2313  2093  1958  1640\n",
       "1  2553  2458  2319  2098  1963  1643\n",
       "3  2554  2459  2320  2099  1963  1643\n",
       "8  2555  2460  2321  2099  1963  1643"
      ]
     },
     "execution_count": 3,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "pca.missingness()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Here the columns indicate progressive removal of the samples with the fewest number of snps. So \"Full\" indicates retention of all samples. \"2E_0\" shows # snps after removing this sample (as it has the most missing data). \"2F_0\" shows the # snps after removing both this sample & \"2E_0\". And so on. You can see as we move from left to right the total number of snps goes down, but also so does the amount of missingness.\n",
    "\n",
    "Rows indicate thresholds for number of allowed missing samples per snp. The \"0\" row shows the condition of allowing 0 missing samples, so this is the complete data matrix. The \"1\" row shows # of snps retained if you allow 1 missing sample. And so on."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Filter by missingness threshold - trim_missing()\n",
    "\n",
    "The `trim_missing()` function takes one argument, namely the maximum number of missing samples per snp. Then it removes all sites that don't pass this threshold."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>Full</th>\n",
       "      <th>1A_0</th>\n",
       "      <th>1B_0</th>\n",
       "      <th>1C_0</th>\n",
       "      <th>2E_0</th>\n",
       "      <th>2F_0</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0</th>\n",
       "      <td>2547</td>\n",
       "      <td>2456</td>\n",
       "      <td>2282</td>\n",
       "      <td>2079</td>\n",
       "      <td>1985</td>\n",
       "      <td>1845</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>1</th>\n",
       "      <td>2553</td>\n",
       "      <td>2462</td>\n",
       "      <td>2286</td>\n",
       "      <td>2083</td>\n",
       "      <td>1989</td>\n",
       "      <td>1849</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "</div>"
      ],
      "text/plain": [
       "   Full  1A_0  1B_0  1C_0  2E_0  2F_0\n",
       "0  2547  2456  2282  2079  1985  1845\n",
       "1  2553  2462  2286  2083  1989  1849"
      ]
     },
     "execution_count": 7,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "pca.trim_missing(1)\n",
    "pca.missingness()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "You can see that this also has the effect of reducing the amount of missingness per sample."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 8,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "1A_0    0\n",
       "1B_0    0\n",
       "1C_0    0\n",
       "1D_0    2\n",
       "2E_0    0\n",
       "2F_0    0\n",
       "2G_0    0\n",
       "2H_0    1\n",
       "3I_0    1\n",
       "3J_0    1\n",
       "3K_0    0\n",
       "3L_0    1\n",
       "dtype: int64"
      ]
     },
     "execution_count": 8,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "pca.get_missing_per_sample()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "__NB:__ This operation is _destructive_ of the data inside the pca object. It doesn't do anything to your data on file, though, so if you want to rewind you can just reload your vcf file."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>Full</th>\n",
       "      <th>2E_0</th>\n",
       "      <th>2F_0</th>\n",
       "      <th>2G_0</th>\n",
       "      <th>1C_0</th>\n",
       "      <th>2H_0</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0</th>\n",
       "      <td>2547</td>\n",
       "      <td>2452</td>\n",
       "      <td>2313</td>\n",
       "      <td>2093</td>\n",
       "      <td>1958</td>\n",
       "      <td>1640</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>1</th>\n",
       "      <td>2553</td>\n",
       "      <td>2458</td>\n",
       "      <td>2319</td>\n",
       "      <td>2098</td>\n",
       "      <td>1963</td>\n",
       "      <td>1643</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>3</th>\n",
       "      <td>2554</td>\n",
       "      <td>2459</td>\n",
       "      <td>2320</td>\n",
       "      <td>2099</td>\n",
       "      <td>1963</td>\n",
       "      <td>1643</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>8</th>\n",
       "      <td>2555</td>\n",
       "      <td>2460</td>\n",
       "      <td>2321</td>\n",
       "      <td>2099</td>\n",
       "      <td>1963</td>\n",
       "      <td>1643</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "</div>"
      ],
      "text/plain": [
       "   Full  2E_0  2F_0  2G_0  1C_0  2H_0\n",
       "0  2547  2452  2313  2093  1958  1640\n",
       "1  2553  2458  2319  2098  1963  1643\n",
       "3  2554  2459  2320  2099  1963  1643\n",
       "8  2555  2460  2321  2099  1963  1643"
      ]
     },
     "execution_count": 9,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "## Voila. Back to the full dataset.\n",
    "pca = ipa.pca(data)\n",
    "pca.missingness()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Imputing missing genotypes\n",
    "McVean (2008) recommends filling missing sites with the average genotype of the population, so that's what we're doing here. For each population, we determine the average genotype at any site with missing data, and then fill in the missing sites with this average. In this case, if the average \"genotype\" is \"./.\", then this is what gets filled in, so essentially any site missing more than 50% of the data isn't getting imputed. If two genotypes occur with equal frequency then the average is just picked as the first one."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 12,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>Full</th>\n",
       "      <th>2E_0</th>\n",
       "      <th>2F_0</th>\n",
       "      <th>2G_0</th>\n",
       "      <th>2H_0</th>\n",
       "      <th>1C_0</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0</th>\n",
       "      <td>2553</td>\n",
       "      <td>2458</td>\n",
       "      <td>2319</td>\n",
       "      <td>2099</td>\n",
       "      <td>1779</td>\n",
       "      <td>1643</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>3</th>\n",
       "      <td>2554</td>\n",
       "      <td>2459</td>\n",
       "      <td>2320</td>\n",
       "      <td>2100</td>\n",
       "      <td>1780</td>\n",
       "      <td>1643</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>8</th>\n",
       "      <td>2555</td>\n",
       "      <td>2460</td>\n",
       "      <td>2321</td>\n",
       "      <td>2100</td>\n",
       "      <td>1780</td>\n",
       "      <td>1643</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "</div>"
      ],
      "text/plain": [
       "   Full  2E_0  2F_0  2G_0  2H_0  1C_0\n",
       "0  2553  2458  2319  2099  1779  1643\n",
       "3  2554  2459  2320  2100  1780  1643\n",
       "8  2555  2460  2321  2100  1780  1643"
      ]
     },
     "execution_count": 12,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "pca.fill_missing()\n",
    "pca.missingness()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "In comparing this missingness matrix with the previous one, you can see that indeed some snps are being recovered (though not many, again because of the clean simulated data). \n",
    "\n",
    "You can also examine the effect of imputation on the amount of missingness per sample. You can see it doesn't have as drastic of an effect as trimming, but it does have some effect, plus you are retaining more data!"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 11,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "1A_0    2\n",
       "1B_0    2\n",
       "1C_0    1\n",
       "1D_0    2\n",
       "2E_0    0\n",
       "2F_0    0\n",
       "2G_0    0\n",
       "2H_0    0\n",
       "3I_0    1\n",
       "3J_0    1\n",
       "3K_0    1\n",
       "3L_0    1\n",
       "dtype: int64"
      ]
     },
     "execution_count": 11,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "pca.get_missing_per_sample()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Dealing with unequal sampling\n",
    "Unequal sampling of populations can potentially distort PC analysis (see for example [Bradburd et al 2016](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005703)). Model based ancestry analysis suffers a similar limitation [Puechmaille 2016](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.12512)). McVean (2008) recommends downsampling larger populations, but nobody likes throwing away data. [Weighted PCA](https://www.asas.org/docs/default-source/wcgalp-proceedings-oral/210_paper_8713_manuscript_220_0.pdf?sfvrsn=2) was proposed, but has not been adopted by the community. "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 18,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "{'cyathophylla': 1,\n",
       " 'cyathophylloides': 2,\n",
       " 'przewalskii': 2,\n",
       " 'rex': 5,\n",
       " 'superba': 1,\n",
       " 'thamno': 2}"
      ]
     },
     "execution_count": 18,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "{x:len(y) for x, y in pca.pops.items()}"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Dealing with linked snps"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "\n",
    "prettier_labels = {\n",
    "    \n",
    "    \"32082_przewalskii\":\"przewalskii\", \n",
    "    \"33588_przewalskii\":\"przewalskii\",\n",
    "    \"41478_cyathophylloides\":\"cyathophylloides\", \n",
    "    \"41954_cyathophylloides\":\"cyathophylloides\", \n",
    "    \"29154_superba\":\"superba\",\n",
    "    \"30686_cyathophylla\":\"cyathophylla\", \n",
    "    \"33413_thamno\":\"thamno\", \n",
    "    \"30556_thamno\":\"thamno\", \n",
    "    \"35236_rex\":\"rex\", \n",
    "    \"40578_rex\":\"rex\", \n",
    "    \"35855_rex\":\"rex\",\n",
    "    \"39618_rex\":\"rex\", \n",
    "    \"38362_rex\":\"rex\"\n",
    "}"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Copying this notebook to your computer/cluster\n",
    "You can easily copy this notebook and then just replace my file names with your filenames to run your analysis. Just click on the [Download Notebook] link at the top of this page. Then run `jupyter-notebook` from a terminal and open this notebook from the dashboard."
   ]
  }
 ],
 "metadata": {
  "celltoolbar": "Raw Cell Format",
  "kernelspec": {
   "display_name": "Python 2",
   "language": "python",
   "name": "python2"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 2
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython2",
   "version": "2.7.15"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 1
}
