# single-cell-dynamics-classifier

This repo contains algorithm I developed for clustering single-cell time-series data at Prof Alexander Loewer's lab. For details, see my [paper](https://doi.org/10.1016/j.celrep.2019.03.031). It clusters single-cell trajactories based on their shapes, and can effectively deal with time-shift and scaling variance issues when comparing time-series signals.

This clustering algorithm was built based on [k-shape clustering](http://www.cs.columbia.edu/~gravano/Papers/2016/sigmod-record16.pdf), which nowadays is widely implemented in time-series analysis. 

The file ‘p21_clustering.m’ will guide you to these algorithms.