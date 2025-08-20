# Neural Networks for Pathogen Phylogeography
This repository is maintained by Lucy O'Brien, Cécile Tran Kiem, and Trevor Bedford.

## Overview
Genetic similarities between pathogens can be indicative of historical disease transmission patterns across geographic regions. By comparing similarities across genomes, we can infer the most likely geographic origin of a pathogen. This process of phylogeographic inference is typically done using tree-based or probabilistic methods. However, these methods are computationally inefficient at large scales and can fail to capture the true evolutionary history of a pathogen. Recent research ([Battey et al., 2020](https://elifesciences.org/articles/54507)) indicates that machine learning approaches to phylogeography can reduce computational cost compared to traditional methods, while maintaining or increasing accuracy. 

This project aimed to further explore the promise of machine learning-based phylogeography in the context of categorical inference. For this project, we developed a multi-class classification multilayer perceptron model for pathogen phylogeography. Code and data for this model, along with figures analyzing its performance, are stored in this repository.

## Repository Structure
This repository includes:
- `calibration-dta`: Code and results for discrete trait analysis 
- `neural-network`: Code relating to the architecture, optimization, training, and results of our neural network. 
- `phylogenetic-inference`: Phylogenetic trees for simulated datasets.
- `simulations`: The simulated data used to train and test the neural network. This folder contains both the simulations and the resulting data.
