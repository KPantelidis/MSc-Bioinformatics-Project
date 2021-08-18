# MSc-Bioinformatics-Project
## Automated pipeline for isotopic charge states in Mass Spectrometry data

This project is my thesis proposal during my MSc Bioinformatics course at Birkbeck University of London

Understanding the interactions, dynamics and structure of the proteins provides a greater level of comprehension of how a protein works, allowing us to create hypotheses about how to affect it, control it, or modify it.

Mass Spectrometry has been a great analytical technique to study the proteome in the past years and provides a lot of information and data via Mass Spectrum results.

This project focuses on the charge states identification in Mass Spectrum fingerprints and the distinction between true peaks and background noise through an automatic pipeline development instead of a visual inspection. 
The algorithm currently searches for +2, +3 and +4 charge states and the results provide the user with the number of isotopic windows for each charge state in the selected dataset as well as creates .txt files that can be used as labelled training datasets for the creation of a deep learning artificial network.

#### Data
The data files used were provided by Dr. Thalassinos, in .txt file format. The data are a subset of an unpublished project and have been uploaded in PRIDE and held privately at the moment.

The files contain Mass Spectra results in a Tabular form, composed of thousands of data points, corresponding to different m/z and their intensities. The Mass Spectra results are peptides or parts of peptides of different proteins and they were used for code testing and development.
The data extraction, manipulation, analysis and visualisation as well as result analysis were all made using Python programming language, Python version 3.7.4.

#### Documentation
You may find a Jupyter Notebook version of the project as well as a .py file containing all the necessary code.
There is also a .txt file one may use as an example to test and run the code and is the one that I used to create the Jupyter Notebook.
The .pdf file is my thesis proposal as this was submitted on the 25/08/21 at Birkbeck University of London
