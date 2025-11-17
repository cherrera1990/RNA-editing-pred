# RNA-editing-pred
RNA-editing prediction using RF and biLSTM
This project contains the source code created and used for the study and prediction of RNA-editing using Random Forest and Neural Networks as described in Zawisza et al. 2023.
* The src folder contains the c++ source code for the programs used to create the datasets for both Random Forest and Neural Networks approaches. (for help, contact: m.zawisza@ub.edu)
* Tha bash folder contains bash scripts that mostly automate the first part of the pipeline (prepare_datasets_1-6.sh) and other auxiliary bash scripts. (for help, contact: m.zawisza@ub.edu)
* biLSTM folder contains code for bidirectional LSTM deep learning approach, developed in python with tensorflow.
* RF folder contains code for Random Forests, developed in R.
* The example_datasets contain a few  RNA-editing datasets, including the original REDIportal human and mouse datasets, as well as a smaller version of the human dataset The small dataset should only be used as a test dataset, it is not expected to obtain good results when training a machine learning model.
