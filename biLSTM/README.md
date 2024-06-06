To install the package with conda do:
pip install -e .
(You have to be inside the same directory than setup.py)

DEMO USING DEEP LEARNING ALGORITHMS
1.	Edit the script    ClipAndCodeGenes.py  and make sure these parameters are set in section INIT:
    seed_number = 12101492 #Change this number if you don't need to replicate the same dataset
    do_splitting = True  #Execute splitting of genes
    parallel_splitting = False #Execute the splitting on windows in parallel
    parallel_coding = True #Execute the coding in parallel
    do_balance = True #Make the datasets balanced in EDITING/NON EDITING samples
    data_dir = '../../data/raw' #Directory where genes with Editing annotations are
    data_result_dir = '../../data/coded'#Directory where the generated datasets will be copied
    pprop = 1 #proportion of positives we take (Adenosines with Editing) [0, 0.1, 0.1,...,0.9,1]
    ntimes = 1 #number of times the negative ones are in proportion to positive ones (in this case 1:1 balanced)
    files = ['hg38_RS_DL']
    #files = ['dummy']
    # window_widths=[20,50,100]
    window_widths = [50] #List of different window sizes
    categories = [['A', 'G', 'C', 'T'], ['s', 'd', 'h', 'i', 'b']] #Two channels: ADN and Secondary Structure
    #categories=[['s','d','h','i','b']] #One channel: Secondary Structure
    #categories=[['A','G','C','T']] #One channel: AND

    padding = True #use padding to reach window width when A is near the end
    select_genes =True #Select a random number of genes 
    if select_genes:
        number_of_random_genes=1000 #number of genes

2.	Open a Python console and execute script ClipAndCodeGenes.py , after execution, you will have the dataset ready to be used for Learning or Prediction. Just be sure to move the dataset to the directory where  Jupyter Notebook can read it. In this demo, dataset have been generated in '../../data/coded'
3.	Open JupyterLab, Jupyter Notebook or Google Colab and load  Demo.ipynb, in this file you only have to configure this parameters, although they are set ready to work with Demo:
########################################
# PARAMS
########################################
fname      = 'Demo.csv' #dataset with coded windows
separator_char=';' #char used to separate fields in CVS
categories=[['A','G','C','T'],['s','d','h','i','b']]
padding=True
ptrain     = 0.7  #proportion of samples used for training
trace_level = 1
model_name = 'Demo_model'
num_lstm_units=256
n_epochs=40
batch_size=32
categories_size=20
random_seed=2022
########################################
If you are executing the notebook from a PC, with a GPU, you’ll need these directories:
    data_dir = '../data/datasets/Demo'
    models_dir = './SAVED_MODELS'

If you are using Google Colab with data in Google Drive you’ll need these ones:
MyDrive/DATASETS/Demo
MyDrive/MODELS'

You can modify the paths as you want as long as you redefine these directories in the note book.

4.	Finally, execute the notebook step by step.  First, it will load the dataset, then it will split the data into train, validation and test, and finally it will make the training and the testing, showing the results numerically and with a confusion matrix.
