# Notebooks

The notebooks used in this study are:

## TFM_LSTMBidiAttention.ipynb

WINDOW_SIZE: 50+1+50
TRAIN: human x_train:(214803, 101) y_train:(214803,) from RANDOM1000
VALIDATE: human x_val: (46029, 101) y_val: (46029,) from RANDOM1000
TEST: human x_test: (46030, 101) y_test: (46030,) from RANDOM1000
MODEL: LSTM Biderectional with custom attention layer
MODEL NAME: TFM_LSTMBidiAttention
ACCURACY: 0.95
KAPPA: 0.896
AUC: 0.97

## TFM_LSTMBidiSelfAttention.ipynb

WINDOW_SIZE: 50+1+50
TRAIN: human x_train:(214803, 101) y_train:(214803,) from RANDOM1000
VALIDATE: human x_val: (46029, 101) y_val: (46029,) from RANDOM1000
TEST: human x_test: (46030, 101) y_test: (46030,) from RANDOM1000
MODEL: LSTM Biderectional with standar self-attention layer
MODEL NAME: TFM_LSTMBidiSelfAttention
ACCURACY: 0.95
KAPPA: 0.896
AUC: 0.97

## TFM_LSTMBidiAttention_1CH_ADN.ipynb

WINDOW_SIZE: 50+1+50
TRAIN: human x_train:(214803, 101) y_train:(214803,) from RANDOM1000
VALIDATE: human x_val: (46029, 101) y_val: (46029,) from RANDOM1000
TEST: human x_test: (46030, 101) y_test: (46030,) from RANDOM1000
MODEL: LSTM Biderectional with custom attention layer
MODEL NAME: TFM_LSTMBidiAttention_1CH_ADN
ACCURACY: 0.95
KAPPA: 0.894
AUC: 0.97

## TFM_LSTMBidiAttention_1CH_STR.ipynb

WINDOW_SIZE: 50+1+50
TRAIN: human x_train:(214803, 101) y_train:(214803,) from RANDOM1000
VALIDATE: human x_val: (46029, 101) y_val: (46029,) from RANDOM1000
TEST: human x_test: (46030, 101) y_test: (46030,) from RANDOM1000
MODEL: LSTM Biderectional with custom attention layer
MODEL NAME: TFM_LSTMBidiAttention_1CH_STR
ACCURACY: 0.86
KAPPA: 0.721
AUC: 0.92

## TFM_Prediction_LSTMBidiAttention.ipynb.ipynb

WINDOW_SIZE: 50+1+50
TRAIN: human x_train:(214803, 101) y_train:(214803,) from RANDOM1000
VALIDATE: human x_val: (46029, 101) y_val: (46029,) from RANDOM1000
TEST: human x_test.shape=(1000000, 101) y_test.shape=(1000000,) TEST LABELS: 0 500225 1  499775  FROM COMPLETE GENOME
MODEL NAME:TFM_LSTMBidiAttention
ACCURACY: 0.91
KAPPA: 0.828
AUC: 0.94

## TFM_Prediction10Genes_LSTMBidiAttention.ipynb

WINDOW_SIZE: 50+1+50
TRAIN: human x_train:(214803, 101) y_train:(214803,) from RANDOM1000
VALIDATE: human x_val: (46029, 101) y_val: (46029,) from RANDOM1000
TEST: human x_test.shape=(17005, 101) y_test.shape=(17005,) TEST LABELS:0 16803 1 202, ONLY 10 COMPLETE GENES
MODEL NAME:TFM_LSTMBidiAttention
ACCURACY: 0.94
KAPPA: 0.281
AUC: NULL

## TFM_LSTMBidiAttention_400K.ipynb

WINDOW_SIZE: 50+1+50
TRAIN: human x_train.shape=(300000, 101) y_train.shape=(300000,) FROM COMPLETE GENOME
VALIDATE:human x_val.shape=(100000, 101) y_val.shape=(100000,) FROM COMPLETE GENOME
TEST: human x_test.shape=(400000, 101) y_test.shape=(400000,) FROM COMPLETE GENOME
MODEL NAME: TFM_LSTMBidiAttention_400K
ACCURACY: 0.93
KAPPA: 0.854
AUC: 0.95

## TFM_LSTMBidiAttention_2CH_MM.ipynb

WINDOW_SIZE: 50+1+50
TRAIN: mouse x_train:(41505, 101) y_train:(41505,)
VALIDATE: mouse x_val: (8894, 101) y_val: (8894,)
TEST: mouse x_test: (8895, 101) y_test: (8895,) TEST LABELS: 0 4509 1  4386
MODEL NAME: TFM_LSTMBidiAttention_MM_ALLMM
ACCURACY: 0.84
KAPPA: 0.671
AUC: 0.89

## TFM_LSTMBidiAttentionW100.ipynb 

WINDOW_SIZE:100+1+100
TRAIN: human x_train:(214803, 101) y_train:(214803,) from RANDOM1000
VALIDATE: human x_val: (46029, 101) y_val: (46029,) from RANDOM1000
TEST: human x_test: (46030, 101) y_test: (46030,) from RANDOM1000
MODEL: LSTM Biderectional with custom attention layer
MODEL NAME: TFM_LSTMBidiAttentionW100
ACCURACY: 0.95
KAPPA: 0.897
AUC: 0.97


## RNNBidiAttention9_Predictor_hg38_to_mm.ipynb

WINDOW_SIZE: 50+1+50
TRAIN: human x_train:(205597, 101) y_train:(205597,)
VALIDATE: human x_val: (101265, 101) y_val: (101265,)
TEST: mouse x_test.shape=(53365, 101) y_test.shape=(53365,)
MODEL NAME: LSTMBidiAttention9
ACCURACY: 0.63
KAPPA: 0.251


## TFM_LSTMBidiAttention_2CH_MM2.ipynb

WINDOW_SIZE: 50+1+50
TRAIN:  Mouse x_train:(41505, 101) y_train:(41505,)
VALIDATE: Mouse x_val: (8894, 101) y_val: (8894,)
TEST: Mouse x_test: (8895, 101) y_test: (8895,)
MODEL NAME: TFM_LSTMBidiAttention_MM_ALLMM2
ACCURACY: 0.83
KAPPA: 0.663
AUC: 0.89

## TFM_CNN_2CH_Hs.ipynb

WINDOW_SIZE: 50+1+50
TRAIN: Human x_train:(378190, 101) y_train:(378190,)
VALIDATE: Human x_val: (81041, 101) y_val: (81041,)
TEST: Human x_test: (81041, 101) y_test: (81041,)
MODEL NAME:TFM_CNN_Hs
ACCURACY: 0.93
KAPPA: 0.85
AUC: 0.95

## TFM_CNN_2CH_Mm.ipynb

WINDOW_SIZE: 50+1+50
DATA: mm10_DL_Pad_W50_BALANCED_1_1_CODED2CH.csv
TRAIN: Mouse x_train:(41505, 101) y_train:(41505,)
VALIDATE: Mouse x_val: (8894, 101) y_val: (8894,)
TEST: Mouse x_test: (8895, 101) y_test: (8895,)
MODEL NAME: TFM_CNN_MM_ALLMM (Convolucional)
ACCURACY: 0.85
KAPPA: 0.702
AUC: 0.91


## TFM_LSTMBidiAttention_400K_hg38Refined.ipynb

Model trained and validated with hg38 refined, without too long genes
WINDOW_SIZE: 50+1+50
DATA: hg38_RS_DL_ROBUST_Pad_W50_BALANCED_1_1_CODED2CH.csv 
TRAIN: human x_train.shape=(300000, 101) y_train.shape=(300000,)
VALIDATE: human x_val.shape=(100000, 101) y_val.shape=(100000,)
TEST: human x_test.shape=(400000, 101) y_test.shape=(400000,)
MODEL NAME: TFM_LSTMBidiAttention_400K_hg38Refined
ACCURACY: 0.93
KAPPA: 0.854
AUC: 0.95

## RNNBidiAttention9_Predictor_mm_to_hg38.ipynb

PREDICTOR FOR HUMAN SAMPLES WITH MOUSE MODEL
WINDOW_SIZE: 50+1+50
DATA: hg38_RS_DL_ROBUST_Pad_W50_BALANCED_1_1_CODED2CH_SMALL.csv
TRAIN: mouse x_train:(41505, 101) y_train:(41505,)
VALIDATE: mouse x_val: (8894, 101) y_val: (8894,)
TEST: human x_test.shape=(79824, 101) y_test.shape=(79824,)
MODEL NAME: TFM_LSTMBidiAttention_MM_ALLMM
ACCURACY: 0.76
KAPPA: 0.513



## TFM_CNN_HYBRID_HsMm.ipynb

Hybrid Model based on Convolutional Neural Network, where outter layers have been trained with Homo sapiens and inner layer with Mus musculus. Final Predictions are made in Mus musculus.
WINDOW_SIZE: 50+1+50
DATA: mm10_DL_Pad_W50_BALANCED_1_1_CODED2CH.csv

LAYERS TRANSFERRED: Human TFM_CNN_Hs Conv1, Conv2
LAYERS RE-TRAINED:  Human TFM_CNN_Hs Conv3
TRAIN: mouse x_train:(41505, 101) y_train:(41505,)
VALIDATE:  mouse x_val: (8894, 101) y_val: (8894,)
TEST:  mouse x_test: (8895, 101) y_test: (8895,)
MODEL NAME: TFM_CNN_HYBRID_HsMm
ACCURACY: 0.76
KAPPA: 0.515
AUC: 0.83


## TFM_CNN_HYBRID_HsMm2.ipynb

Hybrid Model based on Convolutional Neural Network, where outter layer have been trained with Homo sapiens and two inner layers with Mus musculus. Final Predictions are made in Mus musculus.
WINDOW_SIZE: 50+1+50
DATA: mm10_DL_Pad_W50_BALANCED_1_1_CODED2CH.csv
LAYERS TRANSFERRED: Human TFM_CNN_Hs Conv1
LAYERS RE-TRAINED:  Human TFM_CNN_Hs Conv2, Conv3
TRAIN: mouse x_train:(41505, 101) y_train:(41505,)
VALIDATE:  mouse x_val: (8894, 101) y_val: (8894,)
TEST:  mouse x_test: (8895, 101) y_test: (8895,)
MODEL NAME: TFM_CNN_HYBRID_HsMm2
ACCURACY: 0.85
KAPPA: 0.697
AUC: 0.91

## TFM_CNN_2CH_Mm_ActivMap.ipynb

Convolutional Neural Network with Mus musculus including activation map inspection
WINDOW_SIZE: 50+1+50
DATA: mm10_DL_Pad_W50_BALANCED_1_1_CODED2CH.csv
TRAIN: mouse x_train:(41505, 101) y_train:(41505,)
VALIDATE:  mouse x_val: (8894, 101) y_val: (8894,)
TEST:  mouse x_test: (8895, 101) y_test: (8895,)
MODEL NAME: TFM_CNN_MM
ACCURACY: 0.85
KAPPA: 0.692
AUC: 0.91