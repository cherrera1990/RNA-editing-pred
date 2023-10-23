import tensorflow as tf
import numpy as np

# =============================================================================================
# D A T A G E N E R A T O R
# =============================================================================================
"""
DataGenerator and DataGeneratorWithLUT are both classes to use during
machine learning training, for minibatch generation.
The generate batches and code integers to One-Hot on the fly.
DataGenerator uses the function tf.keras.layers.CategoryEncoding which is quite
efficient, but available only from tensorflow 2.7.0.

DataGeneratorWithLUT uses a lookup table instead of CategoryEncoding.

"""


class DataGenerator(tf.keras.utils.Sequence):

    def __init__(self, xdata, ylabels, batch_size=32, seq_size=41, categories_size=20, shuffle=True):
        self.xdata = xdata
        self.ylabels = ylabels
        self.batch_size = batch_size
        self.seq_size = seq_size
        self.categories_size = categories_size
        self.indices = self.xdata.index.tolist()
        self.shuffle = shuffle
        self.OneHotConverter = tf.keras.layers.CategoryEncoding(num_tokens=categories_size, output_mode="one_hot")
        self.on_epoch_end()
        # print("len(indices)={}".format(len(self.indices)))

    def __len__(self):
        'Denotes the number of batches per epoch'
        return (len(self.indices) // self.batch_size)

    def __getitem__(self, index):
        'indexes are the 32 indexes that make up the batch'
        indexes = self.index[index * self.batch_size:(index + 1) * self.batch_size]
        batch = [self.indices[k] for k in indexes]

        X, y = self.__get_data(batch)
        return X, y

    def on_epoch_end(self):
        self.index = np.arange(len(self.indices))
        if self.shuffle == True:
            np.random.shuffle(self.index)

    def __get_data(self, batch):
        'Generates data containing batch_size samples'  # X : (n_samples, 41,20 )
        # Initialization
        X = np.empty((self.batch_size, self.seq_size, self.categories_size))
        y = np.empty((self.batch_size), dtype=int)

        for i, id in enumerate(batch):
            'i is 0...32 and id are the real indexes'
            # print("i={} , id={}".format(i,id))
            tensor = tf.constant(self.xdata.loc[id])
            X[i,] = self.OneHotConverter(tensor).numpy()
            y[i] = self.ylabels[id]

        return X, y


# ============================================
# ============================================
class DataGeneratorWithLUT(tf.keras.utils.Sequence):

    def __init__(self, xdata, ylabels, lookuptable,
                 batch_size=32,
                 seq_size=41,
                 categories_size=20,
                 shuffle=True):
        self.xdata = xdata
        self.ylabels = ylabels
        self.batch_size = batch_size
        self.seq_size = seq_size
        self.categories_size = categories_size
        self.indices = self.xdata.index.tolist()
        self.shuffle = shuffle
        self.lookuptable = lookuptable
        self.on_epoch_end()

    def __len__(self):
        'Denotes the number of batches per epoch'
        return (len(self.indices) // self.batch_size)

    def __getitem__(self, index):
        'indexes are the 32 indexes that make up the batch'
        indexes = self.index[index * self.batch_size:(index + 1) * self.batch_size]
        batch = [self.indices[k] for k in indexes]

        X, y = self.__get_data(batch)
        return X, y

    def on_epoch_end(self):
        self.index = np.arange(len(self.indices))
        if self.shuffle == True:
            np.random.shuffle(self.index)

    def __get_data(self, batch):
        'Generates data containing batch_size samples'  # X : (n_samples, 41,20 )
        # Initialization
        X = np.empty((self.batch_size, self.seq_size, self.categories_size))
        y = np.empty((self.batch_size), dtype=int)

        for i, id in enumerate(batch):
            'i is 0...32 and id are the real indexes'
            line = []
            for j in range(0, len(self.xdata.loc[id])):
                line.append(self.lookuptable[self.xdata.loc[id][j]])
            X[i,] = np.array(line)
            y[i] = self.ylabels[id]

        return X, y

# =============================================================================================

