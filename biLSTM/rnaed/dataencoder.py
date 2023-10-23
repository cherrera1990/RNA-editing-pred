import multiprocessing as mp
import itertools
import pandas as pd
import numpy as np
import concurrent.futures
import os
import math
from rnaed import exceptions

# =============================================================================================
# ENCODING FUNCTIONS TO BE USED IN PARALLEL PROCESSING
# =============================================================================================
def encode_line(row, lut, nchannels):
    channel_selector = 'ADN'  # 'ADN'(3) or 'SECONDARY'(4)
    if nchannels == 1 and channel_selector == 'SECONDARY':
        try:
            # line=[row['GENE'],row['POS'],row['EDITING']]
            line = [row[0], row[1], row[2]]
            for i in range(0, len(row['SECONDARY'])):
                # tupla=(row['SECONDARY'][i],)
                tupla = (row[4][i],)
                code = lut[tupla]
                line.append(code)
            return line
        except KeyError:
            print("Error undefined char in sequence {} pos:{}".format(row[0], row[1]))
        except:
            print("Error encoding to int!")
    if nchannels == 1 and channel_selector == 'ADN':
        try:
            # line=[row['GENE'],row['POS'],row['EDITING']]
            line = [row[0], row[1], row[2]]
            for i in range(0, len(row['ADN'])):
                # tupla=(row['ADN'][i],)
                tupla = (row[3][i],)
                code = lut[tupla]
                line.append(code)
            return line
        except KeyError:
            print("Error undefined char in sequence {} pos:{}".format(row[0], row[1]))
        except:
            print("Error encoding to int!")
    if nchannels == 2:
        try:
            # line=[row['GENE'],row['POS'],row['EDITING']]
            line = [row[0], row[1], row[2]]
            for i in range(0, len(row['ADN'])):
                # tupla=(row['ADN'][i], row['SECONDARY'][i])
                tupla = (row[3][i], row[4][i])
                code = lut[tupla]
                line.append(code)
            return line
        except KeyError:
            print("Error undefined char in sequence {} pos:{}".format(row[0], row[1]))
        except:
            print("Error encoding to int!")
    if nchannels == 3:
        try:
            # line=[row['GENE'],row['POS'],row['EDITING']]
            line = [row[0], row[1], row[2]]
            for i in range(0, len(row['ADN'])):
                # tupla=(row['ADN'][i], row['SECONDARY'][i], row['PAIRS'][i])
                tupla = (row[3][i], row[4][i], row[5][i])
                code = lut[tupla]
                line.append(code)
            return line
        except KeyError:
            print("Error undefined char in sequence {} pos:{}".format(row[0], row[1]))
        except:
            print("Error encoding to int!")


def encode_to_list(df, lut, nchannels):
    coded_list = df.apply(lambda row: encode_line(row, lut, nchannels), axis=1)
    return (coded_list)


# =============================================================================================
# D A T A   E N C O D I N G
# =============================================================================================

class DataEncoder:

    # class Attributes
    # categories=[['A','G','C','T'],['s','d','h','i','b'],['.', '(', ')']]
    # categories=[['A','G','C','T'],['s','d','h','i','b']]
    # categories=[['A','G','C','T']]

    def __init__(self, categories, padding=True, pad_char="*"):

        self.categories = categories
        self.nchannels = len(categories)
        # LUT rapid look up tables for coding
        self.LUT_Tuple_to_Integer = {}
        self.LUT_Tuple_to_OneHot = {}
        self.LUT_Code_to_OneHot = {}
        self.current = 0
        self.target = 0
        self.num_processors = mp.cpu_count()
        counter = 0
        accu_len = 1
        for category in self.categories:
            accu_len *= len(category)
        if padding == True:
            accu_len += 1
        for element in itertools.product(*self.categories):
            vector = [0] * accu_len
            vector[counter] = 1
            self.LUT_Tuple_to_Integer.update({element: counter})
            self.LUT_Tuple_to_OneHot.update({element: vector})
            self.LUT_Code_to_OneHot.update({counter: vector})
            counter += 1
        if padding == True:
            element = [pad_char] * len(categories)
            element = tuple(element)
            vector = [0] * accu_len
            vector[counter] = 1
            self.LUT_Tuple_to_Integer.update({element: counter})
            self.LUT_Tuple_to_OneHot.update({element: vector})
            self.LUT_Code_to_OneHot.update({counter: vector})

    def __serial_encode(self, df, coded_list, lookuptable):
        ''' Basic funtion to encode sequences in dataframe df
            , using wether LUT_Tuple_to_Integer or LUT_Tuple_to_OneHot '''
        if self.nchannels == 3:
            for ind_row, row in df.iterrows():
                try:
                    # line=[row['GENE'],row['POS'],row['EDITING']]
                    line = [row[0], row[1], row[2]]
                    for i in range(0, len(row['ADN'])):
                        # tupla=(row['ADN'][i], row['SECONDARY'][i], row['PAIRS'][i])
                        tupla = (row[3][i], row[4][i], row[5][i])
                        code = lookuptable[tupla]
                        line.append(code)
                    # print("line[{}]->{}\n".format(type(line),line))
                    coded_list.append(line)
                    self.current = ind_row + 1
                except KeyError:
                    print("Error undefined char in sequence {} : {} pos:{}".format(ind_row + 1, row[0], row[1]))
                    continue
                except:
                    print("Error encoding to int!")
        if self.nchannels == 2:
            for ind_row, row in df.iterrows():
                try:
                    # line=[row['GENE'],row['POS'],row['EDITING']]
                    line = [row[0], row[1], row[2]]
                    for i in range(0, len(row['ADN'])):
                        # tupla=(row['ADN'][i], row['SECONDARY'][i])
                        tupla = (row[3][i], row[4][i])
                        code = lookuptable[tupla]
                        line.append(code)
                    # print("line[{}]->{}\n".format(type(line),line))
                    coded_list.append(line)
                    self.current = ind_row + 1
                except KeyError:
                    print("Error undefined char in sequence {} : {} pos:{}".format(ind_row + 1, row[0], row[1]))
                    continue
                except:
                    print("Error encoding to int!")

    def __parallel_encode(self, df, coded_list, lut):
        splitted_df = np.array_split(df, self.num_processors)
        # Preocess Pool
        with concurrent.futures.ProcessPoolExecutor(max_workers=self.num_processors) as executor:
            futures = [executor.submit(encode_to_list, df, lut, self.nchannels) for df in splitted_df]
            # this line waits for all the process to finish
            for future in concurrent.futures.as_completed(futures):
                try:
                    result_series = future.result()
                    result_list = result_series.tolist()
                    # print("result_list[{}]->{}\n".format(type(result_list),result_list))
                    # very important: to use extend not append
                    coded_list.extend(result_list)
                except Exception as ex:
                    print("Exception during __parallel_encode:" + str(ex))
                    pass

    # ------------------------------------------------------------
    def encodeToFile(self, from_dir, from_file, to_dir, to_file, sep=';',
                     parallel=True, encode_to="integer"):
        ''' parameter encode_to has two posible values [integer,one-hot]'''
        lut = None
        if encode_to == "integer":
            print("\nEncoding DataSet to integers with {} channels".format(str(self.nchannels)))
            #self.print_LUT_TupleToInteger()
            lut = self.LUT_Tuple_to_Integer
        elif encode_to == "one-hot":
            print("\nEncoding DataSet to one-hot...")
            lut = self.LUT_Tuple_to_OneHot
        else:
            print("\nError, encode_to parameter is not valid, forced to integer")
            lut = self.LUT_Tuple_to_Integer

        fromfilename = os.path.join(from_dir, from_file)
        tofilename = os.path.join(to_dir, to_file)

        # Load from file to dataframe
        df = pd.read_csv(fromfilename, sep=sep, header=0)
        # Convert all letters to upper/lower case
        df["ADN"] = df["ADN"].str.upper()
        df["SECONDARY"] = df["SECONDARY"].str.lower()
        # Run along the strings encoding
        coded_list = []
        num_rows = df.shape[0]
        self.target = num_rows

        if parallel == True:
            self.__parallel_encode(df, coded_list, lut)
        else:
            self.__serial_encode(df, coded_list, lut)

        df_coded = pd.DataFrame(coded_list)

        # headers
        columns = [df.columns[0], df.columns[1], df.columns[2]]
        num_columns = df_coded.shape[1]
        num_positions = num_columns - 3

        for i in range(1, num_positions + 1):
            columns.append("P" + str(i))
        df_coded.columns = columns

        # Save dataframe to file
        df_coded.to_csv(tofilename, index=False, sep=sep)

    # ------------------------------------------------------------
    def createBalancedDatasetFile(self, from_file, to_file, pprop=1, ntimes=1, sep=";"):
        """ This funtion load the from_file and balances the positive vs negative samples
            in the following way:
            -It takes pprop [0.0,1.0] proportion out of all positive sammples P
            -then it takes up to (ntimes*pprop*P) negative samples.
        """
        print("\nCreating balanced dataset...")
        try:
            np.random.seed(2022)
            # Load from file to dataframe
            df = pd.read_csv(from_file, sep=sep, header=0)
            # Select all the positive sequences
            df_1 = df.loc[df['EDITING'] == 1]
            # Subset the pprop proportion of random positive ones
            df_1 = df_1.sample(frac=pprop)
            # Final number of selected positives
            npos = df_1.shape[0]

            if npos == 0:
                raise exceptions.PositiveCasesNotFound(from_file)

            # Select all the negative sequences
            df_0 = df.loc[df['EDITING'] == 0]
            # subset ntimes*positves out of random negative ones
            df_0 = df_0.sample(n=math.ceil(ntimes * npos))
            # merge selected positive and negative samples
            df_10 = pd.concat([df_0, df_1])
            # copy to dataset file
            df_10.to_csv(to_file, index=False, sep=sep)
        except FileNotFoundError as e:
            print("FileNotFoundError:", e)
        except ValueError as e:
            print("ValueError:", e)
        except exceptions.PositiveCasesNotFound as e:
            print("PositiveCasesNotFound:", e)
            raise


    # ------------------------------------------------------------
    def extractRandomDatasetFile(self, from_file, to_file, nsamples, target=0, sep=";"):
        """ This funtion load the from_file and extract random samples
            from the target class
        """
        print("\nExtracting dataset...")
        try:
            np.random.seed(2022)
            # Load from file to dataframe
            df = pd.read_csv(from_file, sep=sep, header=0)
            # Select all the target sequences
            df_t = df.loc[df['EDITING'] == target]
            # subset nsamples out of random negative ones
            df_t = df_t.sample(n=nsamples)
            # copy to dataset file
            df_t.to_csv(to_file, index=False, sep=sep)
        except FileNotFoundError as e:
            print("FileNotFoundError:", e)
        except ValueError as e:
            print("ValueError:", e)

    # ------------------------------------------------------------
    def progressRequest(self):
        return self.current, self.target;

    # ------------------------------------------------------------
    def getChannels(self):
        return (len(self.categories), self.categories)

    # ------------------------------------------------------------
    def getNumCodes(self):
        return (len(self.LUT_Tuple_to_Integer))

    # ------------------------------------------------------------
    """ Print Look Up Tables"""

    def print_LUT_TupleToInteger(self):
        print("Number of codes={}".format(len(self.LUT_Tuple_to_Integer)))
        for key, value in self.LUT_Tuple_to_Integer.items():
            print("{} / {}".format(key, value))

    def print_LUT_CodeToOneHot(self):
        print("Number of codes={}".format(len(self.LUT_Code_to_OneHot)))
        for key, value in self.LUT_Code_to_OneHot.items():
            string_ints = [str(i) for i in value]
            print("{} / {}".format(key, "".join(string_ints)))

    def print_LUT_IntegerToTuple(self):
        print("Number of codes={}".format(len(self.LUT_Tuple_to_Integer)))
        for key, value in self.LUT_Tuple_to_Integer.items():
            print("{} -> {}".format(value, key))


    def get_LUT_TupleToInteger(self):
        return self.LUT_Tuple_to_Integer

    def get_LUT_CodeToOneHot(self):
        return self.LUT_Code_to_OneHot

