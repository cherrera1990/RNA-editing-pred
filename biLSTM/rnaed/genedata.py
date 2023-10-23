import os
import re #regular expressions
import pandas as pd
import sys
import numpy as np
import multiprocessing as mp
import concurrent.futures
from rnaed import exceptions
import random


#CONSTANTS
IDX_HEADER=0
IDX_SEQ=1
IDX_PAIRS=2
IDX_SECONDARY=3
IDX_EDITING=4

VALID_SEQ_CHARS='^[aAgGcCtT]+$'
VALID_PAIRS_CHARS='^[.()]+$'
VALID_SECONDARY_CHARS='^[sdhib]+$'
VALID_EDITING_CHARS='^[012]+$'

PATTERN_GEN = re.compile(VALID_SEQ_CHARS);
PATTERN_PAIRS = re.compile(VALID_PAIRS_CHARS);
PATTERN_SECONDARY = re.compile(VALID_SECONDARY_CHARS);
PATTERN_EDITING = re.compile(VALID_EDITING_CHARS);


# ------------------------------------------------------------
def splitSequenceParallel(lista, seq, w_width, n_ch, pattern, padding, padding_char):
    seq_length = len(seq[IDX_SEQ])
    header = seq[IDX_HEADER]
    for nucleoside in re.finditer(pattern, seq[IDX_SEQ]):
        try:
            nuc_pos = nucleoside.start()
            # Spliting indexes
            min_idx = max(nuc_pos - w_width, 0)
            max_idx = min(nuc_pos + w_width, seq_length)
            # len_frag=(max_idx-min_idx)+1
            adn = ''.join(seq[IDX_SEQ][min_idx:max_idx + 1])
            secondary = ''.join(seq[IDX_SECONDARY][min_idx:max_idx + 1])
            # We check that symbols are valid
            if checkValidChars(adn, PATTERN_GEN) == False:
                raise exceptions.InvalidCharInGeneSeq(header, nuc_pos)
            if checkValidChars(secondary, PATTERN_SECONDARY) == False:
                raise exceptions.InvalidCharInSecondarySeq(header, nuc_pos)

            if (n_ch == 3):
                pairs = ''.join(seq[IDX_PAIRS][min_idx:max_idx + 1])
                if checkValidChars(pairs, PATTERN_PAIRS) == False:
                    raise exceptions.InvalidCharInPairsSeq(header, nuc_pos)

            if (min_idx == 0 or max_idx == seq_length):  # Nucleoside on extreme
                if (padding == True):
                    if (min_idx == 0):
                        # Prefix padding
                        len_prefix = nuc_pos
                        padding_string = padding_char * (w_width - len_prefix)
                        adn = padding_string + adn
                        secondary = padding_string + secondary
                        if (n_ch == 3):
                            pairs = padding_string + pairs
                    if (max_idx == seq_length):
                        # Sufix padding
                        len_sufix = max_idx - nuc_pos
                        padding_string = padding_char * (w_width - len_sufix + 1)
                        adn = adn + padding_string
                        secondary = secondary + padding_string
                        if (n_ch == 3):
                            pairs = pairs + padding_string

                else:
                    # Ignore nucleoside when no padding
                    continue

            editing = ""
            if (seq[IDX_EDITING][nuc_pos] == '1'):
                editing = '0'
            else:
                if (seq[IDX_EDITING][nuc_pos] == '2'):
                    editing = '1'
                else:
                    raise exceptions.MissingEditingFlag(header, nuc_pos)

            # ROW
            line = []
            # COLUMN 1
            # gene name without spaces
            gene_name = seq[IDX_HEADER]
            gene_name = gene_name.replace("\t", "_")
            gene_name = re.sub('\s+', '_', gene_name)
            line.append(gene_name)
            # COLUMN 2
            # Position at gene
            line.append(nuc_pos)
            # COLUMN 3
            # EDITING Flag
            line.append(editing)
            # COLUMN 4
            # ADN sequence
            line.append(adn)
            # COLUMN 5
            # SECONDARY sequence
            line.append(secondary)
            # COLUMN 6
            # PAIRS sequence
            if (n_ch == 3):
                line.append(pairs)
            # ADD LINE TO DATASET
            lista.append(line)

        except exceptions.InvalidCharInGeneSeq as e:
            if hasattr(e, 'message'):
                print(e.message)
            else:
                print(e)
            continue
        except exceptions.InvalidCharInPairsSeq as e:
            if hasattr(e, 'message'):
                print(e.message)
            else:
                print(e)
            continue
        except exceptions.InvalidCharInSecondarySeq as e:
            if hasattr(e, 'message'):
                print(e.message)
            else:
                print(e)
            continue
        except exceptions.InvalidCharInEditingSeq as e:
            if hasattr(e, 'message'):
                print(e.message)
            else:
                print(e)
            continue
        except exceptions.MissingEditingFlag as e:
            if hasattr(e, 'message'):
                print(e.message)
            else:
                print(e)
            continue
    # RETORNAMOS LA LISTA
    return lista


# ------------------------------------------------------------
def clipProcess(sequences, w_width, n_ch=2,  # number of channels
                pattern=r'[aA]', padding=True, padding_char='*'):
    num_sequences = len(sequences)
    lista = list()

    for i in range(num_sequences):
        try:
            # Checking equal length in all sequences
            seq = sequences[i]
            len_seq = len(seq[IDX_SEQ])
            # if 3 channels
            if (n_ch == 3 and
                    (len(seq[IDX_PAIRS]) != len_seq)):
                raise exceptions.LenSeqError(seq[IDX_HEADER])
            # if 2 channels
            if ((len(seq[IDX_SECONDARY]) != len_seq) or
                    (len(seq[IDX_EDITING]) != len_seq)):
                raise exceptions.LenSeqError(seq[IDX_HEADER])
            # Checking equal number of adenosines and editing flags
            n_aA = seq[IDX_SEQ].count('a') + seq[IDX_SEQ].count('A')
            n_12 = seq[IDX_EDITING].count('1') + seq[IDX_EDITING].count('2')
            if (n_aA != n_12):
                raise exceptions.NumAdeninesDifferent(seq[IDX_HEADER], n_aA, n_12)
            # Checking equal number of left and right parenthesis
            if (n_ch == 3):
                n_lp = seq[IDX_PAIRS].count('(')
                n_rp = seq[IDX_PAIRS].count(')')
                if (n_lp != n_rp):
                    raise exceptions.NumParenthesisDifferent(i)

            # Clipping sequence to DataSet
            # -----------------------------
            lista = splitSequenceParallel(lista, seq, w_width, n_ch, pattern, padding, padding_char)
            # -----------------------------
        except exceptions.LenSeqError as e:
            if hasattr(e, 'message'):
                print(e.message)
            else:
                print(e)
            continue
        except exceptions.NumAdeninesDifferent as e:
            if hasattr(e, 'message'):
                print(e.message)
            else:
                print(e)
            continue
        except exceptions.NumParenthesisDifferent as e:
            if hasattr(e, 'message'):
                print(e.message)
            else:
                print(e)
            continue

    return lista

#------------------------------------------------------------
def checkValidChars(str, pattern):
	if re.search(pattern, str):
		return True
	else:
		return False

# =============================================================================================
# G E N E   D A T A
# =============================================================================================

class GeneData:
    """ Read Sequences from Genes file"""

    # CLASS ATTRIBUTES
    IDX_HEADER = 0
    IDX_SEQ = 1
    IDX_PAIRS = 2
    IDX_SECONDARY = 3
    IDX_EDITING = 4
    VALID_SEQ_CHARS = '^[aAgGcCtT]+$'
    VALID_PAIRS_CHARS = '^[.()]+$'
    VALID_SECONDARY_CHARS = '^[sdhib]+$'
    VALID_EDITING_CHARS = '^[012]+$'

    # ------------------------------------------------------------
    def __init__(self, data_dir, file_name, encoding='utf-8', line_break='\n', lines_per_seq=5):
        self.fname = os.path.join(data_dir, file_name)
        self.file = open(self.fname, encoding=encoding)
        data = self.file.read()
        self.file.close()
        self.sequences = list()
        self.dataset = list()
        self.num_processors = mp.cpu_count()
        lines = data.split(line_break)
        counter = 0
        sequence = []
        for line in lines:
            if (len(line) != 0):
                sequence.append(line)
                if ((counter % 5) == (lines_per_seq - 1)):
                    # trim pairs sequence from last numbers annotated
                    sequence[self.IDX_PAIRS] = sequence[self.IDX_PAIRS][:len(sequence[self.IDX_SEQ])]
                    self.sequences.append(tuple(sequence))
                    sequence = []
                counter += 1
        print("\nFILE: {}  LENGHT: {} lines / {} genes".format(file_name, len(lines), len(self.sequences)))
        self.pattern_gen = re.compile(self.VALID_SEQ_CHARS);
        self.pattern_pairs = re.compile(self.VALID_PAIRS_CHARS);
        self.pattern_secondary = re.compile(self.VALID_SECONDARY_CHARS);
        self.pattern_editing = re.compile(self.VALID_EDITING_CHARS);

    # ------------------------------------------------------------
    def selectRandomSample(self,ngenes):
        try:
            if ngenes <= len(self.sequences):
                self.sequences = random.sample(self.sequences, ngenes)
            else:
                raise exceptions.NotEnoughGenes(ngenes)
        except exceptions.NotEnoughGenes as e:
            if hasattr(e, 'message'):
                print(e.message)
            else:
                print(e)
            raise exceptions.NotEnoughGenes(ngenes)
    # ------------------------------------------------------------
    def _splitSequence(self, seq, w_width, n_ch, pattern, padding, padding_char):
        seq_length = len(seq[self.IDX_SEQ])
        header = seq[self.IDX_HEADER]
        for nucleoside in re.finditer(pattern, seq[self.IDX_SEQ]):
            try:
                nuc_pos = nucleoside.start()
                # Spliting indexes
                min_idx = max(nuc_pos - w_width, 0)
                max_idx = min(nuc_pos + w_width, seq_length)
                # len_frag=(max_idx-min_idx)+1
                adn = ''.join(seq[self.IDX_SEQ][min_idx:max_idx + 1])
                secondary = ''.join(seq[self.IDX_SECONDARY][min_idx:max_idx + 1])
                # We check that symbols are valid
                if self.checkValidChars(adn, self.pattern_gen) == False:
                    raise exceptions.InvalidCharInGeneSeq(header, nuc_pos)
                if self.checkValidChars(secondary, self.pattern_secondary) == False:
                    raise exceptions.InvalidCharInSecondarySeq(header, nuc_pos)

                if (n_ch == 3):
                    pairs = ''.join(seq[self.IDX_PAIRS][min_idx:max_idx + 1])
                    if self.checkValidChars(pairs, self.pattern_pairs) == False:
                        raise exceptions.InvalidCharInPairsSeq(header, nuc_pos)

                if (min_idx == 0 or max_idx == seq_length):  # Nucleoside on extreme
                    if (padding == True):
                        if (min_idx == 0):
                            # Prefix padding
                            len_prefix = nuc_pos
                            padding_string = padding_char * (w_width - len_prefix)
                            adn = padding_string + adn
                            secondary = padding_string + secondary
                            if (n_ch == 3):
                                pairs = padding_string + pairs
                        if (max_idx == seq_length):
                            # Sufix padding
                            len_sufix = max_idx - nuc_pos
                            padding_string = padding_char * (w_width - len_sufix + 1)
                            adn = adn + padding_string
                            secondary = secondary + padding_string
                            if (n_ch == 3):
                                pairs = pairs + padding_string

                    else:
                        # Ignore nucleoside when no padding
                        continue

                editing = ""
                if (seq[self.IDX_EDITING][nuc_pos] == '1'):
                    editing = '0'
                else:
                    if (seq[self.IDX_EDITING][nuc_pos] == '2'):
                        editing = '1'
                    else:
                        raise exceptions.MissingEditingFlag(header, nuc_pos)

                # ROW
                line = []
                # COLUMN 1
                # gene name without spaces
                gene_name = seq[self.IDX_HEADER]
                gene_name = gene_name.replace("\t", "_")
                line.append(gene_name)
                # COLUMN 2
                # Position at gene
                line.append(nuc_pos)
                # COLUMN 3
                # EDITING Flag
                line.append(editing)
                # COLUMN 4
                # ADN sequence
                line.append(adn)
                # COLUMN 5
                # SECONDARY sequence
                line.append(secondary)
                # COLUMN 6
                # PAIRS sequence
                if (n_ch == 3):
                    line.append(pairs)
                # ADD LINE TO DATASET
                self.dataset.append(line)

            except exceptions.InvalidCharInGeneSeq as e:
                if hasattr(e, 'message'):
                    print(e.message)
                else:
                    print(e)
                continue
            except exceptions.InvalidCharInPairsSeq as e:
                if hasattr(e, 'message'):
                    print(e.message)
                else:
                    print(e)
                continue
            except exceptions.InvalidCharInSecondarySeq as e:
                if hasattr(e, 'message'):
                    print(e.message)
                else:
                    print(e)
                continue
            except exceptions.InvalidCharInEditingSeq as e:
                if hasattr(e, 'message'):
                    print(e.message)
                else:
                    print(e)
                continue
            except exceptions.MissingEditingFlag as e:
                if hasattr(e, 'message'):
                    print(e.message)
                else:
                    print(e)
                continue

                # ------------------------------------------------------------

    def checkValidChars(self, str, pattern):
        if re.search(pattern, str):
            return True
        else:
            return False
            # ------------------------------------------------------------

    def clipCenteredSequencesToFile(self,
                                    data_dir,
                                    file_name,
                                    w_width,
                                    n_ch=2,  # number of channels
                                    pattern=r'[aA]',
                                    padding=True,
                                    padding_char='*',
                                    separator=';'):

        print("\nProcessing {} genes in windows W+1+W W={}...".format(len(self.sequences), w_width))
        num_sequences = len(self.sequences)
        for i in range(num_sequences):
            try:
                # Checking equal length in all sequences
                seq = self.sequences[i]
                len_seq = len(seq[self.IDX_SEQ])
                # if 3 channels
                if (n_ch == 3 and
                        (len(seq[self.IDX_PAIRS]) != len_seq)):
                    raise exceptions.LenSeqError(seq[self.IDX_HEADER])
                # if 2 channels
                if ((len(seq[self.IDX_SECONDARY]) != len_seq) or
                        (len(seq[self.IDX_EDITING]) != len_seq)):
                    raise exceptions.LenSeqError(seq[self.IDX_HEADER])
                # Checking equal number of adenosines and editing flags
                n_aA = seq[self.IDX_SEQ].count('a') + seq[self.IDX_SEQ].count('A')
                n_12 = seq[self.IDX_EDITING].count('1') + seq[self.IDX_EDITING].count('2')
                if (n_aA != n_12):
                    raise exceptions.NumAdeninesDifferent(seq[self.IDX_HEADER], n_aA, n_12)
                # Checking equal number of left and right parenthesis
                if (n_ch == 3):
                    n_lp = seq[self.IDX_PAIRS].count('(')
                    n_rp = seq[self.IDX_PAIRS].count(')')
                    if (n_lp != n_rp):
                        raise exceptions.NumParenthesisDifferent(i)

                # Reporting activity on console
                # -----------------------------
                sys.stdout.write("\r Clipping gene: {} of {} \r".format(i + 1, num_sequences))
                # print("Secuencia: {} de {}".format(i+1,num_sequences), flush=True)
                # sys.stdout.flush()

                # Clipping sequence to DataSet
                # -----------------------------
                self._splitSequence(seq, w_width, n_ch, pattern, padding, padding_char)
                # -----------------------------
            except exceptions.LenSeqError as e:
                if hasattr(e, 'message'):
                    print(e.message)
                else:
                    print(e)
                continue
            except exceptions.NumAdeninesDifferent as e:
                if hasattr(e, 'message'):
                    print(e.message)
                else:
                    print(e)
                continue
            except exceptions.NumParenthesisDifferent as e:
                if hasattr(e, 'message'):
                    print(e.message)
                else:
                    print(e)
                continue

            # except:
            # print("Error: Something else went wrong when saving sequences!")
            # break
        # saving dataset to dataframe
        # ---------------------------
        if (n_ch == 3):
            columns = ['GENE', 'POS', 'EDITING', 'ADN', 'SECONDARY', 'PAIRS']
        else:
            columns = ['GENE', 'POS', 'EDITING', 'ADN', 'SECONDARY']
        df = pd.DataFrame(self.dataset, columns=columns)

        # saving dataframe to CSV
        # --------------------------
        fname = os.path.join(data_dir, file_name)
        df.to_csv(fname, index=False, sep=separator)

    # ------------------------------------------------------------
    def clipCenteredSequencesToFileParallel(self,
                                            data_dir,
                                            file_name,
                                            w_width,
                                            n_ch=2,  # number of channels
                                            pattern=r'[aA]',
                                            padding=True,
                                            padding_char='*',
                                            separator=';'):

        print("\nProcessing {} genes in windows W+1+W W={}...".format(len(self.sequences), w_width))
        clipped_sequence_list = list()
        # split sequences list into num_processors list
        splitted_list_sequences = np.array_split(self.sequences, self.num_processors)
        del self.sequences
        # Process Pool
        with concurrent.futures.ProcessPoolExecutor(max_workers=self.num_processors) as executor:
            futures = [executor.submit(clipProcess, list_seq, w_width, n_ch, pattern, padding, padding_char) for
                       list_seq in splitted_list_sequences]
            # this line waits for all the process to finish
            for future in concurrent.futures.as_completed(futures):
                try:
                    result_series = future.result()
                    # result_list=result_series
                    # print("result_list[{}]->{}\n".format(type(result_list),result_list))
                    # very important: to use extend not append
                    clipped_sequence_list.extend(result_series)
                except Exception as ex:
                    print("Exception during parallel clipping:" + str(ex))
                    pass

        # saving dataset to dataframe
        # ---------------------------
        if (n_ch == 3):
            columns = ['GENE', 'POS', 'EDITING', 'ADN', 'SECONDARY', 'PAIRS']
        else:
            columns = ['GENE', 'POS', 'EDITING', 'ADN', 'SECONDARY']
        df = pd.DataFrame(clipped_sequence_list, columns=columns)

        # saving dataframe to CSV
        # --------------------------
        fname = os.path.join(data_dir, file_name)
        df.to_csv(fname, index=False, sep=separator)

    # ------------------------------------------------------------
    def checkDataConsistency(self):

        lengthOK = True
        adenineNumOk = True
        parenthesisOK = True
        charsInGeneOK = True
        charsInPairsOK = True
        charsInSecondaryOK = True
        charsInEditingOK = True

        print("\nChecking data consistency in {} sequences...".format(len(self.sequences)))

        for i in range(len(self.sequences)):
            try:
                # Checking equal length in all sequences
                seq = self.sequences[i]
                len_seq = len(seq[self.IDX_SEQ])
                if ((len(seq[self.IDX_PAIRS]) != len_seq) or
                        (len(seq[self.IDX_SECONDARY]) != len_seq) or
                        (len(seq[self.IDX_EDITING]) != len_seq)):
                    raise exceptions.LenSeqError(seq[self.IDX_HEADER])
                # Checking equal number of adenines and editing flags
                n_aA = seq[self.IDX_SEQ].count('a') + seq[self.IDX_SEQ].count('A')
                n_12 = seq[self.IDX_EDITING].count('1') + seq[self.IDX_EDITING].count('2')
                if (n_aA != n_12):
                    raise exceptions.NumAdeninesDifferent(seq[self.IDX_HEADER], n_aA, n_12)
                # Checking equal number of left and right parenthesis
                n_lp = seq[self.IDX_PAIRS].count('(')
                n_rp = seq[self.IDX_PAIRS].count(')')
                if (n_lp != n_rp):
                    raise exceptions.NumParenthesisDifferent(i)

                # Checking if any character is gene sequence is not allowed
                if self.checkValidChars(seq[self.IDX_SEQ], self.pattern_gen) == False:
                    raise exceptions.InvalidCharInGeneSeq(i)
                # Checking if any character is pairs sequence is not allowed
                if self.checkValidChars(seq[self.IDX_PAIRS], self.pattern_pairs) == False:
                    raise exceptions.InvalidCharInPairsSeq(i)
                # Checking if any character is secondary sequence is not allowed
                if self.checkValidChars(seq[self.IDX_SECONDARY], self.pattern_secondary) == False:
                    raise exceptions.InvalidCharInSecondarySeq(i)
                # Checking if any character is editing sequence is not allowed
                if self.checkValidChars(seq[self.IDX_EDITING], self.pattern_editing) == False:
                    raise exceptions.InvalidCharInEditingSeq(i)

            except exceptions.LenSeqError as e:
                lengthOK = False
                if hasattr(e, 'message'):
                    print(e.message)
                else:
                    print(e)
            except exceptions.NumAdeninesDifferent as e:
                adenineNumOk = False
                if hasattr(e, 'message'):
                    print(e.message)
                else:
                    print(e)
            except exceptions.NumParenthesisDifferent as e:
                parenthesisOK = False
                if hasattr(e, 'message'):
                    print(e.message)
                else:
                    print(e)
            except exceptions.InvalidCharInGeneSeq as e:
                charsInGeneOK = False
                if hasattr(e, 'message'):
                    print(e.message)
                else:
                    print(e)
            except exceptions.InvalidCharInPairsSeq as e:
                charsInPairsOK = False
                if hasattr(e, 'message'):
                    print(e.message)
                else:
                    print(e)
            except exceptions.InvalidCharInSecondarySeq as e:
                charsInSecondaryOK = False
                if hasattr(e, 'message'):
                    print(e.message)
                else:
                    print(e)
            except exceptions.InvalidCharInEditingSeq as e:
                charsInEditingOK = False
                if hasattr(e, 'message'):
                    print(e.message)
                else:
                    print(e)
            except:
                print("Something else went wrong!")
            else:
                continue

        print("All sequences are consistent in length... {}".format(lengthOK))
        print("All sequences are consistent in adenine number... {}".format(adenineNumOk))
        print("All sequences are consistent in parenthesis left/right... {}".format(parenthesisOK))
        print("All sequences are consistent in valid gene chars... {}".format(charsInGeneOK))
        print("All sequences are consistent in valid pairs chars... {}".format(charsInPairsOK))
        print("All sequences are consistent in valid secondary chars... {}".format(charsInSecondaryOK))
        print("All sequences are consistent in valid editing chars... {}".format(charsInEditingOK))

    # ------------------------------------------------------------
    """Print sequence n"""

    def printSeq(self, n):
        seq = self.sequences[n]
        print("\nSEQUENCE: {}          \n=============".format(n))
        print("\nHEADER:              \n----------\n{}".format(seq[self.IDX_HEADER]))
        print("\nSEQUENCE: [{} chars] \n----------\n{}".format(len(seq[self.IDX_SEQ]), seq[self.IDX_SEQ]))
        print("\nPAIRS: [{} chars]    \n----------\n{}".format(len(seq[self.IDX_PAIRS]), seq[self.IDX_PAIRS]))
        print("\nSECONDARY: [{} chars]\n----------\n{}".format(len(seq[self.IDX_SECONDARY]), seq[self.IDX_SECONDARY]))
        print("\nADENINE: [{} chars]  \n----------\n{}".format(len(seq[self.IDX_EDITING]), seq[self.IDX_EDITING]))
        print("\nnumber of adenines: [{}]".format(seq[self.IDX_SEQ].count('a') + seq[self.IDX_SEQ].count('A')))
        print("\nnumber of 1,2 flags: [{}]".format(seq[self.IDX_EDITING].count('1') + seq[self.IDX_EDITING].count('2')))
        print("\nnumber of left parenthesis: [{}]".format(seq[self.IDX_PAIRS].count('(')))
        print("\nnumber of right parenthesis: [{}]".format(seq[self.IDX_PAIRS].count(')')))
        # positions with Adenosine
        adenosinePattern = r'[aA]'
        for aA in re.finditer(adenosinePattern, seq[self.IDX_SEQ]):
            print('aA found at pos:', aA.start())

    # ------------------------------------------------------------
    def getSeqList(self):
        return (self.sequences)

