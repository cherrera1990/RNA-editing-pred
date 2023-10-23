import rnaed as rna
import os
import time
import datetime
import multiprocessing as mp
import sys
import traceback
import logging
'''
files = ['hg38_RS_DL_ES_first200',
          'hg38_RS_DL_ES_first1000',
          'hg38_RS_DL_ES_NES_rand200',
          'hg38_RS_DL_ES_NES_rand1000',
          'hg38_RS_DL',
          'hg38_RS_DL_ROBUST',
          'hg38_RS_DL_old',
          'hg38_RS_NES_DL',
          'hg38_10_GENES_COMPLETOS',
          'mm10_DL_ES_NES_rand1000',
          'mm10_DL_ES_rand1000',
          'mm10_DL',
          'ttr_DL',
          'ttr_NES_DL',
          'dummy']
'''

def main():
    # ------ INIT ---------------------------
    parallel_spliting = False
    parallel_coding = True
    do_spliting = True
    do_balance = True
    data_dir = '../../data/raw'
    data_result_dir = '../../data/coded'
    pprop = 1 #proportion of positives
    ntimes = 1 #times the positive ones
    files = ['hg38_RS_DL_PRE']
    # window_widths=[20,50,100]
    window_widths = [50]
    categories = [['A', 'G', 'C', 'T'], ['s', 'd', 'h', 'i', 'b']]
    #categories=[['s','d','h','i','b']]
    #categories=[['A','G','C','T']]
    padding = True
    select_genes =True
    if select_genes:
        number_of_random_genes=1000

    # =========================================================
    # M U L T I P R O C E S S I N G
    # =========================================================

    print("Number of processors: {}".format(mp.cpu_count()))

    start = time.perf_counter()
    encoding = rna.DataEncoder(categories, padding=padding, pad_char="*")
    n_ch, cat = encoding.getChannels()
    print("Number of channels: {} /categories: {}".format(n_ch, cat))
    encoding.print_LUT_TupleToInteger()

    for file in files:
        for width in window_widths:
            # Sequences
            # -------------------------------
            seq = rna.GeneData(data_dir, file)
            # seq.checkDataConsistency()
            fname1 = file + '_Pad_W' + str(width)
            start_clipping = time.perf_counter()
            if select_genes:
                try:
                    seq.selectRandomSample(number_of_random_genes)
                    fname1 = fname1 + '_' + str(number_of_random_genes) + "GENES"
                except Exception as e:
                    sys.exit(1)
            if do_spliting == True:
                print("Spliting genes...")
                if parallel_spliting == True:
                    seq.clipCenteredSequencesToFileParallel(data_dir=data_result_dir,
                                                            file_name=fname1 + ".csv",
                                                            w_width=width,
                                                            n_ch=n_ch,
                                                            padding=padding,
                                                            padding_char='*')
                else:
                    seq.clipCenteredSequencesToFile(data_dir=data_result_dir,
                                                    file_name=fname1 + ".csv",
                                                    w_width=width,
                                                    n_ch=n_ch,
                                                    padding=padding,
                                                    padding_char='*')
            else:
                print("Genes are not split in this execution!")
            end_clipping = time.perf_counter()
            elapsed_seconds_clipping = end_clipping - start_clipping
            conversion_clipping = datetime.timedelta(seconds=elapsed_seconds_clipping)

            # BalancedDataset
            # -------------------------------
            start_balancing = time.perf_counter()

            if do_balance == True:
                print("Balancing sequences by positive/negative ...")
                fname2 = fname1 + '_BALANCED_' + str(pprop) + "_" + str(ntimes)
                from_file = os.path.join(data_result_dir, fname1)
                to_file = os.path.join(data_result_dir, fname2)
                try:
                    encoding.createBalancedDatasetFile(from_file + ".csv",
                                                       to_file + ".csv",
                                                       pprop, ntimes)
                except Exception as e:
                    logging.error(traceback.format_exc())
                    sys.exit(1)
            else:
                fname2 = fname1
                print("Balancing sequences not done!")

            end_balancing = time.perf_counter()
            elapsed_seconds_balancing = end_balancing - start_balancing
            conversion_balancing = datetime.timedelta(seconds=elapsed_seconds_balancing)

            # Coding
            # -------------------------------
            fname3 = fname2 + "_CODED"
            if n_ch == 1:
                fname3 += "1CH"
            elif n_ch == 2:
                fname3 += "2CH"
            elif n_ch == 3:
                fname3 += "3CH"
            if parallel_coding == False:
                myProgressReport = rna.ProgressReporter(2, "ProgressReporter2", encoding.progressRequest, delay=1)
                myProgressReport.daemon = True  # to stop thread at main termination
                myProgressReport.start()
            start_encoding = time.perf_counter()
            print("Encoding file {}.csv".format(fname2))
            encoding.encodeToFile(from_dir=data_result_dir,
                                  from_file=fname2 + ".csv",
                                  to_dir=data_result_dir,
                                  to_file=fname3 + ".csv",
                                  parallel=parallel_coding,
                                  encode_to="integer")
            end_encoding = time.perf_counter()
            elapsed_seconds_encoding = end_encoding - start_encoding
            conversion_encoding = datetime.timedelta(seconds=elapsed_seconds_encoding)

    # =========================================================
    end = time.perf_counter()
    elapsed_seconds = end - start
    conversion = datetime.timedelta(seconds=elapsed_seconds)
    print("TIEMPO EMPLEADO CLIPPING: {} SEGUNDOS / {}".format(elapsed_seconds_clipping, str(conversion_clipping)))
    print("TIEMPO EMPLEADO BALANCING: {} SEGUNDOS / {}".format(elapsed_seconds_balancing, str(conversion_balancing)))
    print("TIEMPO EMPLEADO ENCODING: {} SEGUNDOS / {}".format(elapsed_seconds_encoding, str(conversion_encoding)))
    print("TIEMPO TOTAL EMPLEADO: {} SEGUNDOS / {}".format(elapsed_seconds, str(conversion)))
    rna.stop_threads = True
    if parallel_coding == False:
        myProgressReport.join()


if __name__ == "__main__":
    main()
    sys.exit()