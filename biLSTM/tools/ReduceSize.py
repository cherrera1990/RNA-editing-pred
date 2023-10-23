import rnaed as rna
import os
import time
import datetime
import sys




def main():
    #------ INIT ---------------------------
    data_dir =      '../data/datasets/ALL_HG38/W50'
    data_result_dir='../data/datasets/ALL_HG38/W50'
    pprop=0.0625
    ntimes=1
    categories=[['A','G','C','T'],['s','d','h','i','b']]
    encoding=rna.DataEncoder(categories, padding=True, pad_char="*")
    #BalancedDataset
    #-------------------------------
    fname1 = "hg38_RS_DL_Pad_W50_BALANCED_1_1_CODED_REDUCED.csv"
    fname2 = "hg38_RS_DL_Pad_W50_BALANCED_1_1_CODED_SMALL3.csv"
    from_file = os.path.join(data_dir,fname1)
    to_file   = os.path.join(data_result_dir,fname2)
    start_balancing=time.perf_counter()
    encoding.createBalancedDatasetFile(from_file,
                                       to_file,
                                       pprop,ntimes)
    end_balancing=time.perf_counter()
    elapsed_seconds_balancing=end_balancing-start_balancing
    conversion_balancing=datetime.timedelta(seconds=elapsed_seconds_balancing)
            
    #=========================================================
    print("TIEMPO EMPLEADO BALANCING: {} SEGUNDOS / {}".format(elapsed_seconds_balancing,str(conversion_balancing)))  


if __name__ == "__main__":
    main()
    sys.exit()      