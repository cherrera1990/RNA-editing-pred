import concurrent.futures
import numpy as np


# =============================================================================================
# L I S T   P A R A L L E L   P R O C E S S O R
# =============================================================================================

class ListParallelProcessor():

    def __init__(self, plist):
        self.input_list = plist
        self.num_processors = mp.cpu_count()
        self.output_list = list()

    # -----------------------------------------------
    def parallelProcessList(self, processor, *args, **kwargs):
        # split input list into num_processors list
        splitted_list = np.array_split(self.input_list, self.num_processors)
        # Process Pool
        with concurrent.futures.ProcessPoolExecutor(max_workers=self.num_processors) as executor:
            futures = [executor.submit(processor, elist, *args, **kwargs) for elist in splitted_list]
            # this line waits for all the process to finish
            for future in concurrent.futures.as_completed(futures):
                try:
                    result_series = future.result()
                    result_list = result_series.tolist()
                    # print("result_list[{}]->{}\n".format(type(result_list),result_list))
                    # very important: to use extend not append
                    self.output_list.extend(result_list)
                except Exception as ex:
                    print("Exception during parallelProcessList:" + str(ex))
                    pass

        return (self.output_list)
    # -----------------------------------------------
