import threading
import sys
import time

# =============================================================================================
# L O G   R E P O R T I N G
# =============================================================================================

stop_threads = False


class ProgressReporter(threading.Thread):
    def __init__(self, threadID, name, progressRequest, delay=5):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.progressRequest = progressRequest
        self.delay = delay

    def run(self):
        self.report_progress()

    def report_progress(self):
        current_value, target_value = self.progressRequest()
        while (current_value == 0 and target_value == 0):
            sys.stdout.write("\r Progress Reporter Waiting...")
            current_value, target_value = self.progressRequest()
            time.sleep(2)
            # global stop_threads
            if stop_threads:
                break
        while (current_value < target_value):
            current_value, target_value = self.progressRequest()
            sys.stdout.write("\r Processing unit: {} of {} \r".format(current_value, target_value))
            time.sleep(self.delay)
            # global stop_threads
            if stop_threads:
                break
        print("\r ProgressReporter Ending \r")

# =============================================================================================