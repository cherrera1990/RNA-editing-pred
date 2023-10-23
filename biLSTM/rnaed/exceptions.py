# =============================================================================================
# E X C E P T I O N S
# =============================================================================================
# define Python user-defined exceptions
class Error(Exception):
    """Base class for other exceptions"""
    pass


class LenSeqError(Error):
    def __init__(self, seq_header):
        self.seq_header = seq_header
        self.message = "Ignoring gene " + str(self.seq_header)
        self.message += " because channels are different in length"

    def __str__(self):
        return ("Error length in sequence")


class NumAdeninesDifferent(Error):
    def __init__(self, seq_header, n_aA, n_12):
        self.seq_header = seq_header
        self.n_aA = n_aA
        self.n_12 = n_12
        self.message = "Ignoring gene " + str(self.seq_header)
        self.message += " because  num adenosines, a or A=" + str(self.n_aA)
        self.message += ", not equal 1's or 2's=" + str(self.n_12)

    def __str__(self):
        return ("Error num adenines in sequence")


class NumParenthesisDifferent(Error):
    def __init__(self, seq_num):
        self.seq_num = seq_num
        self.message = "Ignoring gene " + str(self.seq_num)
        self.message += " because number of left parenthesis different from right parenthesis"

    def __str__(self):
        return ("Error parenthesis in sequence")


class MissingEditingFlag(Error):
    def __init__(self, seq_name, adenosine_num):
        self.adenosine_num = adenosine_num
        self.seq_name = seq_name
        self.message = "Error: missing editing flag for sequence " + str(self.seq_name)
        self.message += " adenosine number " + str(self.adenosine_num)

    def __str__(self):
        return ("Error editing Flag")


class InvalidCharInGeneSeq(Error):
    def __init__(self, seq_name, adenosine_num):
        self.adenosine_num = adenosine_num
        self.seq_name = seq_name
        self.message = "Ignoring gene fragment " + str(self.seq_name)
        self.message += " adenosine position " + str(self.adenosine_num)
        self.message += " because of invalid char"

    def __str__(self):
        return ("Error gene sequence char not valid")


class InvalidCharInPairsSeq(Error):
    def __init__(self, seq_name, adenosine_num):
        self.adenosine_num = adenosine_num
        self.seq_name = seq_name
        self.message = "Ignoring pairs fragment " + str(self.seq_name)
        self.message += " adenosine position " + str(self.adenosine_num)
        self.message += " because of invalid char"

    def __str__(self):
        return ("Error pairs sequence char not valid")


class InvalidCharInSecondarySeq(Error):
    def __init__(self, seq_name, adenosine_num):
        self.adenosine_num = adenosine_num
        self.seq_name = seq_name
        self.message = "Ignoring secondary fragment " + str(self.seq_name)
        self.message += " adenosine position " + str(self.adenosine_num)
        self.message += " because of invalid char"

    def __str__(self):
        return ("Error secondary sequence char not valid")


class InvalidCharInEditingSeq(Error):
    def __init__(self, seq_num):
        self.seq_num = seq_num
        self.message = "Ignoring invalid char for editing fragment " + str(self.seq_num)

    def __str__(self):
        return ("Error editing sequence char not valid")

class PositiveCasesNotFound(Error):
    def __init__(self, file_name):
        self.file_name = file_name
        self.message = "Positive cases of editing not found at file " + str(self.file_name)

    def __str__(self):
        return ("Error Positive cases of editing not found")


class NotEnoughGenes(Error):
    def __init__(self, ngenes):
        self.ngenes = ngenes
        self.message = "Error!, There are less than " + str(self.ngenes) + " genes at dataset"

    def __str__(self):
        return ("Error!, there are less genes than expected")
# =============================================================================================