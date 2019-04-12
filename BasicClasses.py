from random import randint
import logging


class Table(object):
    """Table file class for complex indel simulation."""

    __reference = ''

    def __init__(self, fapath, tablepath, isvariant, num):
        """Initialize an instance.
        :param fapath: A string represents the path of fasta file.
        :param tablepath: A string represents the path of table file.
        :param isvariant: A boolean suggests whether to generate variants.
        :param num: An integer count of variants. Needless if you don't generate variants.
        """

        with open(fapath, 'r') as file:
            file.readline()
            s = file.read()
            s = s.split()
            for ss in s:
                self.__reference += ss

        # Table file will be empty if it's not used for variants generation.
        with open(tablepath, 'w') as file:
            ss = ""
            if isvariant:
                ss = self.__genVariants(num)
            file.write(ss)

    def __genVariants(self, num):
        """Generate variants information.
        :param num: An integer indicating number of variants.
        :return: A string contains variants information.
        """
        pos = 0
        ss = ''

        """
        Total length of reference is 1000,000.
        Step value between complex indel is 100000/2/n.
        TODO: Automatically set step.
        """
        def __genpos():
            """Generate position of complex indel.
            :return: An integer indicates the position of complex indel.
            """
            return randint(600, 900)

        def __genlen():
            """Generate length of deletion or insertion.
            :return: An integer indicates the length.
            """
            return randint(200, 800)

        # Generate num complex indels for each kind of cindel.
        # Totally generate 2*num complex indels.
        # Note: two lines in table file provide one complex indel information.
        try:
            for i in range(1, 2*num + 1, 4):
                # Homozygote.
                pos += __genpos()
                length = __genlen()
                ss += (str(i) + '\t' + str(pos) + '\t-' + str(length) + '\t3' +
                       '\t1\t0\t' + self.__reference[pos - 1:pos - 1 + length] + '\n')
                length = __genlen()
                ss += (str(i+1) + '\t' + str(pos) + '\t' + str(length) + '\t3' +
                       '\t1\t0\t' + self.__gensequence(length) + '\n')

                # Heterozygote.
                pos += __genpos()
                length = __genlen()
                ss += (str(i+2) + '\t' + str(pos) + '\t-' + str(length) + '\t1' +
                       '\t1\t0\t' + self.__reference[pos - 1:pos - 1 + length] + '\n')
                length = __genlen()
                ss += (str(i+3) + '\t' + str(pos) + '\t' + str(length) + '\t1' +
                       '\t1\t0\t' + self.__gensequence(length) + '\n')
        except Exception as e:
            logging.exception(e)

        return ss

    def __gensequence(self, num):
        """Random generation of gene sequence.
        :param num: An integer indicating the length of gene sequence.
        :return: A string contains generated gene sequence.
        """
        ss = ['A', 'G', 'C', 'T']
        sq = ''
        for i in range(num):
            sq += ss[randint(0, 3)]
        return sq


class VCF(object):
    """Complex indel information for feature collection."""

    __vcindel = []

    def __init__(self, tablepath, isvariant, num):
        """Initialize an instance.
        :param tablepath: A string represents the path of table file.
        :param isvariant: A boolean suggests whether to generate simulation data.
        :param num: A integer count of simulation data we want to generate.
        """

        if isvariant:      # Load complex indels information from table file.
            with open(tablepath, 'r') as file:
                lines = file.readlines()[::2]
                lines = [line.split('\t') for line in lines]
                for line in lines:
                    self.__vcindel.append(self.__splitline(line))
        else:       # Randomly generate normal data, add them to complex indels.
            for i in range(num):
                brkpntl = randint(10000, 1000000 - 10000)
                brkpntr = brkpntl + randint(500, 1500)
                self.__vcindel.append([-1, brkpntl, brkpntr])

    def getdata(self):
        """Return the list of all complex indels' feature.
        :return: A list contains the feature of all complex indels.
        """
        return self.__vcindel

    def __splitline(self, line):
        """Split line from table files and transform them to the features of complex indel.
        :param line: A string contains the information of complex indel.
        :return: A list contains the transformed feature of complex indel.
        """
        brkpntl = int(line[1])
        brkpntr = int(brkpntl + abs(int(line[2])))
        if int(line[3]) == 3:
            label = 0       # Homozygote.
        else:
            label = 1       # Heterozygote.
        return [label, brkpntl, brkpntr]


class CIndel(object):
    """Complex indel class.
    Attributes:
        label: An integer indicates complex indel's label.
        brkpntl: An integer indicates the left margin of complex indel.
        brkpntr: An integer indicates the right margin of complex indel.
    """
    label = 0
    brkpntl = 0
    brkpntr = 0

    def __init__(self,label, brkpntl, brkpntr):
        """ Initialize an instance.
        :param label: An integer indicates complex indel's label.
        :param brkpntl: An integer indicates the left margin of complex indel.
        :param brkpntr: An integer indicates the right margin of complex indel.
        """

        # CIndel attributes
        self.label = label
        self.brkpntl = brkpntl
        self.brkpntr = brkpntr


class Read(object):
    """Read class.
    Attributes:
        pos: An integer indicates the position of first read of paired-end reads.
        mapq: An integer indicates the mapped quality of read to the reference sequence.
        pnext: An integer indicates the position of second read of paired-end reads.
        mate: A instance of Read, indicates the mate of current read
    """
    pos = 0
    mapq = 0
    pnext = 0
    mate = None

    def __init__(self, pos, mapq, pnext):
        """Initialize an instance.
        :param pos: An integer indicates the position of first read of paired-end reads.
        :param mapq: An integer indicates the mapped quality of read to the reference sequence.
        :param pnext: An integer indicates the position of second read of paired-end reads.
        """
        self.pos = pos
        self.mapq = mapq
        self.pnext = pnext

    def setmate(self, mate):
        """ Set paired-end read's mate.
        :param mate: A instance of Read, indicates the mate of current read
        :return: None
        """
        self.mate = mate