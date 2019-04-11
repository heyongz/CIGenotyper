from random import randint


class Table(object):
    """Table file class for complex indel simulation."""

    __reference = ''
    __idx = 0

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
                num = int(num / 50) * 2
                for i in range(num):
                    ss += self.__genVariants(0)       # Generate homozygote.
                    ss += self.__genVariants(1)       # Generate heterozygote.
            file.write(ss)

    def __genVariants(self, type):
        """Generate variants information.
        :param type: An integer (0 or 1) indicating type of variants.
        :return: A string contains variants information.
        """
        pos = 0
        ss = ''

        # Generate 50 variants every call to ensure the mutation rate below 5%.
        # Every 2 lines providing an complex indel information.
        for i in range(1, 50 + 1, 2):
            pos += randint(15000, 20000)

            if type == 0:
                label = 3       # Homozygote.
            else:
                label = 1       # Heterozygote.

            num = randint(500, 1500)
            ss += (str(self.__idx + i) + '\t' + str(pos) + '\t-' + str(num) + '\t' +
                   str(label) + '\t1\t0\t' + self.__reference[pos - 1:pos - 1 + num] + '\n')
            num = randint(500, 1500)
            ss += (str(self.__idx + i + 1) + '\t' + str(pos) + '\t' + str(num) + '\t' +
                   str(label) + '\t1\t0\t' + self.__gensequence(num) + '\n')
        self.__idx += 50
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
        brkpntL: An integer indicates the left margin of complex indel.
        brkpntR: An integer indicates the right margin of complex indel.
    """
    label = 0
    brkpntl = 0
    brkpntr = 0

    def __init__(self,label, brkpntL, brkpntR):
        """ Initialize an instance.
        :param label: An integer indicates complex indel's label.
        :param brkpntL: An integer indicates the left margin of complex indel.
        :param brkpntR: An integer indicates the right margin of complex indel.
        """

        # CIndel attributes
        self.label = label
        self.brkpntl = brkpntL
        self.brkpntr = brkpntR


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