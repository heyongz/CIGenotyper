from BasicClasses import Table, VCF
from _pickle import dump, load
from FeatureCollection import *
from SparseBayes import RVC
import numpy as np
import sys
import os


def collectfeature(num, isvariant):
    """Generate complex indels, call BWA software, and features of collect complex indels.
    :param num: An integer indicates the number of complex indels/normal sequence we want to generate.
    :param isvariant: A boolean indicates whether to generate complex indels or not.
    :return: A list contains extracted label and feature information.
    """

    Table(fapath, tablepath, isvariant, num)
    os.system('./'+bashpath)
    vcindel = VCF(tablepath, isvariant, num)
    fc = FeatureCollection(vcindel.getdata())
    fc.setparams(sampath, insertsize=500, stdvar=50, readlength=100)
    return fc.collectfeature()


def gendata(num):
    """Generate dataset which includes normal sequence and complex indels.
    :param num: An integer indicates the number of complex indels and normal sequence to generate.
    :return: A list contains mixed label and feature of both complex indels and normal sequence.
    """
    [label1, feature1] = collectfeature(num, True)
    [label2, feature2] = collectfeature(num, False)
    return [np.array(label1+label2).astype(np.float), np.array(feature1+feature2).astype(float)]


def writedata(num):
    """Write extracted feature to file.
    :param num: An integer indicates the scale of dataset.
    :return: None
    """
    dump(gendata(num), open("FILES/DATA/data_"+str(num)+".txt", 'wb'))


def loaddata(num):
    """Load features from file.
    :param num: An integer indicates the scale of dataset.
    :return: A N*M list contains extracted features, N samples and M features.
    """
    return load(open("FILES/DATA/data_"+str(num)+".txt", 'rb'))


if __name__ == "__main__":

    fapath = r'FILES/chr19_100w.fa'
    sampath = r'FILES/result.sam'
    tablepath = r'FILES/table.txt'
    modelpath = r'FILES/MODEL/model_200.rvc'
    bashpath = r'bash.sh'

    assert sys.platform, "linux"
    assert os.path.exists(fapath), "Fasta file not found."
    assert os.path.exists(sampath), "Sam file not found."
    assert os.path.exists(tablepath), "Table file not found."
    assert os.path.exists(bashpath), "Bash file not found."

    # [label, feature] = loaddata(200)
    # clf = RVC().fit(feature, label)
    # dump(clf, open(modelpath, 'wb'))

    clf = load(open(modelpath, 'rb'))
    [label, feature] = loaddata(100)
    print("Accuracy: ", clf.score(feature, label)*100, "%", end='')