import skrvm
from os import system
from sklearn.preprocessing import StandardScaler


class RVC(object):
    """A simple wrapper for skrvm.RVC()"""

    __scaler = None
    __clf = None

    # TODO: Provide posterior probability.

    def fit(self, feature, label):
        """Generate RVM model.
        :param feature: A N*M list contains N samples and M features.
        :param label: A 1*N list contains labels of input samples.
        :return: Return self for cascading calling.
        """
        self.__scaler = StandardScaler()
        feature = self.__scaler.fit_transform(feature)
        self.__clf = skrvm.RVC().fit(feature, label)
        return self

    def predict(self, feature):
        """Predict labels using calculated model.
        :param feature: A N*M list contains N samples and M features.
        :return: A 1*N list contains predicted labels.
        """
        self.__check()
        feature = self.__scaler.transform(feature)
        return self.__clf.predict(feature)

    def score(self, feature, label):
        """Calculate the prediction accuracy.
        :param feature: A N*M list contains N samples and M features.
        :param label: A 1*N list contains labels of input samples.
        :return: A floating-number indicates the accuracy of prediction.
        """
        self.__check()
        feature = self.__scaler.transform(feature)
        return self.__clf.score(feature, label)

    def __check(self):
        """Make sure RVC initialized before use.
        :return: None
        """
        assert self.__scaler is not None, "Scaler not initialized!"
        assert self.__clf is not None, "Classifier not initialized!"


class SVC(object):
    """For test purpose only.
    To compare the result with relevant vector machine classifier.
    """
    pass

    # from main import gendata
    #
    # [label, feature] = gendata(200)
    # with open("FILES/train.txt", 'w') as file:
    #     ss = ""
    #     for i in range(len(label)):
    #         ss += "{0} 1:{1} 2:{2} 3:{3} 4:{4}\n".format(label[i], feature[i][0], feature[i][1],
    #                                                      feature[i][2], feature[i][3])
    #     file.write(ss)
    #
    # [label, feature] = gendata(200)
    # with open("FILES/test.txt", 'w') as file:
    #     ss = ""
    #     for i in range(len(label)):
    #         ss += "{0} 1:{1} 2:{2} 3:{3} 4:{4}\n".format(label[i], feature[i][0], feature[i][1],
    #                                                      feature[i][2], feature[i][3])
    #     file.write(ss)
    # system("python3 Tools/libsvm/tools/easy.py FILES/train.txt FILES/test.txt")
