import numpy as np

class QoI(object):
    def __init__(self):
        pass

    def __call__(self):
        return None

class QoIXsquared(QoI):
    def __init__(self):
        pass

    def __call__(self,X):
        return X[0]**2
