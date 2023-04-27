import numpy as np
import math as ma

class RBF:
    # J, xDim
    # bias, w[J], m[J][xDim], sigma[J]: 16 parameters
    def __init__(self, J=3, xDim=3):
        self.J = J
        self.xDim = xDim
        self.bias = 0
        self.w = np.zeros((J), np.float32)
        self.m = np.zeros((J, xDim), np.float32)
        self.sigma = np.zeros((J), np.float32)

    def calculateOutput(self, x):
        f = 0.0
        for i in range(0, self.J):
            f += self.w[i] * self.getPhi(x, i)
        f += self.bias
        return f

    def getPhi(self, x, index):
        return ma.exp(- self.getSquare(x, self.m[index]) / (2 * (self.sigma[index] ** 2)))

    def getSquare(self, x, mean):
        ret = 0.0
        for i in range(0, self.xDim):
            ret += (x[i] - mean[i]) ** 2
        return ret
