import RBF as rbf
import numpy as np
import random
import re


class Gene:
    J = 3
    xDim = 3
    DNALength = 1 + J + J * xDim + J
    f = 0.0
    bias_min = 0
    bias_max = 1
    w_min = 0
    w_max = 1
    m_min = 0
    m_max = 30
    sigma_min = 1e-6
    sigma_max = 10

    def __init__(self):
        self.myrbf = rbf.RBF(self.J, self.xDim)
        self.DNA = np.zeros((self.DNALength), np.float64)

    def generate(self):
        for i in range(0, 1):
            self.DNA[i] = random.uniform(self.bias_min, self.bias_max)

        for i in range(1, 1 + self.J):
            self.DNA[i] = random.uniform(self.w_min, self.w_max)

        for i in range(1 + self.J, 1 + self.J + self.J * self.xDim):
            self.DNA[i] = random.uniform(self.m_min, self.m_max)

        for i in range(1 + self.J + self.J * self.xDim, 1 + self.J + self.J * self.xDim + self.J):
            self.DNA[i] = random.uniform(self.sigma_min, self.sigma_max)

    def clone(self):
        g = Gene()
        for i in range(0, self.DNALength):
            g.DNA[i] = self.DNA[i]
        g.f = self.f
        g.setrbf()
        return g

    # input x[N][xDim], output[N] 計算適應函數直
    def calculateFitness(self, input, output):  
        self.setrbf()
        f = 0.0
        for i in range(0, len(output)):
            f += (output[i] - self.myrbf.calculateOutput(input[i])) ** 2
        f /= 2
        self.f = f

    def setrbf(self):
        for i in range(0, 1):
            self.myrbf.bias = min(max(self.DNA[i], self.bias_min), self.bias_max)
            self.DNA[i] = self.myrbf.bias

        j = 0
        for i in range(1, 1 + self.J):
            self.myrbf.w[j] = min(max(self.DNA[i], self.w_min), self.w_max)
            self.DNA[i] = self.myrbf.w[j]
            j += 1

        j = 0
        for i in range(1 + self.J, 1 + self.J + self.J * self.xDim):
            self.myrbf.m[int(j / self.xDim)][j % self.xDim] = min(max(self.DNA[i], self.m_min), self.m_max)
            self.DNA[i] = self.myrbf.m[int(j / self.xDim)][j % self.xDim]
            j += 1

        j = 0
        for i in range(1 + self.J + self.J * self.xDim, 1 + self.J + self.J * self.xDim + self.J):
            self.myrbf.sigma[j] = min(max(self.DNA[i], self.sigma_min), self.sigma_max)
            self.DNA[i] = self.myrbf.sigma[j]
            j += 1

    def setGene(self, geneStr):
        geneList = re.findall(r'[\w.+-]+', geneStr)
        for i in range(0, len(geneList)):
            self.DNA[i] = geneList[i]
        self.setrbf()

    def getTheta(self, input):
        ret = self.myrbf.calculateOutput(input)
        ret *= 80
        ret -= 40
        return ret