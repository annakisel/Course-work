import math
import numpy as np
import matplotlib.pyplot as pl


class RandomProcesses:

    def __init__(self):
        self.n = 100
        self.q = 5
        self.sv = np.empty(shape=self.n)
        self.matrixRP = [[0] * self.n] * self.q
        self.r = np.array([1.0, 2.0, 3.0, 4.0, 5.0], float)
        self.b = np.array([1.0, 2.0, 3.0, 4.0, 5.0], float)

    def R(self, h, r):
        return math.exp(-1 / r * h)

    def write_rv_in_file(self):
        self.sv = np.random.normal(size=self.n)

        f = open('random-variables.txt', 'w+')
        for i in range(0, self.n):
            f.write(str(self.sv[i]) + ' ')
        f.close()

    def read_rv_from_file(self):
        self.sv = np.fromfile('random-variables.txt', float, self.n, ' ')

    def creating__random_process(self):
        for i in range(0, self.q):
            ro = math.exp(-1.0 / self.r[i])
            b1 = -ro
            a0 = math.sqrt(1 - ro * ro)
            self.matrixRP[i][0] = a0 * self.sv[0]
            for j in range(1, self.n):
                self.matrixRP[i][j] = a0 * self.sv[j] - b1 * self.matrixRP[i][j - 1]

            x = list(range(0, 100))
            pl.ylabel(r'$eta(t)$')
            pl.xlabel(r'$t$')
            pl.plot(x, self.matrixRP[i])
            pl.axis([0, 110, -3, 3])
            pl.title('Случайный процесс ' + str(i))
            pl.show()


rp = RandomProcesses()
# rp.write_rv_in_file()
rp.read_rv_from_file()
rp.creating__random_process()
