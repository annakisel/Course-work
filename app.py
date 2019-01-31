import math
import numpy as np
import matplotlib.pyplot as pl


class RandomProcesses:

    def __init__(self):
        self.n = 100
        self.sv = np.empty(shape=self.n)
        self.process = []
        self.r = np.array([1.0, 2.0, 3.0, 4.0], float)
        self.b = np.array([1.0, 2.0, 3.0, 4.0], float)

    def R(self, h, r):
        return math.exp(-1 / r * h)

    def write_rv_in_file(self):
        self.sv = np.random.normal(size=self.n)

        f = open("random-variables.txt", "w+")
        for i in range(0, self.n):
            f.write(str(self.sv[i]) + " ")
        f.close()

    def read_rv_from_file(self):
        self.sv = np.fromfile("random-variables.txt", float, self.n, ' ')

    def creating__random_process(self):
        ro = math.exp(-1.0 / self.r[0])
        b1 = -ro
        a0 = math.sqrt(1 - ro * ro)
        result = list(range(0, self.n))
        result[0] = a0 * self.sv[0]
        for i in range(1, self.n):
            result[i] = a0 * self.sv[i] - b1 * result[i - 1]
        print('result\n' + str(result))

        x = list(range(0, 100))
        pl.ylabel(r'$eta(t)$')
        pl.xlabel(r'$t$')
        pl.plot(x, result)
        pl.axis([0, 110, -3, 3])
        pl.title('Случайный процесс')
        pl.show()


rp = RandomProcesses()
# rp.write_rv_in_file()
rp.read_rv_from_file()
rp.creating__random_process()
