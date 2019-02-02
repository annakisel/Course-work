import math
import numpy as np
import matplotlib.pyplot as pl


class RandomProcesses:

    def __init__(self):
        self.n = 100
        self.q = 5
        self.sv = np.empty(shape=self.n)
        self.matrixRP = [[0.] * self.n for _ in range(self.q)]
        self.r = np.array([1.5, 1.0, 1.0 / 2.0, 1.0 / 4.0, 1.0 / 8.0], float)
        self.z = [0.] * self.n
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
        pl.ylabel(r'$eta(t)$')
        pl.xlabel(r'$t$')
        x = np.linspace(0, 100, self.n).reshape(-1, 1)
        for i in range(0, self.q):
            ro = math.exp(-1.0 / self.r[i])
            b1 = -ro
            a0 = math.sqrt(1 - ro * ro)
            self.matrixRP[i][0] = a0 * self.sv[0]
            for j in range(0, self.n):
                self.matrixRP[i][j] = a0 * self.sv[j] - b1 * self.matrixRP[i][j - 1]
            pl.ylabel(r'$eta(t)$')
            pl.xlabel(r'$t$')
            pl.plot(x, self.matrixRP[i])
            pl.axis([0, 110, -3, 3])
            pl.title('Случайный процесс, r = ' + str(self.r[i]))
        pl.show()

    def semivarams_and_estimates(self):
        x = np.linspace(0, 100, self.n).reshape(-1, 1)
        estimate = [0.] * self.n
        semivar = [0.] * self.n
        for j in range(0, self.q):
            for h in range(0, self.n):
                semivar[h] = self.R(0, self.r[j]) - self.R(h, self.r[j])
                summ = 0
                for i in range(0, self.n - h):
                    summ += (self.matrixRP[j][i] - self.matrixRP[j][i + h]) * (
                            self.matrixRP[j][i] - self.matrixRP[j][i + h])
                estimate[h] = summ / 2.0 / (self.n - h)
            pl.title('Семивариагарамма и ее оценка, r = ' + str(self.r[j]))
            pl.ylabel(r'$y(t)$')
            pl.xlabel(r'$t$')
            pl.plot(x, estimate)
            pl.plot(x, semivar)
        pl.show()

    def plots_of_cov_and_d(self):
        x = np.linspace(0, 100, self.n).reshape(-1, 1)
        R = [0.] * self.n
        for i in range(0, self.q):
            for h in range(0, self.n):
                R[h] = self.R(h, self.r[i])
            pl.ylabel(r'$R(t)$')
            pl.xlabel(r'$t$')
            pl.plot(x, R)
            pl.title('Ковариационная функция, r = ' + str(self.r[i]))
            pl.show()
            print(str(i) + ' process\nсреднее значение: ' + str(self.X_(self.matrixRP[i])) + '\nдисперсия: ' + str(
                self.Dx(self.matrixRP[i])) + '\n')

    def X_(self, x):
        s = 0
        for i in range(0, self.n):
            s += x[i]
        return s / self.n

    def Dx(self, x):
        d = 0
        x_ = self.X_(x)
        for i in range(0, self.n):
            d += (x[i] - x_) * (x[i] - x_)
        return d / self.n

    def z_modelling(self):
        x = np.linspace(0, 100, self.n).reshape(-1, 1)
        for h in range(0, self.n):
            for j in range(0, self.q):
                self.z[h] += self.b[j] * self.matrixRP[j][h]
        pl.ylabel(r'$eta(t)$')
        pl.xlabel(r'$t$')
        pl.plot(x, self.z)
        pl.title('Случайный процесс z ')
        pl.show()
        print('-------z process -----------')
        print('среднее значение: ' + str(self.X_(self.z)) + '\nдисперсия: ' + str(
            self.Dx(self.z)) + '\n')


rp = RandomProcesses()
# rp.write_rv_in_file()
rp.read_rv_from_file()
rp.creating__random_process()
rp.semivarams_and_estimates()
rp.plots_of_cov_and_d()
rp.z_modelling()
