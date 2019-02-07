import math
import numpy as np
import matplotlib.pyplot as pl


class RandomProcesses:

    def __init__(self):
        self.n = 100
        self.q = 5
        self.sv = []
        self.matrixRP = [[0.] * self.n for _ in range(self.q)]
        self.w = np.array([1.0, 1/2.0, 1/3.0, 1/4.0, 1/5.0], float)
        self.z = [0.] * self.n
        self.b = np.array([0.1, 0.15, 0.2, 0.3, 0.25], float)
        # self.b = np.array([1.0, 2.0, 3.0, 4.0, 5.0], float)
        self.s_z = [0.] * int(self.n * 2 / 3)

    def R(self, h, w):
        return math.exp(-w * h)

    def write_rv_in_file(self):
        f = open('random-variables.txt', 'w+')
        for j in range(0, self.q):
            rv = np.random.normal(size=self.n)
            for i in range(0, self.n):
                if i != self.n - 1:
                    f.write(str(rv[i]) + ',')
                else:
                    f.write(str(rv[i]))
            f.write('\n')
        f.close()

    def read_rv_from_file(self):
        f = open('random-variables.txt', 'r')
        if f.mode == 'r':
            lines = f.readlines()
            for line in lines:
                self.sv.append([float(x) for x in line.split(",")])
        else:
            print('file not opened')
        print('----SV----')
        for i in range(0, self.q):
            print(self.sv[i])

    def creating__random_process(self):
        pl.ylabel(r'$eta(t)$')
        pl.xlabel(r'$t$')
        x = np.linspace(0, 100, self.n).reshape(-1, 1)
        for i in range(0, self.q):
            ro = math.exp(-self.w[i])
            b1 = -ro
            a0 = math.sqrt(1 - ro * ro)
            self.matrixRP[i][0] = a0 * self.sv[i][0]
            for j in range(1, self.n):
                self.matrixRP[i][j] = a0 * self.sv[i][j] - b1 * self.matrixRP[i][j - 1]
            pl.ylabel('$Y$' + str(i) + '$(t)$')
            pl.xlabel(r'$t$')
            pl.plot(x, self.matrixRP[i])
            pl.axis([0, 110, -3, 3])
            pl.title('Случайный процесс, w = ' + str(self.w[i]))
            pl.show()

    def semivarams_and_estimates(self):
        x = np.linspace(0, int(self.n * 2 / 3), int(self.n * 2 / 3)).reshape(-1, 1)
        estimate = [0.] * int(self.n * 2 / 3)
        semivar = [0.] * int(self.n * 2 / 3)
        for j in range(0, self.q):
            for h in range(0, int(self.n * 2 / 3)):
                semivar[h] = self.R(0, self.w[j]) - self.R(h, self.w[j])
                summ = 0
                for i in range(0, self.n - h):
                    summ += (self.matrixRP[j][i] - self.matrixRP[j][i + h]) * (
                            self.matrixRP[j][i] - self.matrixRP[j][i + h])
                estimate[h] = summ / 2.0 / (self.n - h)
            pl.title('Семивариагарамма и ее оценка, w = ' + str(self.w[j]))
            pl.ylabel('$y$' + str(j) + '$(t)$')
            pl.xlabel(r'$t$')
            pl.plot(x, estimate, 'C1')
            pl.plot(x, semivar, 'C2')
            pl.show()

    def plots_of_cov_and_d(self):
        x = np.linspace(0, 100, self.n).reshape(-1, 1)
        R = [0.] * self.n
        for i in range(0, self.q):
            for h in range(0, self.n):
                R[h] = self.R(h, self.w[i])
            pl.ylabel('$R$' + str(i) + '$(t)$')
            pl.xlabel(r'$t$')
            pl.plot(x, R)
            pl.title('Ковариационная функция, w = ' + str(self.w[i]))
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
                self.z[h] += (self.b[j] * self.matrixRP[j][h])
        pl.ylabel('Z(t)')
        pl.xlabel(r'$t$')
        pl.plot(x, self.z)
        pl.title('Случайный процесс Z(t) ')
        pl.axis([0, 110, None, None])
        pl.show()
        print('-------z process -----------')
        print('среднее значение: ' + str(self.X_(self.z)) + '\nдисперсия: ' + str(
            self.Dx(self.z)) + '\n')

    def S_z(self):
        for h in range(0, int(self.n * 2 / 3)):
            for j in range(0, self.q):
                self.s_z[h] += (self.b[j] * self.b[j] * (1.0 - self.R(h, self.w[j])))

    def R_z(self):
        r_z = [0.] * self.n
        x = np.linspace(0, 100, self.n).reshape(-1, 1)
        for h in range(0, self.n):
            for j in range(0, self.q):
                r_z[h] += (self.b[j] * self.b[j] * self.R(h, self.w[j]))
        pl.ylabel(r'$R(t)$')
        pl.xlabel(r'$t$')
        pl.plot(x, r_z)
        pl.title('Ковариационная функция процесса Z(t) ')
        pl.show()

    def estimate_sem_z(self):
        self.S_z()
        estimate_s_z = [0.] * int(self.n * 2 / 3)
        x = np.linspace(0, int(self.n * 2 / 3), int(self.n * 2 / 3)).reshape(-1, 1)
        for h in range(0, int(self.n * 2 / 3)):
            for j in range(0, self.n - h):
                estimate_s_z[h] += ((self.z[j] - self.z[j + h]) * (self.z[j] - self.z[j + h]))
            estimate_s_z[h] /= 2.0
            estimate_s_z[h] /= (self.n - h)
        pl.ylabel(r'$y(t)$')
        pl.xlabel(r'$t$')
        pl.plot(x, self.s_z, 'C1')
        pl.plot(x, estimate_s_z, 'C2')
        pl.title('Семивариограмма процесса Z(t) и ее оценка ')
        pl.show()


rp = RandomProcesses()
# rp.write_rv_in_file()
rp.read_rv_from_file()
rp.creating__random_process()
rp.semivarams_and_estimates()
rp.plots_of_cov_and_d()
rp.z_modelling()
rp.R_z()
rp.estimate_sem_z()
print('сумма квадратов')
summ = 0
for i in range(0, rp.q):
    summ += rp.b[i] * rp.b[i]
print(summ)
