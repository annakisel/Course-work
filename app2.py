import math
import numpy as np
import matplotlib.pyplot as pl


class RandomProcesses2:

    def __init__(self):
        self.n = 100
        self.q = 5
        self.sv = []
        self.matrixRP = [[0.] * self.n for _ in range(self.q)]
        self.w = np.array([1.0, 1/2.0, 1/4.0, 1/8.0, 1/16.0], float)
        # self.w = np.array([1.0, 2.0, 3.0, 4.0, 5.0], float)
        self.z = [0.] * self.n
        self.b = np.array([0.1, 0.15, 0.2, 0.3, 0.25], float)
        # self.b = np.array([1.0, 2.0, 3.0, 4.0, 5.0], float)
        self.s_z = [0.] * int(self.n * 2 / 3)

    def R(self, h, w):
        return (1 - w * abs(h)) if abs(h) <= 1. / w else 0

    def write_rv_in_file(self):
        f = open('random-variables_for_second_process.txt', 'w+')
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
        f = open('random-variables_for_second_process.txt', 'r')
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
            N = int(1 / self.w[i])
            for j in range(1, self.n):
                for k in range(0, N):
                    self.matrixRP[i][j] += 1. / math.sqrt(N) * self.sv[i][j - k]

            pl.ylabel('$Y$' + str(i) + '$(t)$')
            pl.xlabel(r'$t$')
            pl.plot(x, self.matrixRP[i], marker='.')
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
        # x = np.linspace(0, 100, self.n).reshape(-1, 1)
        R = [0.] * self.n
        for i in range(0, self.q):
            for h in range(0, self.n):
                R[h] = self.R(h, self.w[i])
            # pl.ylabel('$R$' + str(i) + '$(t)$')
            # pl.xlabel(r'$t$')
            # pl.plot(x, R)
            # pl.title('Ковариационная функция, w = ' + str(self.w[i]))
            # pl.show()
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
        pl.plot(x, self.z, marker='.')
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

    def D_before_beta(self, w_j, w_p, t1, t2, h):
        R = self.R
        return (4 * R(t1 - t2, w_j) * R(t1 - t2, w_p) - 2 * R(t1 - t2, w_j) * R(t1 - t2 - h, w_p) + R(t1 - t2 - h, w_j)
                * R(t1 - t2 - h, w_p) - 2 * R(t1 - t2, w_j) * R(t1 - t2 + h, w_p) + 2 * R(t1 - t2 - h, w_j) * R(t1 -
                                                                                                                t2 + h,
                                                                                                                w_p) - 2 * R(
                    t1 - t2 - h, w_j) * R(t1 - t2, w_p) + R(t1 - t2 + h, w_j) * R(t1 - t2 + h, w_p) - 2 *
                R(t1 - t2 + h, w_j) * R(t1 - t2, w_p))

    def D_gamma_z(self):
        h = 30
        summ = 0
        numbers = [50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500]
        d_gamma_estimate = [0.] * len(numbers)
        for count in range(0, len(numbers)):
            for t1 in range(1, numbers[count] - h):
                for t2 in range(1, numbers[count] - h):
                    for j in range(0, self.q):
                        for p in range(0, self.q):
                            summ += ((self.b[j] ** 2) * (self.b[p] ** 2) * self.D_before_beta(self.w[j],
                                                                                              self.w[p], t1, t2,
                                                                                              h))
            d_gamma_estimate[count] = 1 / (2 * (numbers[count] - h) ** 2) * summ
            summ = 0
        print(d_gamma_estimate)
        pl.ylabel(r'$Dy\^ (h)$')
        pl.xlabel(r'$t$')
        pl.plot(numbers, d_gamma_estimate)
        pl.title('Зависимость Dy^(h) и t, h =' + str(h))
        pl.show()


rp = RandomProcesses2()
rp.write_rv_in_file()
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
# rp.D_gamma_z()
