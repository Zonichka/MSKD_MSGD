import sympy as sym


N = 4  # номер в списке
A = sym.Matrix([[4, 1, 1], [1, (3+0.1*N)*2, -1], [1, -1, (4+0.1*N)*2]])  # задание матрицы квадратичной формы для f
b = sym.Matrix([1, -2, 3])
e = 10**(-6)  # погрешность
en = sym.Matrix.eye(3)  # матрица ортов


def f(x):  # задание функции
    return 2 * x[0]**2 + (3 + 0.1*N)*x[1]**2 + (4+0.1*N)*x[2]**2 + x[0]*x[1] - x[1]*x[2] + x[0]*x[2]\
           + x[0] - 2*x[1] + 3*x[2] + N


def x_kk(x):  # рассчет следующего элемента последовательности x_k для МНГС
    q_k = A * x + b
    m_k = - (q_k[0] * q_k[0] + q_k[1] * q_k[1] + q_k[2] * q_k[2]) / \
          (q_k[0] * (A * q_k)[0] + q_k[1] * (A * q_k)[1] + q_k[2] * (A * q_k)[2])

    return x + m_k * q_k


def x_k(x, en):  # рассчет следующего элемента последовательности x_k для МНКС
    m_k = - (en[0] * (A * x + b)[0] + en[1] * (A * x + b)[1] + en[2] * (A * x + b)[2]) / \
          (en[0] * (A * en)[0] + en[1] * (A * en)[1] + en[2] * (A * en)[2])

    return x + m_k * en


def MSKD():  # метод наискорейшего покоординатного спуска
    x = sym.Matrix([0, 0, 0])
    f0 = f(x) - 1
    fk = f(x)
    n = 0
    while abs(fk - f0) >= e:
        f0 = fk
        for i in range(3):
            xk = x_k(x, en[:, i % 3])
            if f(xk) < f0:
                fk = f(xk)
                x = xk
        print(fk, x)
    return fk


def MSGD():  # метод наискорейшего градиентного спуска
    e1 = sym.Matrix([1, 0, 0])
    f0 = f(e1) - 1
    fk = f(e1)
    x = e1
    n = 1
    while abs(fk - f0) >= e:
        f0 = f(x_kk(x))
        fk = f(x_kk(x_kk(x)))
        x = x_kk(x_kk(x))
        print(n, fk, x)
        n += 1
    return fk


print("MSGD:")
print(MSGD())
print("MSKD:")
print(MSKD())
