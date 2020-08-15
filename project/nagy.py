import numpy as np
import matplotlib.pyplot as plt


NAGY_KAPPA15 = "../data/15M_sun.csv"
NAGY_KAPPA30 = "../data/30M_sun.csv"
t_k15, kappa_t15 = np.loadtxt(NAGY_KAPPA15, unpack=True, delimiter=", ")

times = np.linspace(4, 365, 1000)


def smooth_min(a, b, k):
    a = pow(a, k); b = pow(b, k)
    return pow((a*b)/(a+b), 1/k)


# def smooth_poly(a, b, k):
#     h = max(k-abs(a-b), 0.0)/k
#     return min(a, b) - h*h*h*k*(1.0/6.0)


def part1(t):
    h = 4.47
    j = 0.613
    k = 0.1
    g = 0.264
    return g + h * (t**k) * (j**t)


def part2(t):
    a = -0.013
    b = 1.1
    c = 3.4
    d = 0.8
    n = 1.39
    m = -0.15
    z = -1
    return a + b * np.sqrt(t - c) / (d * t) + (n * np.e ** (m * t + z))


def part3(t):
    return 0.1207 +0*t  # 0*t to work with ndarrays


def tsp_nagy(t):
    k1 = 30  # part1-part2 smoothing value
    k2 = -50  # part2-part3 smoothing value | negative means use maximum
    if isinstance(t, (list, np.ndarray)):
        y = []
        for i in t:
            # up until intersect between part2 and part3:
            y.append(smooth_min(smooth_min(part1(i), part2(i), k1), part3(i), k2))
        return y
    elif isinstance(t, (int, float)):
        return smooth_min(smooth_min(part1(t), part2(t), k1), part3(t), k2)
    # elif isinstance(t, (float, int)):
    #     return smooth_min(smooth_min(smooth_min(part1(t), part2(t), 12), part3(t), k2))
    else:
        raise(TypeError("argument t must be a number or a list/numpy array of numbers"))


plt.semilogx(times, part1(times), color="green", linestyle='--', label="function 1")
plt.semilogx(times, part2(times), color="orange", linestyle='--', label="function 2")
plt.semilogx(times, part3(times), color="purple", linestyle='--', label="function 3")
plt.semilogx(times, tsp_nagy(times), color="red", label="Combined function")
plt.scatter(t_k15, kappa_t15, color="blue", label="Nagy's data")
plt.legend()
plt.xlim(3.9,400)
plt.ylim(-0.1,0.9)
plt.xlabel("time (days)")
plt.ylabel("kappa")
plt.show()
