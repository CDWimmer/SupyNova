import numpy as np


def smooth_min(a, b, k):
    """Smoothly joins two values. Shapes to min(a, b) when k>0, max when k<0"""
    a = pow(a, k); b = pow(b, k)
    return pow((a*b)/(a+b), 1/k)


# def smooth_poly(a, b, k):
#     h = max(k-abs(a-b), 0.0)/k
#     return min(a, b) - h*h*h*k*(1.0/6.0)


def part1(t):
    a = 0.583
    b = 1.72
    c = 4.34
    return a - np.exp((- b * t) + c)


def part2(t):
    h = 4.47
    j = 0.613
    k = 0.1
    g = 0.264
    return g + h * (t**k) * (j**t)


def part3(t):
    a = -0.013
    b = 1.1
    c = 3.4
    d = 0.8
    n = 1.39
    m = -0.15
    z = -1
    # This is horrid and a source of slowdown, I need to find another function to describe this part that doesn't use
    # a square root.
    # a + b * np.sqrt(abs(t - c)) / (d * t) + (n * np.e ** (m * t + z))
    if isinstance(t, (list, np.ndarray)):
        res = []
        for i in t:
            if i <= c:
                res.append(part4(i))
            else:
                res.append(a + b * np.sqrt(abs(i - c)) / (d * i) + (n * np.e ** (m * i + z)))
        print(res)
        return res
    else:  # t is number
        if t <= c:
            return part4(t)
        else:
            return a + b * np.sqrt(abs(t - c)) / (d * t) + (n * np.e ** (m * t + z))


def part4(t):
    """Just a constant value"""
    return 0.1207 +0*t  # 0*t to work with ndarrays


def kappa_nagy(t):
    k1 = 30  # part2-part3 smoothing value
    k2 = -50  # part3-part4 smoothing value | negative means use maximum
    k3 = 50  # part1-part2 smoothing value
    if isinstance(t, (list, np.ndarray)):
        y = []
        for i in t:
            # up until intersect between part2 and part3:
            # y.append(smooth_min(part1(i), smooth_min(smooth_min(part2(i), part3(i), k1), part4(i), k2), k3))
            y.append(smooth_min(smooth_min(smooth_min(part2(i), part3(i), k1), part1(i), k3), part4(i), k2))
        return y
    elif isinstance(t, (int, float)):
        return smooth_min(part1(t), smooth_min(smooth_min(part2(t), part3(t), k1), part4(t), k2), k3)
    else:
        raise(TypeError("argument t must be a number or a list/numpy array of numbers"))


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    NAGY_KAPPA15 = "../data/15M_sun.csv"
    #NAGY_KAPPA30 = "../data/30M_sun.csv"
    t_k15, kappa_t15 = np.loadtxt(NAGY_KAPPA15, unpack=True, delimiter=", ")

    times = np.linspace(1, 365, 1000)
    results = kappa_nagy(times)
    print(results)
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(8, 6)
    ax.semilogx(times, part1(times), color="yellow", linestyle='--', label="function 1")
    ax.semilogx(times, part2(times), color="green", linestyle='--', label="function 2")
    ax.semilogx(times, part3(times), color="orange", linestyle='--', label="function 3")
    ax.semilogx(times, part4(times), color="purple", linestyle='--', label="function 4")
    ax.semilogx(times, results, color="red", label="Combined function")
    ax.scatter(t_k15, kappa_t15, color="blue", label="Nagy's data")
    ax.legend()
    ax.set_xlim(1, 400)
    ax.set_ylim(-0.1, 0.9)
    ax.set_xlabel("time (days)")
    ax.set_ylabel("kappa")
    fig.savefig("nagy_func.png")
    plt.show()

