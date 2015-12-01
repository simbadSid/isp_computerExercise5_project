import numpy as np
from scipy.signal import correlate, correlate2d
import matplotlib.pyplot as plt
import re


def create_impulse_signal(N, a):
    sig = np.zeros(N)
    sig[a] = 1
    return sig


def energy(signal):
    e = 0
    for s in signal:
        e += s ** 2
    return e


def create_one_period_sin_signal(T0):
    x = np.arange(0, T0)
    return np.sin(2 * np.pi * x / T0)


def p(signal, T0):
    return energy(signal) / T0


def SNR(SNR_DB):
    return 10 ** (SNR_DB / 10.)


def createX(T0, k0, N, SNR):
    sig = create_one_period_sin_signal(T0)
    var = p(sig, T0) / SNR
    noise = np.random.normal(0, scale=var, size=N)
    x = np.zeros(N)
    for k in range(N):
        if k - k0 >= 0 and k - k0 < T0:
            x[k] = sig[k - k0] + noise[k]
        else:
            x[k] = noise[k]
    return x


def correlation(u, v):
    # function assumes that that len(v) < len(u)
    # if not we swap them
    if len(v) > len(u):
        temp = u
        u = v
        v = temp
    corr = np.zeros(len(u) + len(v) - 1)
    for k in range(len(u) + len(v) - 1):
        val = 0
        if k < len(v):
            for i in range(k + 1):
                val += u[i] * v[k - i]
        if k >= len(v) and k < len(u):
            for i in range(k - len(v) + 1, k + 1):
                val += u[i] * v[k - i]
        if k >= len(u):
            for i in range(k - len(v) + 1, len(u)):
                val += u[i] * v[k - i]
        corr[k] = val
    return corr


def M_corr(x, s):
    M = np.zeros(len(x) - len(s))
    # max_val = -1
    # max_ind = 0
    for l in range(0, len(x) - len(s)):
        val = np.dot(x[l : l + len(s)], s) / (np.sqrt(np.sum(np.square(x[l : l + len(s)]))) * np.sqrt(np.sum(np.square(s))))
        # if M > max_val:
            # max_val = M
            # max_ind = l
        # print(M)
        M[l] = val
    # return (max_val, max_ind)
    return M


def max_correlation(C):
    max = -1
    k = 0
    for i in range(len(C)):
        if C[i] > max:
            max = C[i]
            k = i
    return (max, k)


def read_pgm(filename, byteorder='>'):
    """Return image data from a raw PGM file as numpy array.
 
    Format specification: http://netpbm.sourceforge.net/doc/pgm.html
 
    """
    with open(filename, 'rb') as f:
        buffer = f.read()
    try:
        header, width, height, maxval = re.search(
            b"(^P5\s(?:\s*#.*[\r\n])*"
            b"(\d+)\s(?:\s*#.*[\r\n])*"
            b"(\d+)\s(?:\s*#.*[\r\n])*"
            b"(\d+)\s(?:\s*#.*[\r\n]\s)*)", buffer).groups()
    except AttributeError:
        raise ValueError("Not a raw PGM file: '%s'" % filename)
    return np.frombuffer(buffer,
                            dtype='u1' if int(maxval) < 256 else byteorder+'u2',
                            count=int(width)*int(height),
                            offset=len(header)
                            ).reshape((int(height), int(width)))


def display_plot(x, y, title="Plot", xlabel="x", ylabel="y", limits=[], legend=False):
    plt.figure()
    if isinstance(y, list):
        for i in range(len(y)):
            plt.plot(x, y[i], lebel=str(i))
    else:
        plt.plot(x, y)
    plt.grid()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if legend:
        plt.legend()
    if len(limits) > 0:
        plt.axis(limits)
    # if SHOW_TITLE:
    plt.title(title)
    # if SAVE_IMAGE:
        # plt.savefig('../output/' + title.replace(' ', '_') + '.png')
    plt.show()


def display_image(image):
    plt.imshow(image, plt.cm.gray, vmin=0, vmax=255)
    # if SAVE_IMAGE:
        # plt.savefig('../output/img.png')
    plt.show()

if __name__ == '__main__':
# Questions 2.1 and 2.2
    x = create_impulse_signal(100, 50)
    sT = create_impulse_signal(40, 20)
    subx = correlate(sT, x, mode = 'same')
    print("Maximum correlation in subvector of x with length " + str(len(subx)))
    print(subx)

# Question 2.3
    print("E(x) = " + str(energy(x)))
    print("E(sT) = " + str(energy(sT)))

# Question 3.4
    val = 10 * np.log(5) / np.log(10)
    # val = SNR(-5)
    x = createX(20, 60, 200, val)
    period = create_one_period_sin_signal(20)
    # C = correlation(x, period)
    C = correlate(x, period, mode = "same")
    # display_plot(range(len(x)), x, title="Noisy signal", xlabel="k", ylabel="x[k]")
    # display_plot(range(len(period)), period, title="One period sinus signal", xlabel="k", ylabel="period[k]")
    # display_plot(range(len(C)), C, title="Correlation of period and x", xlabel="k", ylabel="C[k]")

# Question 3.5
    max_val = np.max(C)
    max_index = np.argmax(C)
    # (max_val, max_index) = max_correlation(C)
    print(max_val, max_index) # WTF?

    M = M_corr(x, period)
    max_val = np.max(M)
    max_index = np.argmax(M)
    print(max_val, max_index)

# Question 4.1
    img = read_pgm("../input/FindWaldo1.pgm")
    pattern = read_pgm("../input/WadoTarget1.pgm")
    display_image(img)
    display_image(pattern)
    corr_2d = correlate2d(img, pattern, mode="same")
    display_image(corr_2d)
    print(np.max(corr_2d), np.argmax(corr_2d))
    # [y,x] = np.unravel_index(np.argmax(map), corr_2d.map)
    # print(y)
    # print(x)








