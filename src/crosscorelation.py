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


def create_x(T0, k0, N, SNR):
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


# useless function, depricated
def corr_cosinus(x, s):
    M = np.zeros(len(x) - len(s))
    for l in range(0, len(x) - len(s)):
        denominator = (np.sqrt(np.sum(np.square(x[l : l + len(s)]))) * np.sqrt(np.sum(np.square(s))))
        if denominator != 0:
            M[l] = np.dot(x[l : l + len(s)], s) / denominator
        else:
            M[l] = -float('Inf')
    return M


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


def show_complex_plot(ordinates, titles=None, xlabels=None, ylabels=None, ylim=None, show_complex=False):
    plt.figure()
    for i in range(len(ordinates)):
        xAxis = np.arange(0, len(ordinates[i]))
        # if i == 1: # here only
            # xAxis = np.append(np.zeros(k), np.arrange)
        plot = plt.subplot(len(ordinates) * 100 + 11 + i)
        plot.set_xlim([0, len(ordinates[0])])
        if show_complex:
            plt.plot(xAxis, np.abs(ordinates[i]), 'g', label='abs')
            plt.plot(xAxis, np.real(ordinates[i]), 'r', label='real')
            plt.plot(xAxis, np.imag(ordinates[i]), 'b', label='imag')
            plt.plot(xAxis, np.angle(ordinates[i]), 'y', label='angle')
        else:
            plt.plot(xAxis, ordinates[i])
        if ylim:
            plot.set_ylim(ylim)
        if titles:
            plt.title(titles[i])
        if xlabels:
            plt.xlabel(xlabels[i])
        if ylabels:
            plt.ylabel(ylabels[i])
        plt.grid()

    plt.tight_layout()
    # if SAVE_IMAGE:
        # plt.savefig('../output/img_complex.png')
    plt.show()


if __name__ == '__main__':
    print("Question 2")
# Questions 2.1 and 2.2
    x = create_impulse_signal(100, 50)
    sT = create_impulse_signal(40, 20)
    subx = correlate(x, sT, mode = 'valid')
    print("Maximum correlation in subvector of x with length " + str(len(subx))) # shit here
    print(subx)
    show_complex_plot([x, sT, subx],
                       titles=["Impulse signal of length 100 and peak on 50", "Impulse signal of length 40 and peak on 20", "Correlation"],
                       xlabels=["k", "k", "k"],
                       ylabels=["x[k]", "sT[k]", "C[k]"],
                       ylim=[-0.1, 1.1])

# Question 2.3
    print("E(x) = " + str(energy(x)))
    print("E(sT) = " + str(energy(sT)))

# Question 2.4
    M = correlate(x, sT, mode='valid')
    print("l* = max(M(x_l, sT)) = " + str(np.max(M)))
    print("index of l* = argmax(M(x_l, sT)) = " + str(np.argmax(M)))

    print("Question 3")
# Question 3.4
    sigma = SNR(-5)
    x = create_x(20, 60, 200, sigma)
    period = create_one_period_sin_signal(20)
    C = correlate(x, period, mode = "valid")
    show_complex_plot([x, period, C],
                      titles=["Noisy signal", "One period sinus signal", "Correlation of period and x"],
                      xlabels=["k", "k", "k"],
                      ylabels=["x[k]", "period[k]", "C[k]"])

# Question 3.5
    M = correlate(x / np.sqrt(energy(x)), period / np.sqrt(energy(period)), "valid")
    print("C* = max(M(x_l, sT)) = " + str(np.max(M)))
    print("index of C* = argmax(M(x_l, sT)) = " + str(np.argmax(M)))

    print("Question 4 - from here code is not good now")
# Question 4.1
    img = read_pgm("../input/FindWaldo1.pgm")
    pattern = read_pgm("../input/WadoTarget1.pgm")
    # display_image(img)
    # display_image(pattern)
    corr_2d = correlate2d(img, pattern, mode="same")
    # display_image(corr_2d)
    print(np.max(corr_2d), np.argmax(corr_2d))
    # [y,x] = np.unravel_index(np.argmax(map), corr_2d.map)
    # print(y)
    # print(x)








