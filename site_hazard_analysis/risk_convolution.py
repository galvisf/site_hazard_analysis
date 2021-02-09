from .base import *

def risk_convolution(im_exceedance_frequency, im_list, fragilities):
    medians = fragilities['Median']
    betas = fragilities['Beta']

    freq = im_exceedance_frequency
    dim = im_list[1] - im_list[0]

    dfreq_im = [np.abs((freq[i + 1] - freq[i]) / dim) for i in range(len(im_list) - 1)]
    dfreq_im = np.append(dfreq_im, dfreq_im[-1])

    if type(medians) is not list:
        medians = [medians]
        betas = [betas]

    freq_collapse = np.zeros(len(medians))
    for median, beta, i in zip(medians, betas, range(len(medians))):
        p_collapse_im = stats.lognorm(beta, scale=median).cdf(im_list)
        y = p_collapse_im * dfreq_im
        freq_collapse[i] = np.trapz(y, dx=dim)

    return freq_collapse


def fun_to_fit3Lin(x, k_max, a_min, a, b1, b2):
    y = np.zeros(x.size)

    for i in range(x.size):
        if x[i] < a_min:
            y[i] = k_max
        elif x[i] < a:
            y[i] = k_max + b1 * (x[i] - a_min)
        else:
            y[i] = k_max + b1 * (a - a_min) + b2 * (x[i] - a)

    return y