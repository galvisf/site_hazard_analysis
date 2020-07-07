from .base import *
from .eq_rates import *


def cy_2008_nga(t, m, r_rup, r_jb, r_x, delta, rake_angle, z_tor, f_as, vs30, f_vs30, z_10):
    periods = ['pga', 'pgv', 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5,
               2.0, 3.0, 4.0, 5.0, 7.5, 10.0]

    cy2008_parameters = load_CY2008_parameters(display_tables=False)

    # set up a new function sets all parameters other than t
    cy_2008_nga_partial = partial(cy_2008_nga_sub, m, r_rup, r_jb, r_x,
                                  delta, rake_angle, z_tor,
                                  f_as, vs30, f_vs30, z_10,
                                  cy2008_parameters)

    # check the t input
    if type(t) == str:
        # 'pga' or 'pgv' is a viable string
        if t in ['pga', 'pgv']:
            t = [t]
        # any other string will return all the T values from the model
        else:
            t = periods[2:]  # remove pga and pgv
    # if t is a scalar value, put it in list form
    if (type(t) == float) | (type(t) == int):
        t = [t]

    # if all the desired T values are available in the model, run it directly
    if all(t_i in periods for t_i in t):
        gmpe_results = [cy_2008_nga_partial(t_i) for t_i in t]
        median = [x[0] for x in gmpe_results]
        sigma = [x[1] for x in gmpe_results]

    else:
        # check that all the periods are within the period range
        periods = np.array(periods[2:])
        if all((t >= min(periods)) & (t <= max(periods))):

            # if less than n periods need to be interplolated
            n_periods = 5
            if n_periods >= np.sum([t_i not in periods for t_i in t]):

                median = np.zeros(len(t))
                sigma = np.zeros_like(median)
                # loop over periods
                for i in range(len(t)):
                    t_i = t[i]

                    # if period does not need to be interpolated
                    if t_i in periods:
                        median[i], sigma[i] = cy_2008_nga_partial(t_i)

                    # if period does to be interpolated
                    else:
                        idx = np.searchsorted(periods, t_i)

                        t_low = periods[idx - 1]
                        median_low, sigma_low = cy_2008_nga_partial(t_low)

                        t_hi = periods[idx]
                        median_hi, sigma_hi = cy_2008_nga_partial(t_hi)

                        # interpolate in log space
                        log_ti = np.log(t_i)
                        x = np.log([t_low, t_hi])
                        medians = np.log([median_low, median_hi])
                        sigmas = [sigma_low, sigma_hi]

                        median[i] = np.exp(np.interp(log_ti, x, medians))
                        sigma[i] = np.interp(log_ti, x, sigmas)

            # if more than n periods needs to be interpolated
            else:
                # retrieve only the periods that are necessary for interpolation
                idx_min = max(np.searchsorted(periods, min(t)) - 1, 0)
                idx_max = np.searchsorted(periods, max(t)) + 1
                periods = periods[idx_min:idx_max]

                # retrieve all the available T values
                gmpe_results = [cy_2008_nga_partial(t_i) for t_i in periods]
                median = [x[0] for x in gmpe_results]
                sigma = [x[1] for x in gmpe_results]

                # interpolate in log space
                median = np.exp(np.interp(np.log(t), np.log(periods), np.log(median)))
                sigma = np.interp(np.log(t), np.log(periods), sigma)

        else:
            raise ValueError('Periods do not fall within the GMPE range')

    return median, sigma, t


def load_CY2008_parameters(display_tables):
    periods = ['pga', 'pgv', 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5,
               2.0, 3.0, 4.0, 5.0, 7.5, 10.0]

    c1 = [-1.2687, 2.2884, -1.2687, -1.2515, -1.1744, -1.0671, -0.9464, -0.7051, -0.5747, -0.5309, -0.6352, -0.7766,
          -0.9278, -1.2176, -1.4695, -1.9278, -2.2453, -2.7303, -3.1413, -3.7413, -4.1814, -4.5178, -5.1224, -5.5872]
    c1a = [0.1, 0.1094, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.0999, 0.0997, 0.0991, 0.0936, 0.0766,
           0.0022, -0.0591, -0.0931, -0.0982, -0.0994, -0.0999, -0.1000]
    c1b = [-0.2550, -0.0626, -0.2550, -0.2550, -0.2550, -0.2550, -0.2550, -0.2540, -0.2530, -0.2500, -0.2449, -0.2382,
           -0.2313, -0.2146, -0.1972, -0.1620, -0.1400, -0.1184, -0.1100, -0.1040, -0.1020, -0.1010, -0.1010, -0.1000]
    c2 = 1.06
    c3 = 3.45
    cn = [2.996, 1.648, 2.996, 3.292, 3.514, 3.563, 3.547, 3.448, 3.312, 3.044, 2.831, 2.658, 2.505, 2.261, 2.087,
          1.812, 1.648, 1.511, 1.470, 1.456, 1.465, 1.478, 1.498, 1.502]
    cm = [4.1840, 4.2979, 4.1840, 4.1879, 4.1556, 4.1226, 4.1011, 4.0860, 4.1030, 4.1717, 4.2476, 4.3184, 4.3844,
          4.4979, 4.5881, 4.7571, 4.8820, 5.0697, 5.2173, 5.4385, 5.5977, 5.7276, 5.9891, 6.1930]
    c4 = -2.1
    c4a = -0.5
    crb = 50
    c5 = [6.1600, 5.1760, 6.1600, 6.1580, 6.1550, 6.1508, 6.1441, 6.1200, 6.0850, 5.9871, 5.8699, 5.7547, 5.6527,
          5.4997, 5.4029, 5.2900, 5.2480, 5.2194, 5.2099, 5.2040, 5.2020, 5.2010, 5.2000, 5.2000]
    c6 = [0.4893, 0.4407, 0.4893, 0.4892, 0.4890, 0.4888, 0.4884, 0.4872, 0.4854, 0.4808, 0.4755, 0.4706, 0.4665,
          0.4607, 0.4571, 0.4531, 0.4517, 0.4507, 0.4504, 0.4501, 0.4501, 0.4500, 0.4500, 0.4500]
    chm = 3.0
    c7 = [0.0512, 0.0207, 0.0512, 0.0512, 0.0511, 0.0508, 0.0504, 0.0495, 0.0489, 0.0479, 0.0471, 0.0464, 0.0458,
          0.0445, 0.0429, 0.0387, 0.0350, 0.0280, 0.0213, 0.0106, 0.0041, 0.0010, 0.0000, 0.0000]
    c7a = [0.0860, 0.0437, 0.0860, 0.0860, 0.0860, 0.0860, 0.0860, 0.0860, 0.0860, 0.0860, 0.0860, 0.0860, 0.0860,
           0.0850, 0.0830, 0.0690, 0.0450, 0.0134, 0.0040, 0.0010, 0.0000, 0.0000, 0.0000, 0.0000]
    c9 = [0.7900, 0.3079, 0.7900, 0.8129, 0.8439, 0.8740, 0.8996, 0.9442, 0.9677, 0.9660, 0.9334, 0.8946, 0.8590,
          0.8019, 0.7578, 0.6788, 0.6196, 0.5101, 0.3917, 0.1244, 0.0086, 0.0000, 0.0000, 0.0000]
    c9a = [1.5005, 2.6690, 1.5005, 1.5028, 1.5071, 1.5138, 1.5230, 1.5597, 1.6104, 1.7549, 1.9157, 2.0709, 2.2005,
           2.3886, 2.5000, 2.6224, 2.6690, 2.6985, 2.7085, 2.7145, 2.7164, 2.7172, 2.7177, 2.7180]
    c10 = [-0.3218, -0.1166, -0.3218, -0.3323, -0.3394, -0.3453, -0.3502, -0.3579, -0.3604, -0.3565, -0.3470, -0.3379,
           -0.3314, -0.3256, -0.3189, -0.2702, -0.2059, -0.0852, 0.0160, 0.1876, 0.3378, 0.4579, 0.7514, 1.1856]
    cy1 = [-0.00804, -0.00275, -0.00804, -0.00811, -0.00839, -0.00875, -0.00912, -0.00973, -0.00975, -0.00883, -0.00778,
           -0.00688, -0.00612, -0.00498, -0.00420, -0.00308, -0.00246, -0.00180, -0.00147, -0.00117, -0.00107, -0.00102,
           -0.00096, -0.00094]
    cy2 = [-0.00785, -0.00625, -0.00785, -0.00792, -0.00819, -0.00855, -0.00891, -0.00950, -0.00952, -0.00862, -0.00759,
           -0.00671, -0.00598, -0.00486, -0.00410, -0.00301, -0.00241, -0.00176, -0.00143, -0.00115, -0.00104, -0.00099,
           -0.00094, -0.00091]
    cy3 = 4.0
    phi1 = [-0.4417, -0.7861, -0.4417, -0.4340, -0.4177, -0.4000, -0.3903, -0.4040, -0.4423, -0.5162, -0.5697, -0.6109,
            -0.6444, -0.6931, -0.7246, -0.7708, -0.7990, -0.8382, -0.8663, -0.9032, -0.9231, -0.9222, -0.8346, -0.7332]
    phi2 = [-0.1417, -0.0699, -0.1417, -0.1364, -0.1403, -0.1591, -0.1862, -0.2538, -0.2943, -0.3113, -0.2927, -0.2662,
            -0.2405, -0.1975, -0.1633, -0.1028, -0.0699, -0.0425, -0.0302, -0.0129, -0.0016, 0.0000, 0.000, 0.000]
    phi3 = [-0.007010, -0.008444, -0.007010, -0.007279, -0.007354, -0.006977, -0.006467, -0.005734, -0.005604,
            -0.005845, -0.006141, -0.006439, -0.006704, -0.007125, -0.007435, -0.008120, -0.008444, -0.007707,
            -0.004792, -0.001828, -0.001523, -0.01440, -0.001369, -0.001361]
    phi4 = [0.102151, 5.41000, 0.102151, 0.108360, 0.119888, 0.133641, 0.148927, 0.190596, 0.230662, 0.266468, 0.255253,
            0.231541, 0.207277, 0.165464, 0.133828, 0.085153, 0.058595, 0.031787, 0.019716, 0.009643, 0.005379,
            0.003223, 0.001134, 0.000515]
    phi5 = [0.2289, 0.2899, 0.2289, 0.2289, 0.2289, 0.2289, 0.2290, 0.2292, 0.2297, 0.2326, 0.2386, 0.2497, 0.2674,
            0.3120, 0.3610, 0.4353, 0.4629, 0.4756, 0.4785, 0.4796, 0.4799, 0.4799, 0.4800, 0.4800]
    phi6 = [0.014996, 0.006718, 0.014996, 0.014996, 0.014996, 0.014996, 0.014996, 0.014996, 0.014996, 0.014988,
            0.014964, 0.014881, 0.014639, 0.013493, 0.011133, 0.006739, 0.005749, 0.005544, 0.005521, 0.005517,
            0.005517, 0.005517, 0.005517, 0.005517]
    phi7 = [580.0, 459.0, 580.0, 580.0, 580.0, 579.9, 579.9, 579.6, 579.2, 577.2, 573.9, 568.5, 560.5, 540.0, 512.9,
            441.9, 391.8, 348.1, 332.5, 324.1, 321.7, 320.9, 320.3, 320.1]
    phi8 = [0.0700, 0.1138, 0.0700, 0.0699, 0.0701, 0.0702, 0.0701, 0.0686, 0.0646, 0.0494, -0.0019, -0.0479, -0.0756,
            -0.0960, -0.0998, -0.0765, -0.0412, 0.0140, 0.0544, 0.1232, 0.1859, 0.2295, 0.2660, 0.2682]
    tau1 = [0.3437, 0.2539, 0.3437, 0.3471, 0.3603, 0.3718, 0.3848, 0.3878, 0.3835, 0.3719, 0.3601, 0.3522, 0.3438,
            0.3351, 0.3353, 0.3429, 0.3577, 0.3769, 0.4023, 0.4406, 0.4784, 0.5074, 0.5328, 0.5542]
    tau2 = [0.2637, 0.2381, 0.2637, 0.2671, 0.2803, 0.2918, 0.3048, 0.3129, 0.3152, 0.3128, 0.3076, 0.3047, 0.3005,
            0.2984, 0.3036, 0.3205, 0.3419, 0.3703, 0.4023, 0.4406, 0.4784, 0.5074, 0.5328, 0.5542]
    sigma1 = [0.4458, 0.4496, 0.4458, 0.4458, 0.4535, 0.4589, 0.4630, 0.4702, 0.4747, 0.4798, 0.4816, 0.4815, 0.4801,
              0.4758, 0.4710, 0.4621, 0.4581, 0.4493, 0.4459, 0.4433, 0.4424, 0.4420, 0.4416, 0.4414]
    sigma2 = [0.3459, 0.3554, 0.3459, 0.3459, 0.3537, 0.3592, 0.3635, 0.3713, 0.3769, 0.3847, 0.3902, 0.3946, 0.3981,
              0.4036, 0.4079, 0.4157, 0.4213, 0.4213, 0.4213, 0.4213, 0.4213, 0.4213, 0.4213, 0.4213]
    sigma3 = [0.8000, 0.7504, 0.8000, 0.8000, 0.8000, 0.8000, 0.8000, 0.8000, 0.8000, 0.8000, 0.8000, 0.7999, 0.7997,
              0.7988, 0.7966, 0.7792, 0.7504, 0.7136, 0.7035, 0.7006, 0.7001, 0.7000, 0.7000, 0.7000]
    sigma4 = [0.0663, 0.0133, 0.0663, 0.0663, 0.0663, 0.0663, 0.0663, 0.0663, 0.0663, 0.0612, 0.0530, 0.0457, 0.0398,
              0.0312, 0.0255, 0.0175, 0.0133, 0.0090, 0.0068, 0.0045, 0.0034, 0.0027, 0.0018, 0.0014]

    # table 1
    c_constants = dict()
    c_constants['c2'] = c2
    c_constants['c3'] = c3
    c_constants['c4'] = c4
    c_constants['c4a'] = c4a
    c_constants['crb'] = crb
    c_constants['chm'] = chm
    c_constants['cy3'] = cy3

    # table 2
    df = pd.DataFrame(index=periods)
    df.index.name = 'Spectra Period (sec)'
    df['c1'] = c1
    df['c1a'] = c1a
    df['c1b'] = c1b
    df['cn'] = cn
    df['cm'] = cm
    df['c5'] = c5
    df['c6'] = c6
    df['c7'] = c7
    df['c7a'] = c7a
    df['c9'] = c9
    df['c9a'] = c9a
    df['c10'] = c10
    df['cy1'] = cy1
    df['cy2'] = cy2
    c_dict = df.to_dict()

    # table 3
    df = pd.DataFrame(index=periods)
    df.index.name = 'Spectra Period (sec)'
    df['phi1'] = phi1
    df['phi2'] = phi2
    df['phi3'] = phi3
    df['phi4'] = phi4
    df['phi5'] = phi5
    df['phi6'] = phi6
    df['phi7'] = phi7
    df['phi8'] = phi8
    phi_dict = df.to_dict()

    # table 4
    df = pd.DataFrame(index=periods)
    df.index.name = 'Spectra Period (sec)'
    df['tau1'] = tau1
    df['tau2'] = tau2
    df['sigma1'] = sigma1
    df['sigma2'] = sigma2
    df['sigma3'] = sigma3
    df['sigma4'] = sigma4
    variance_dict = df.to_dict()

    if display_tables:
        display(pd.DataFrame(c_constants, index=[0]))
        display(pd.DataFrame(c_dict))
        display(pd.DataFrame(phi_dict))
        display(pd.DataFrame(variance_dict))

    return [c_constants, c_dict, phi_dict, variance_dict]



def cy_2008_nga_sub(m, r_rup, r_jb, r_x, delta, rake_angle, z_tor, f_as, vs30, f_vs30, z_10, cy2008_parameters, i_period):
    delta = delta * np.pi / 180  # convert to radians
    f_rv = (rake_angle >= 30) & (rake_angle <= 150)  # reverse fault flag, 1 for lambda between 30 and 150, 0 otherwise
    f_nm = (rake_angle >= -120) & (
    rake_angle <= -60)  # normal fault flag,  1 for lambda between -120 and -60, 0 otherwise
    f_hw = r_x >= 0  # hanging wall flag

    f_inferred = (f_vs30 == 1)  # 1 if inferred
    f_measured = (f_vs30 == 0)  # 1 if measured

    ## Evaluate Equation 13a for the reference intensity at Vs30=1130

    # Retrieve the coefficient that do not depend on the period of interest
    c_constants = cy2008_parameters[0]
    c2  = c_constants['c2']
    c3  = c_constants['c3']
    c4  = c_constants['c4']
    c4a = c_constants['c4a']
    crb = c_constants['crb']
    chm = c_constants['chm']
    cy3 = c_constants['cy3']

    # Retrieve the coefficients that depend on the period of interest
    c_dict = cy2008_parameters[1]
    c1 = c_dict['c1'][i_period]
    c1a = c_dict['c1a'][i_period]
    c1b = c_dict['c1b'][i_period]
    cn = c_dict['cn'][i_period]
    cm = c_dict['cm'][i_period]
    c5 = c_dict['c5'][i_period]
    c6 = c_dict['c6'][i_period]
    c7 = c_dict['c7'][i_period]
    c7a = c_dict['c7a'][i_period]
    c9 = c_dict['c9'][i_period]
    c9a = c_dict['c9a'][i_period]
    c10 = c_dict['c10'][i_period]
    cy1 = c_dict['cy1'][i_period]
    cy2 = c_dict['cy2'][i_period]

    term1 = c1

    # reverse faulting, normal faulting, depth to top of rupture
    # AfterShock flag selects either term2 or term3, other term goes to zero
    term2 = (c1a * f_rv + c1b * f_nm + c7 * (z_tor - 4)) * (1 - f_as)
    term3 = (c10 + c7a * (z_tor - 4)) * f_as

    # magnitude and distance
    term4 = c2 * (m - 6)
    term5 = ((c2 - c3) / cn) * np.log(1 + np.exp(cn * (cm - m)))
    term6 = c4 * np.log(r_rup + c5 * np.cosh(c6 * max(m - chm, 0)))
    term7 = (c4a - c4) * np.log(np.sqrt(r_rup ** 2 + crb ** 2))
    term8 = (cy1 + cy2 / (np.cosh(max(m - cy3, 0)))) * r_rup
    term9 = c9 * f_hw * np.tanh((r_x * np.cos(delta) ** 2) / c9a) * (
    1 - np.sqrt(r_jb ** 2 + z_tor ** 2) / (r_rup + 0.001))

    y_ref = np.exp(term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9)

    ## Apply soil effects via Equation 13b

    # Retrieve the coefficients that depend on the period of interest
    phi_dict = cy2008_parameters[2]
    phi1 = phi_dict['phi1'][i_period]
    phi2 = phi_dict['phi2'][i_period]
    phi3 = phi_dict['phi3'][i_period]
    phi4 = phi_dict['phi4'][i_period]
    phi5 = phi_dict['phi5'][i_period]
    phi6 = phi_dict['phi6'][i_period]
    phi7 = phi_dict['phi7'][i_period]
    phi8 = phi_dict['phi8'][i_period]

    term1 = phi1 * min(np.log(vs30 / 1130), 0)
    term2 = phi2 * (np.exp(phi3 * (min(vs30, 1130) - 360)) - np.exp(phi3 * (1130 - 360))) * np.log(
        (y_ref + phi4) / phi4)
    term3 = phi5 * (1 - (1 / np.cos(phi6 * max(0, z_10 - phi7))))
    term4 = phi8 / (np.cosh(0.15 * max(0, z_10 - 15)))

    y = np.exp(np.log(y_ref) + term1 + term2 + term3 + term4)

    ## Evaluate the variability

    # Retrieve the coefficients that depend on the period of interest
    variance_dict = cy2008_parameters[3]
    tau1 = variance_dict['tau1'][i_period]
    tau2 = variance_dict['tau2'][i_period]
    sigma1 = variance_dict['sigma1'][i_period]
    sigma2 = variance_dict['sigma2'][i_period]
    sigma3 = variance_dict['sigma3'][i_period]
    sigma4 = variance_dict['sigma4'][i_period]

    # b and c, from Equation 10
    b = phi2 * (np.exp(phi3 * (min(vs30, 1130) - 360)) - np.exp(phi3 * (1130 - 360)))
    c = phi4
    nl0 = b * (y_ref / (y_ref + c))

    # Equation 19
    tau = tau1 + (tau2 - tau1) / 2 * (min(max(m, 5), 7) - 5)
    # Equation 20
    sigma = (sigma1 + (sigma2 - sigma1) / 2 * (min(max(m, 5), 7) - 5) + sigma4 * f_as) * \
            (np.sqrt((sigma3 * f_inferred + 0.7 * f_measured) + (1 + nl0) ** 2))
    # Equation 21
    sigma = np.sqrt((1 + nl0) ** 2 * tau ** 2 + sigma ** 2)

    return y, sigma


def baker_jayaram_correlation(t1, t2):
    """
    Computes the correlation of epsilons for the NGA ground motion models

    Created by Anne Hulsey, 2/20/2019,
        modeled after Matlab code by Jack Baker, available at:
        (https://web.stanford.edu/~bakerjw/GMPEs/baker_jayaram_correlation.m)


    The function is strictly empirical, fitted over the range 0.01s <= T1, T2 <= 10s

    Documentation is provided in the following document:
    Baker, J.W. and Jayaram, N. (2008), "Correlation of spectral acceleration
    values from NGA ground motion models," Earthquake Spectra, 24 (1), 299-317.

    Parameters
    ----------
    t1, t2: float
        the two periods of interest. The periods may be equal,
        with no restriction on which one is larger

    Returns
    -------
    rho: float
        the predicted correlation coefficent between the two periods
    """

    # order the input periods
    t_min = min(t1, t2)
    t_max = max(t1, t2)

    # calculate c1
    c1 = 1 - np.cos((np.pi / 2) - 0.366 * np.log(t_max / (max(t_min, 0.109))))

    # calculate c2
    if t_max < 0.2:
        c2 = 1 - 0.105 * (1 - 1 / (1 + np.exp(100 * t_max - 5))) * ((t_max - t_min) / (t_max - 0.0099))
    else:
        c2 = 0

    # calculate c3
    if t_max < 0.109:
        c3 = c2
    else:
        c3 = c1

    # calculate c4
    c4 = c1 + 0.5 * (np.sqrt(c3) - c3) * (1 + np.cos((np.pi * t_min) / 0.109))

    # calculate rho
    if t_max < 0.109:
        rho = c2
    elif t_min > 0.109:
        rho = c1
    elif t_max < 0.2:
        rho = min(c2, c4)
    else:
        rho = c4

    return rho


def set_sa_avg_parameters(t, t_min_factor, t_max_factor, dt):
    sa_avg_parameters = dict()

    # create a vector of period values
    t_min = t_min_factor * t
    t_max = t_max_factor * t
    t = np.arange(t_min, t_max, dt)
    if t[-1] != t_max:
        t = np.append(t, t_max)
    sa_avg_parameters['t_range'] = t

    # calculate the empirical correlation between each period pair
    n_t = len(t)
    rho = np.zeros((n_t, n_t))
    for i_t1, t1 in enumerate(t):
        for i_t2, t2 in enumerate(t[t >= t1]):
            rho[i_t1, i_t2 + i_t1] = baker_jayaram_correlation(t1, t2)

    # mirror the upper triangle to make a symmetric correlation matrix
    i_lower = np.tril_indices(n_t, -1)
    rho[i_lower] = rho.T[i_lower]
    sa_avg_parameters['rho'] = rho

    return sa_avg_parameters


def sa_avg_gmpe(t, sa_avg_parameters, m, r_rup, r_jb, r_x, delta, rake_angle, z_tor, f_as, vs30, f_vs30, z_10):

    # if the sa_avg_parameters were not passed, set them based on Eads 2015
    if sa_avg_parameters is None:
        t_min_factor = 0.2
        t_max_factor = 3
        dt = 0.01
        sa_avg_parameters = set_sa_avg_parameters(t, t_min_factor, t_max_factor, dt)

    t = sa_avg_parameters['t_range']
    rho = sa_avg_parameters['rho']

    # retrieve median and standard deviations for each period
    y, sigma, t = cy_2008_nga(t, m, r_rup, r_jb, r_x, delta, rake_angle, z_tor, f_as, vs30, f_vs30, z_10)

    # calculate the median for sa_avg
    n_t = len(t)
    y_sa_avg = np.exp((1 / n_t) * np.sum(np.log(y)))

    # calculate the sigma for sa_avg
    sigma_sa_avg = np.sqrt((1/n_t**2) * np.matmul(np.matmul(sigma, rho), sigma))

    return y_sa_avg, sigma_sa_avg


def prob_exceeding_im(t, m, r, f_as, vs30, f_vs30, z_10, im_type, sa_avg_parameters, im_list):

    # set all distances equal to r
    r_rup = r  # closest distance to rupture plane (km)
    r_jb = r  # Joyner-Boore distance to rupture plane (km)
    r_x = r  # Site coordinate perpendicular to fault strike (km)

    # assume a vertical fault
    delta = 90  # fault dip angle

    # assume rupture reaches the surface
    z_tor = 0  # depth to top of rupture (km)

    # assume strike-slip (see Chiou and Youngs 2008 for classification of normal and reverse)
    rake_angle = 0  # fault rake angle ('lambda' in Chiou and Youngs 2008)

    # if z_10 is not specified, use Eqn 1 from Chiou and Youngs 2008
    # depth to shear wave velocity of 1 km/s (m)
    if z_10 is None:
        z_10 = np.exp(28.5 - (3.82 / 8) * np.log(vs30 ** 8 + 378.7 ** 8))

    # call the gmpe based on the intensity measure
    if im_type == 'sa_t':
        y, sigma, t = cy_2008_nga(t, m, r_rup, r_jb, r_x, delta, rake_angle, z_tor, f_as, vs30, f_vs30, z_10)
        y = y[0]
        sigma = sigma[0]
    elif im_type == 'sa_avg':
        y, sigma = sa_avg_gmpe(t, sa_avg_parameters, m, r_rup, r_jb, r_x, delta, rake_angle, z_tor, f_as, vs30, f_vs30, z_10)

    # calculate the probability of exceeding each im value
    p = 1 - stats.lognorm.cdf(im_list, s=sigma, scale=y)
    epsilon = (np.log(im_list) - np.log(y)) / sigma

    return p, epsilon


def hazard_curve_integration(t, rate_r_m, im_type, sa_avg_parameters, im_list, site_parameters, f_as, im_plot):
    # retrieve the list of m and r values in the rate matrix
    m_list = rate_r_m.index
    r_list = rate_r_m.columns
    rate_r_m = rate_r_m.to_numpy()

    # retrieve relevant site parameters
    vs30 = site_parameters['vs30']
    f_vs30 = site_parameters['f_vs30']
    z_10 = site_parameters['z_10']

    # initialize an array for storing the probabilities of exceedance
    prob_im_given_r_m = np.zeros([len(m_list), len(r_list), len(im_list)])

    # loop over all m and r values for the p(exceedance) of each im level
    for m, i_m in zip(m_list, range(len(m_list))):
        for r, i_r in zip(r_list, range(len(r_list))):
            p_im, epsilon_im = prob_exceeding_im(t, m, r, f_as, vs30, f_vs30, z_10, im_type, sa_avg_parameters, im_list)
            prob_im_given_r_m[i_m, i_r, :] = p_im

    # combine the rate of each m and r pair with the p(exceedance)
    rate_im_given_r_m = np.repeat(np.expand_dims(rate_r_m, 2), len(im_list), 2) * prob_im_given_r_m
    rate_im = np.sum(np.sum(rate_im_given_r_m, axis=0), axis=0)

    # find the probability of each m and r pair, given the im was exceeded
    prob_m_r_given_im = np.zeros_like(rate_im_given_r_m)
    for i_im in range(len(im_list)):
        prob_m_r_given_im[:, :, i_im] = rate_im_given_r_m[:, :, i_im] / rate_im[i_im]

    # plot figures for an im level, if desired
    if im_plot is not None:
        i_im = np.where(im_list == im_plot)[0][0]

        # set the intensity measure label
        if im_type == 'sa_t':
            label_tag = ')'
        elif im_type == 'sa_avg':
            label_tag = ')$_{avg}$'
        label_tag = 'Sa(T=' + '{0:.1f}'.format(t) + label_tag + ' >= ' + '{0:.2g}g'.format(im_list[i_im])

        # probability of exceedance, given m and r pairs
        fig, ax = plt.subplots(1, 1, figsize=(15, 5))
        _ = sns.heatmap(prob_im_given_r_m[:, :, i_im], linewidth=0.5)
        _ = ax.set_xlabel('Site to Source Distance, R')
        _ = ax.set_xticks(np.arange(len(r_list)))
        _ = ax.set_xticklabels(r_list, rotation=90)
        _ = ax.set_ylabel('Magnitude, M')
        _ = ax.set_yticklabels(m_list, rotation=0)
        _ = ax.set_title('P(' + label_tag + ') | M,R')
        _ = plt.show()

        # contribution of each m and r pair to the current exceedance level
        fig, ax = plt.subplots(1, 1, figsize=(15, 5))
        _ = sns.heatmap(prob_m_r_given_im[:, :, i_im], linewidth=0.5)
        _ = ax.set_xlabel('Site to Source Distance, R')
        _ = ax.set_xticks(np.arange(len(r_list)))
        _ = ax.set_xticklabels(r_list, rotation=90)
        _ = ax.set_ylabel('Magnitude, M')
        _ = ax.set_yticklabels(m_list, rotation=0)
        _ = ax.set_title('Contribution to ' + label_tag)
        _ = plt.show()

        # mean annual frequency of exceeding the current im level from each m and r pair
        fig, ax = plt.subplots(1, 1, figsize=(15, 5))
        _ = sns.heatmap(rate_im_given_r_m[:, :, i_im], linewidth=0.5)
        _ = ax.set_xlabel('Site to Source Distance, R')
        _ = ax.set_xticks(np.arange(len(r_list)))
        _ = ax.set_xticklabels(r_list, rotation=90)
        _ = ax.set_ylabel('Magnitude, M')
        _ = ax.set_yticklabels(m_list, rotation=0)
        _ = ax.set_title(
            'MAF of ' + label_tag + ' | M,R')
        _ = plt.show()

    return rate_im, rate_im_given_r_m


def aftershock_im_exceedance_probability(m_mainshock, aftershock_parameters, site_fault_geometry, site_parameters,
                                         bldg_t, im_type, sa_avg_parameters, im_list, plot_flag):
    f_as = 1

    m_max = m_mainshock
    m_min = aftershock_parameters['m_min']
    b = aftershock_parameters['b']
    dm = 0.25
    dr = 5

    r_m_distribution = aftershock_r_and_m_distribution(site_fault_geometry, m_min, m_max, b, dm, dr, plot_flag=False)

    im_plot = None
    as_im_probability, prob_im_given_r_m = hazard_curve_integration(bldg_t, r_m_distribution, im_type,
                                                                    sa_avg_parameters, im_list, site_parameters, f_as,
                                                                    im_plot)

    if plot_flag:
        if im_type == 'sa_t':
            label_tag = ')'
        elif im_type == 'sa_avg':
            label_tag = ')$_{avg}$'

        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        _ = ax.loglog(im_list, as_im_probability)

        _ = ax.set_xlabel('Sa(T=' + '{0:.1f}'.format(bldg_t) + label_tag + ' [g]')
        _ = ax.set_ylabel('Probability of Exceedance,\nGiven the Occurance\nof an Aftershock')
        _ = ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y, _: '{:.2g}'.format(y)))
        _ = ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y, _: '{:.2g}'.format(y)))

        dim = min(im_list)
        _ = ax.set_ylim([0.00001, 1])
        _ = ax.set_xlim([dim, 1])
        _ = ax.set_xticks([dim * 5, 0.005, 0.05, 0.5], minor=True)
        _ = ax.grid(axis='y', which='major')
        _ = ax.grid(axis='x', which='both')

        _ = ax.set_title('Mainshock Magnitude: ' + str(m_max) + 'Mw')

    return as_im_probability
