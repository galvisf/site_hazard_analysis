from .base import *


def jprint(obj):
    # create a formatted string of the Python JSON object
    text = json.dumps(obj, sort_keys=True, indent=4)
    print(text)


def getHCfromUSGS(parameters):
    # INPUT
    #   parameters = dict with info of the ste
    #   parameters = {
    #                   "edition":  ["E2008", "E2014"]
    #                   "region":   ["COUS", "WUS", "CEUS"]
    #                   "latitude": XXX,
    #                   "longitude": XXXX,
    #                   "imt": (intensity measure type) [PGA, SA0P2, SA1P0]
    #                   "vs30": [180, 259, 360, 537, 760, 1150, 2000] m/s
    #                }
    #
    # OUTPUTS
    #   hc  = dict with the hazard curve
    #         ims = list with values of the IM of interest
    #         maf = Mean annual frequency of exceedence of each possible value of IM
    #

    url = "https://earthquake.usgs.gov/nshmp-haz-ws/hazard"

    status_code = 'busy'
    while status_code == 'busy':
        response = requests.get(url, params=parameters)
        response = response.json()
        status_code = response['status']

    if status_code != 'success':
        print('Not Found.')
        jprint(response)

    ims = response['response'][0]['metadata']['xvalues']
    maf = response['response'][0]['data'][0]['yvalues']

    hc = {'ims': ims, 'maf': maf}
    return hc


def getUHSfromUSGS(parameters, Tr):
    # INPUT
    #   parameters = dict with info of the ste
    #                     {
    #                           "edition":  ["E2008", "E2014"]
    #                           "region":   ["COUS", "WUS", "CEUS"]
    #                           "latitude": XXX,
    #                           "longitude": XXXX,
    #                           "vs30": [180, 259, 360, 537, 760, 1150, 2000] m/s
    #                     }
    #  Tr = reutrn period desired for UHS
    #
    # OUTPUTS
    #   spectra = dict list with values of the IM of interest
    #          Ts  = list o periods
    #          ims = Sa(Ts) value for the return period selected
    #

    imt = ['PGA', 'SA0P1', 'SA0P2', 'SA0P3', 'SA0P5', 'SA0P75', 'SA1P0', 'SA2P0', 'SA3P0', 'SA4P0', 'SA5P0']
    T = np.array([0, 0.1, 0.2, 0.3, 0.5, 0.75, 1, 2, 3, 4, 5])

    # Read hazard curves
    hc = []
    for im_i in imt:
        parameters['imt'] = im_i
        hc.append(getHCfromUSGS(parameters))

    # Interpolate return period (in log space)
    maf_tgt = 1 / Tr
    #     print(maf_tgt)
    ims = np.zeros(len(imt))

    for i, im_i in enumerate(imt):
        maf_list = hc[i]['maf']

        if maf_tgt <= max(maf_list) and maf_tgt >= min(maf_list):
            #             print(maf_list)
            #             print(hc[i]['ims'])
            f = sp.interpolate.interp1d(np.log(maf_list), np.log(hc[i]['ims']))
            ims[i] = np.exp(f(np.log(maf_tgt)))
        else:
            print('Target return period out of range of hazard curve: T = ' + str(Ts[i]))
            return

    spectra = {'T': T, 'ims': ims.flatten()}
    return spectra


def DE_NEHRP(Sds, Sd1, TL, T):
    # DE_NEHRP computes the Design Spectrum according to NEHRP (2009) in g

    # By. Francisco Galvis & Anne Hulsey, Stanford University (October, 05, 2017)

    # Input variables
    #  Sds: Spectral ordinate at 0.2s
    #  Sd1: Spectral ordinate at 1.0s
    #  TL: Long period limit response
    #  T: vector o scalar of desired periods

    # Output variables
    #  PS_t: Spectral ordinates of design spectrum

    To = 0.2 * Sd1 / Sds;
    Ts = Sd1 / Sds;
    n = len(T);
    PS_t = np.zeros(n);

    for i in range(n):
        if T[i] <= To:
            PS_t[i] = Sds * (0.4 + 0.6 * T[i] / To);
        elif T[i] <= Ts:
            PS_t[i] = Sds;
        elif T[i] <= TL:
            PS_t[i] = Sd1 / T[i];
        else:
            PS_t[i] = Sd1 * TL / (T[i] ^ 2);

    return PS_t


def getASCE7Spectra(parameters, dT, Tmax, spectraType):
    # INPUT
    #   parameters = dict with info of the ste
    #                     {
    #                           "edition":  ["E2008", "E2014"]
    #                           "region":   ["COUS", "WUS", "CEUS"]
    #                           "latitude": XXX,
    #                           "longitude": XXXX,
    #                           "riskCategory": ['I', 'II', 'III', 'IV'],
    #                           "siteClass: ['A', 'B', 'C', 'D', 'E'],
    #                           "title":''
    #                     }
    #  Tr = return period desired for UHS
    #
    # OUTPUTS
    #   spectra = dict list with values of the IM of interest
    #          Ts  = list o periods
    #          ims = Sa(Ts) value for the return period selected
    #

    url = "https://earthquake.usgs.gov/ws/designmaps/asce7-16.json"

    response = requests.get(url, params=parameters)
    response = response.json()

    site_parameters = {}

    #     parameters.PGAuh = response['response']['data']['pgauh']
    #     parameters.PGAd = response['response']['data']['pgad']
    site_parameters['PGA'] = response['response']['data']['pga']
    site_parameters['Fpga'] = response['response']['data']['fpga']

    site_parameters['Ssuh'] = response['response']['data']['ssuh']
    site_parameters['Ssd'] = response['response']['data']['ssd']
    site_parameters['Ss'] = response['response']['data']['ss']
    site_parameters['Fa'] = response['response']['data']['fa']

    site_parameters['S1uh'] = response['response']['data']['s1uh']
    site_parameters['S1d'] = response['response']['data']['s1d']
    site_parameters['S1'] = response['response']['data']['s1']
    site_parameters['Fv'] = response['response']['data']['fv']

    if site_parameters['Fv'] is None:
        site_parameters['Fv'] = 1.7
        print('WARNING: Requiere site response study, using Fv = 1.7')

    site_parameters['SeismicDesignCategory'] = response['response']['data']['sdc']
    site_parameters['TL'] = 8  # response['response']['data']['t_sub_l']

    # Compute spectra
    T = np.arange(round(Tmax / dT)) * dT
    if spectraType == 'DBE':
        Sds = 2 / 3 * site_parameters['Ss'] * site_parameters['Fa']
        Sd1 = 2 / 3 * site_parameters['S1'] * site_parameters['Fv']
        Sa = DE_NEHRP(Sds, Sd1, site_parameters['TL'], T)
    else:
        Sds = site_parameters['Ss'] * site_parameters['Fa']
        Sd1 = site_parameters['S1'] * site_parameters['Fv']
        Sa = DE_NEHRP(Sds, Sd1, site_parameters['TL'], T)

    spectra = {'T': T, 'ims': Sa}

    return spectra, site_parameters