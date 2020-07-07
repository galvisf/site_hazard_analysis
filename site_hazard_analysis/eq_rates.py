from .base import *


def m_and_r_rates_from_OpenSHA(eq, all_seismicity, r_list, m_list, plot_flag, save_flag):

    distances = all_seismicity['DistanceRup']
    magnitudes = all_seismicity['Magnitude']
    rates = all_seismicity['MeanAnnualRate']

    rate_r_m = np.zeros([len(m_list), len(r_list)])
    n_ruptures_r_m = np.zeros_like(rate_r_m)
    for i_m in range(len(m_list ) -1):
        m_min = m_list[i_m]
        m_max = m_list[ i_m +1]

        magnitude_idx = (magnitudes >= m_min) & (magnitudes < m_max)
        for i_r in range(len(r_list ) -1):
            r_min = r_list[i_r]
            r_max = r_list[ i_r +1]

            distance_idx = (distances >= r_min) & (distances < r_max)
            combined_idx = magnitude_idx & distance_idx

            n_ruptures_r_m[i_m, i_r] = np.sum(combined_idx)
            rate_r_m[i_m, i_r] = np.sum(rates[combined_idx])

    rates = pd.DataFrame(rate_r_m, index=m_list, columns=r_list)

    if save_flag:
        rates.to_excel( eq +'_rates.xlsx')

    if plot_flag:
        fig, ax = plt.subplots(1 ,1 ,figsize=(12 ,5))
        _ = sns.heatmap(n_ruptures_r_m, linewidth=0.5)
        _ = ax.set_xlabel('Site to Source Distance, R')
        _ = ax.set_xticks(np.arange(len(r_list)))
        _ = ax.set_xticklabels(r_list, rotation=90)
        _ = ax.set_ylabel('Magnitude, M')
        _ = ax.set_yticklabels(m_list, rotation=0)
        _ = ax.set_title('OpenSHA Rupture List, Magnitude and Distance')

        if save_flag:
            _ = plt.tight_layout()
            filename = eq + '_Magnitude and Distance Counts.png'
            plt.savefig(filename)
        _ = plt.show()


        fig, ax = plt.subplots(1 ,1 ,figsize=(12 ,5))
        _ = sns.heatmap(rate_r_m, linewidth=0.5)
        _ = ax.set_xlabel('Site to Source Distance, R')
        _ = ax.set_xticks(np.arange(len(r_list)))
        _ = ax.set_xticklabels(r_list, rotation=90)
        _ = ax.set_ylabel('Magnitude, M')
        _ = ax.set_yticklabels(m_list, rotation=0)
        _ = ax.set_title('Annualized Rates for Magnitude and Distance')

        if save_flag:
            _ = plt.tight_layout()
            filename = eq + '_Magnitude and Distance Rates.png'
            plt.savefig(filename)
        _ = plt.show()

    return rates


def get_site_fault_geometry(sites, faults, fault_list, plot):
    site_points = sites['x,y'].values

    site_names = list()
    fault_names = list()
    segment_idx = list()
    d_values = list()
    p_values = list()
    l1_values = list()
    l2_values = list()

    for i_site in range(len(sites)):
        if plot:
            fig, ax = plt.subplots(1, 1)
            _ = plt.scatter(site_points[i_site][0], site_points[i_site][1], color='k')

        for fault in fault_list:
            fault_points = faults.loc[fault, 'x,y'].values
            n_segments = len(fault_points) - 1

            for i_segment in range(n_segments):
                a = fault_points[i_segment]
                b = fault_points[i_segment + 1]
                c = site_points[i_site]
                [d, p, l1, l2] = perpendicular_distance(a, b, c)

                site_names.append(sites.index[i_site])
                fault_names.append(fault)
                segment_idx.append(i_segment)
                p_values.append(p)

                d_values.append(d)
                l1_values.append(l1)
                l2_values.append(l2)

                if plot:
                    x = [a[0], b[0]]
                    y = [a[1], b[1]]
                    _ = plt.plot(x, y, color='C' + str(i_segment))
                    _ = plt.scatter(p[0], p[1], color='C' + str(i_segment))
                    _ = plt.text((x[0] + p[0]) / 2, (y[0] + p[1]) / 2, '{0:.0f}'.format(l1))
                    _ = plt.text((x[1] + p[0]) / 2, (y[1] + p[1]) / 2, '{0:.0f}'.format(l2))

            if plot:
                _ = plt.axis('equal')
                _ = plt.show()

    site_fault_geometry = pd.DataFrame()
    site_fault_geometry['Site'] = site_names
    site_fault_geometry['Fault'] = fault_names
    site_fault_geometry['Segment'] = segment_idx
    site_fault_geometry['Closest Point on Fault Vector'] = p_values
    site_fault_geometry['d'] = d_values
    site_fault_geometry['l1'] = l1_values
    site_fault_geometry['l2'] = l2_values

    return site_fault_geometry


def perpendicular_distance(a, b, c):
    # project c onto line spanned by a,b but consider the end points
    # should the projection fall "outside" of the segment
    n, v = b - a, c - a

    # the projection q of c onto the infinite line defined by points a,b
    # can be parametrized as q = a + t*(b - a). In terms of dot-products,
    # the coefficient t is (c - a).(b - a)/( (b-a).(b-a) ). If we want
    # to restrict the "projected" point to belong to the finite segment
    # connecting points a and b, it's sufficient to "clip" it into
    # interval [0,1] - 0 corresponds to a, 1 corresponds to b.

    clip = False
    if clip:
        t = max(0, min(np.dot(v, n) / np.dot(n, n), 1))
    else:
        t = np.dot(v, n) / np.dot(n, n)

    # distance to the line
    d = np.linalg.norm(v - t * n)

    # nearest point on the line
    p = a + n * t

    n_unit = n / np.linalg.norm(n)
    l1 = n * t
    l1 = np.dot(l1, n_unit)
    l2 = n * (1 - t)
    l2 = np.dot(l2, n_unit)

    return d, p, l1, l2

def cut_line(line, distance):
    # Cuts a line in two at a distance from its starting point
    if distance <= 0.0 or distance >= line.length:
        return [LineString(line)]
    coords = list(line.coords)
    for i, p in enumerate(coords):
        pd = line.project(Point(p))
        if pd == distance:
            return [
                LineString(coords[:i+1]),
                LineString(coords[i:])]
        if pd > distance:
            cp = line.interpolate(distance)
            return [
                LineString(coords[:i] + [(cp.x, cp.y)]),
                LineString([(cp.x, cp.y)] + coords[i:])]


def m_to_length(m, fault_type):
    # calculate rupture length based on Kramar 1999 empirical relationships
    # current implementation assumes a deterministic value for rup_length
    if fault_type == 'strike-slip':
        rup_length = 10 ** (0.74 * m - 3.55)
        sigma = 0.23
    elif fault_type == 'reverse':
        rup_length = 10 ** (0.63 * m - 2.86)
        sigma = 0.20
    elif fault_type == 'normal':
        rup_length = 10 ** (0.50 * m - 2.01)
        sigma = 0.21
    else:
        'unknown fault type'

    return rup_length


def m_to_area(m, fault_type):
    # calculate rupture area based on Kramar 1999 empirical relationships
    # current implementation assumes a deterministic value for rup_area
    if fault_type == 'strike-slip':
        rup_area = 10 ** (0.90 * m - 3.42)
        sigma = 0.22
    elif fault_type == 'reverse':
        rup_area = 10 ** (0.98 * m - 3.99)
        sigma = 0.26
    elif fault_type == 'normal':
        rup_area = 10 ** (0.82 * m - 2.87)
        sigma = 0.22
    else:
        'unknown fault type'

    return rup_area


def prob_distance(m, fault_segment, fault_type, dr, plot_flag):

    rup_length = m_to_length(m, fault_type)

    # calculate the range of relevant distances
    d, l1, l2 = fault_segment[['d', 'l1', 'l2']].values[0]

    l_max = max(max(l1, l2) - rup_length, 0)
    r_max = np.sqrt(d ** 2 + l_max ** 2)

    if (l1 >= 0) & (l2 >= 0):
        r_min = d
    else:
        l_min = -min(l1, l2)
        r_min = np.sqrt(d ** 2 + l_min ** 2)

        r_max = max(r_min, r_max)

    # calculate the probability of each distance
    r_list = np.arange(0, r_max + 2 * dr, dr)

    # if the first r within range is only fraction*dr over r_min,
    #    increase r_min in order to skip to the next value of r
    fraction = 0.25
    idx = np.searchsorted(r_list, r_min)
    if r_min + fraction * dr >= r_list[idx]:
        r_min = r_min + dr
        # reset r_max if it is now less than r_min
        #    this keeps the r_min equivalent if the magnitude ruptures the whole fault
        if r_min > r_max:
            r_max = r_min
    # similar correction, for case where the site is offset from the fault
    if (l1 < 0) | (l2 < 0):
        l_min = -min(l1, l2) + (1 - fraction) * dr
        r_min = np.sqrt(l_min ** 2 + d ** 2)
        # reset r_max if it is now less than r_min
        #    this keeps the r_min equivalent if the magnitude ruptures the whole fault
        if r_min > r_max:
            r_max = r_min

    cdf_dist = list()
    for r in r_list:
        if r < r_min:
            cdf_dist.append(0)
        elif r >= r_max:
            cdf_dist.append(1)
        else:
            len_within_r = np.sqrt(r ** 2 - d ** 2)

            if (l1 < 0) | (l2 < 0):
                total_available = max(l1 + l2 - rup_length, 0)
                total_contribution = max(len_within_r + min(l1, l2), 0)

            else:
                l1_available = max(l1 - rup_length, 0)
                l1_contribution = min(len_within_r, l1_available)

                l2_available = max(l2 - rup_length, 0)
                l2_contribution = min(len_within_r, l2_available)

                overlap_available = min(rup_length, min(l1, l2))
                overlap_contribution = min(rup_length, overlap_available)

                total_available = l1_available + l2_available + overlap_available
                total_contribution = l1_contribution + l2_contribution + overlap_contribution

            cdf_dist.append(total_contribution / total_available)

    cdf_dist = np.array(cdf_dist)
    p_dist = cdf_dist[1:] - cdf_dist[:-1]
    p_dist = np.append(p_dist, 0)

    if plot_flag:
        fig, ax = plt.subplots(1, 1)
        color = 'C0'
        if len(r_list) > 1:
            _ = ax.plot(r_list, p_dist, color=color, marker='o')
        else:
            _ = ax.scatter(r_list, p_dist, color=color)
        _ = ax.set_xlim(left=0)
        _ = ax.set_xlabel('r')
        _ = ax.set_ylim(bottom=0)
        _ = ax.set_ylabel('PDF', color=color)

        ax2 = ax.twinx()
        color = 'C1'
        if len(r_list) > 1:
            _ = ax2.plot(r_list, cdf_dist, color=color)
        else:
            _ = ax2.scatter(r_list, cdf_dist, color=color)
        _ = ax2.set_ylabel('CDF', color=color)
        _ = ax2.set_ylim([0, 1])

        _ = plt.show()

    return r_list, p_dist


def prob_magnitude(m_max, m_min, b, dm, plot_flag):
    m_list = np.arange(m_min, m_max + 0.1, dm)

    cdf_mag = np.array([(1 - 10 ** (-b * (m - m_min))) / (1 - 10 ** (-b * (m_max - m_min))) for m in m_list])
    p_mag = cdf_mag[1:] - cdf_mag[:-1]
    p_mag = np.append(p_mag, 0)

    if plot_flag:
        fig, ax = plt.subplots(1, 1)
        _ = plt.bar(m_list, p_mag, width=dm, edgecolor='k')
        _ = plt.xlabel('Magnitude')
        _ = plt.ylabel('P(M = m)')
        _ = plt.show()

    return m_list, p_mag


def select_mainshock_rupture(site, ruptures, m_mainshock, dm, fault):
    m_min = m_mainshock - dm / 2
    m_max = m_mainshock + dm / 2

    magnitude_idx = (ruptures['Magnitude'] > m_min) & (ruptures['Magnitude'] <= m_max)
    fault_idx = ruptures['Name'].str.contains(fault)
    idx = magnitude_idx & fault_idx
    r_min = ruptures[idx]['DistanceRup'].min()
    distance_idx = (ruptures['DistanceRup'] >= r_min) & (ruptures['DistanceRup'] < r_min + 0.1)
    idx = idx & distance_idx
    idx = np.around(np.median(np.where(idx)[0])).astype('int')

    return idx


def rupture_segment_to_line(ruptures, idx):
    crs = ruptures.crs['init']
    if crs != 'epsg:4326':  # wgs84
        raise Warning('ruptures must be in WGS84 coordinates')

    points = list(ruptures['geometry'].iloc[idx].coords)
    points = [points[0], points[-1]]
    df_points = gpd.GeoDataFrame(geometry=[Point(xy) for xy in points])
    df_points['Label'] = 'Mainshock'
    df_points.set_index('Label', inplace=True)

    longitude = np.array([x for x, y in points])
    latitude = np.array([y for x, y in points])
    df_points['Latitude'] = latitude
    df_points['Longitude'] = longitude

    [x, y] = utm_conversion(latitude, longitude)
    x = x / 1000
    y = y / 1000
    df_points['x'] = x
    df_points['y'] = y
    df_points['x,y'] = [np.array([i, j]) for i, j in zip(x, y)]
    df_points = df_points[['Latitude', 'Longitude', 'x', 'y', 'x,y', 'geometry']]

    df_line = df_points.copy()
    df_line = df_line.groupby(['Label'])['geometry'].apply(lambda x: LineString(x.tolist()))
    df_line = gpd.GeoDataFrame(df_line, geometry='geometry')
    df_line['Name'] = ruptures['Name'][idx]
    df_line['Magnitude'] = ruptures['Magnitude'][idx]
    df_line['Distance'] = ruptures['DistanceRup'][idx]
    df_line['Mean Annual Frequency'] = ruptures['MeanAnnualRate'][idx]
    df_line = df_line[['Name', 'Magnitude', 'Distance', 'Mean Annual Frequency', 'geometry']]

    df_points.crs = crs
    df_line.crs = crs

    return df_points, df_line


def instantaneous_rate_of_aftershocks(t, m_max, aftershock_parameters):
    a = aftershock_parameters['a']
    b = aftershock_parameters['b']
    p = aftershock_parameters['p']
    c = aftershock_parameters['c']

    m_min = aftershock_parameters['m_min']

    mu = (10 ** (a + b * (m_max - m_min)) - 10 ** a) / ((t + c) ** p)

    return mu


def mean_number_of_aftershocks(t, T, m_max, aftershock_parameters):
    a = aftershock_parameters['a']
    b = aftershock_parameters['b']
    p = aftershock_parameters['p']
    c = aftershock_parameters['c']

    m_min = aftershock_parameters['m_min']

    mu = ((10 ** (a + b * (m_max - m_min)) - 10 ** a) / (p - 1)) * ((t + c) ** (1 - p) - (t + T + c) ** (1 - p))

    return mu


def mean_occurance_surface(m_max, aftershock_parameters, t_values, T_values, plot_flag):

    mean_occurance = np.zeros([len(T_values), len(t_values)])
    for i_T, T in enumerate(T_values):
        mean_occurance[i_T,:] = [mean_number_of_aftershocks(t=t, T=T, m_max=m_max, aftershock_parameters=aftershock_parameters) for t in t_values]
    mean_occurance = pd.DataFrame(mean_occurance, index=T_values, columns=t_values)

    if plot_flag:
        fig = plt.figure(figsize=(10,5))
        elev = 40
        azim = 40
        ax = fig.add_subplot(111, projection='3d', elev=elev, azim=azim)

        t_grid,T_grid = np.meshgrid(t_values, T_values)
        _ = ax.plot_surface(np.log(t_grid), np.log(T_grid), np.log(mean_occurance.values), alpha=0.5)

        zticks = np.array([0.0001, 0.001, 0.01, 0.1, 1, 10])
        _ = ax.set_zticks(np.log(zticks))
        _ = ax.set_zticklabels(zticks)
        _ = ax.set_zlabel('Mean Number of Aftershocks\n')

        _ = ax.set_xticks(np.log(t_values))
        _ = ax.set_xticklabels(t_values)
        _ = ax.set_xlim([min(np.log(t_values)), max(np.log(t_values))])
        _ = ax.set_xlabel('\nElapsed Time in Days, t')

        _ = ax.set_yticks(np.log(T_values))
        _ = ax.set_yticklabels(T_values)
        _ = ax.set_ylim([min(np.log(T_values)), max(np.log(T_values))])
        _ = ax.set_ylabel('\nTime Frame Considered, T')

        _ = ax.set_title('Mainshock Magnitude: ' + str(m_max) + 'Mw')
        _ = plt.tight_layout()
        _ = plt.show()

    return mean_occurance


def aftershock_r_and_m_distribution(site_fault_geometry, m_min, m_max, b, dm, dr, plot_flag):
    fault_type = 'strike-slip'

    [m_list, p_mag] = prob_magnitude(m_max, m_min, b, dm, plot_flag=False)

    for m, i_m in zip(m_list, range(len(m_list))):

        [r_list, p_dist] = prob_distance(m, site_fault_geometry, fault_type, dr, plot_flag=False)

        if i_m == 0:
            prob_r_given_m = np.zeros([len(m_list), len(r_list)])
            prob_r_given_m[i_m, :] = p_dist
            prob_r_m = np.zeros([len(m_list), len(r_list)])
            prob_r_m[i_m, :] = p_dist * p_mag[i_m]
            r_labels = r_list
        else:
            n_r = len(r_list)
            prob_r_given_m[i_m, :n_r] = p_dist
            prob_r_m[i_m, :n_r] = p_dist * p_mag[i_m]
        r_list = r_labels

    if plot_flag:
        fig, ax = plt.subplots(1, 1, figsize=(12, 5))
        _ = sns.heatmap(prob_r_given_m, linewidth=0.5)
        _ = ax.set_xlabel('Site to Source Distance, R')
        _ = ax.set_xticks(np.arange(len(r_labels)))
        _ = ax.set_xticklabels(r_labels, rotation=90)
        _ = ax.set_ylabel('Magnitude, M')
        _ = ax.set_yticklabels(m_list, rotation=0)
        _ = ax.set_title('PMF Distance | Magnitude')
        _ = plt.show()

        fig, ax = plt.subplots(1, 1, figsize=(12, 5))
        _ = sns.heatmap(prob_r_m, linewidth=0.5)
        _ = ax.set_xlabel('Site to Source Distance, R')
        _ = ax.set_xticks(np.arange(len(r_labels)))
        _ = ax.set_xticklabels(r_labels, rotation=90)
        _ = ax.set_ylabel('Magnitude, M')
        _ = ax.set_yticklabels(m_list, rotation=0)
        _ = ax.set_title('PMF Magnitude and Distance')
        _ = plt.show()

    return pd.DataFrame(prob_r_m, index=m_list, columns=r_list)