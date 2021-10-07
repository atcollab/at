import matplotlib.pyplot as plt
from at.physics import Acceptance6D
import numpy as np
from datetime import datetime
import os

def plot(Acc6D, ax=None, file_name_save=None):
    """
    plots a figure of the acceptance for one of the defined modes

    :param Acc6D: instance of the class Acceptance6D
    :param ax: figure axes
    :param file_name_save: if given, save figure to file
    :return figure axes
    """

    if Acc6D.mode:
        print('auto plot {m}'.format(m=Acc6D.mode))

        h, v, sel = Acc6D.select_test_points_based_on_mode()
        ax = plot_base(Acc6D, h, v, sel, Acc6D.mode.split('-', 2), ax=ax, file_name_save=file_name_save)
    else:
        print('mode is None')

    return ax


def plot_base(Acc6D, h, v, sel=None, pl=('x', 'y'), ax=None, file_name_save=None, font_size=10):
    """
    plots results of acceptance scan in a given 2D plane

    :param Acc6D: instance of Acceptance6D class (in pyat/at/physics/dynamic_aperture.py)
    :param h: list of x coordinates of test points
    :param v: list of y coordinates of test points
    :param sel: boolean array to mark if the specified coordinate is lost or not
    :param pl: 2 element tuple of planes default ('x', 'y')
    :param ax: figure axes
    :param file_name_save: if given, save figure to file
    :return figure axes
    """
    if ax is None:
        fig, ax = plt.subplots()

    cols = []
    for s in Acc6D.survived:
        if s:
            cols.append('deepskyblue')
        else:
            cols.append('gainsboro')

    if (Acc6D.mode is None) or (Acc6D.mode == '6D'):
        if type(pl[0]) != str:
            pl_h = Acc6D.planes[pl[0]]
        else:
            pl_h = pl[0]

        if type(pl[1]) != str:
            pl_v = Acc6D.planes[pl[1]]
        else:
            pl_v = pl[1]
    else:
        # get planes from mode
        pl_h, pl_v = Acc6D.mode.split('-', 2)

    if not sel:
        sel = range(len(h))

    num_sel = [int(s) for s in sel]
    survsel = [int(s) for s in sel if s]
    if Acc6D.verbose:
        print('{p1}-{p2} survived in {n} cases'.format(n=len(survsel), p1=pl_h, p2=pl_v))

    # apply scale factors for plot
    hs = np.array([h_ * Acc6D.dict_units[pl_h][0] for h_ in h])
    vs = np.array([v_ * Acc6D.dict_units[pl_v][0] for v_ in v])

    ax.scatter(hs, vs, s=20, c=cols, label='tested', facecolors='none')
    ax.plot(hs[sel], vs[sel], 'x', color='royalblue', markersize=3, label='survived')
    ax.tick_params(axis='both', which='major', labelsize=font_size)
    try:
        cs = ax.tricontour(hs, vs, sel, linewidths=2)
        for c in cs.collections:
            c.set_edgecolor("face")

        if len(cs.allsegs) > 3:
            dat0 = cs.allsegs[-2][0]
            ax.plot(dat0[:, 0], dat0[:, 1], ':', label='limit')
        else:
            if Acc6D.verbose:
                print('DA limit could not be computed (probably no closed contour)')
    except Exception:
        print(' limits not computed ')

    ax.set_xlabel(pl_h + ' [' + Acc6D.dict_units[pl_h][1] + ']', fontsize=font_size)
    ax.set_ylabel(pl_v + ' [' + Acc6D.dict_units[pl_v][1] + ']', fontsize=font_size)
    ax.set_xlim([r * Acc6D.dict_units[pl_h][0] for r in Acc6D.dict_def_range[pl_h]])
    ax.set_ylim([r * Acc6D.dict_units[pl_v][0] for r in Acc6D.dict_def_range[pl_v]])

    ax.set_title('{m} for {t} turns\n at {ll} \n dp/p= {dpp}%, rad: {rad}'.format(
        m=Acc6D.mode, t=Acc6D.number_of_turns, ll=Acc6D.ring[0].FamName, dpp=Acc6D.dp * 100,
        rad=Acc6D.ring.radiation), fontsize=font_size)

    ax.legend(prop={'size': font_size})
    plt.tight_layout()

    if file_name_save:
        path, file = os.path.split(file_name_save)
        fns = file.split('.',2)

        if len(fns)>1:
            plt.savefig(path + '/' + fns[0] + '_' + Acc6D.mode.replace('-', '_') + '.' + fns[1], dpi=600)
        else:
            plt.savefig(path + '/' + fns[0] + '_' + Acc6D.mode.replace('-', '_'), dpi=600)

    return ax


def plot6d(Acc6D, axs_top=None, axs_bottom=None, file_name_save='./test.png'):
    """
    plot DA

    :param Acc6D: instance of Acceptance6D class (in pyat/at/physics/dynamic_aperture.py)

    :return:
    """

    # get 6D space center
    min_pl = {'x': 0.0,
              'xp': 0.0,
              'y': 0.0,
              'yp': 0.0,
              'delta': 0.0,
              'ct': 0.0,
              }

    for p in Acc6D.planes:
        min_pl[p]=min(abs(Acc6D.coordinates[p]))

    if not axs_top:
        fig, (axs_top, axs_bottom) = plt.subplots(2, 3)

    cols = []
    for s in Acc6D.survived:
        if s:
            cols.append('royalblue')
        else:
            cols.append('white')

    # plot X - Y
    axnum = 0
    h = Acc6D.coordinates['x']
    v = Acc6D.coordinates['y']

    sel = [a and b == min_pl['delta'] and c == min_pl['xp'] and d == min_pl['yp'] and e == min_pl['ct']
           for a, b, c, d, e in zip(
        Acc6D.survived,
        Acc6D.coordinates['delta'],
        Acc6D.coordinates['xp'],
        Acc6D.coordinates['yp'],
        Acc6D.coordinates['ct'])]

    axs_top[axnum] = plot_base(Acc6D, h, v, sel, ['x', 'y'], ax=axs_top[axnum])
    axs_top[axnum].get_legend().remove()

    # plot xp yp
    axnum = 1
    h = Acc6D.coordinates['xp']
    v = Acc6D.coordinates['yp']
    sel = [a and b == min_pl['delta'] and c == min_pl['x'] and d == min_pl['y'] and e == min_pl['ct']
           for a, b, c, d, e in zip(
        Acc6D.survived,
        Acc6D.coordinates['delta'],
        Acc6D.coordinates['x'],
        Acc6D.coordinates['y'],
        Acc6D.coordinates['ct'])]

    axs_top[axnum] = plot_base(Acc6D, h, v, sel, ['xp', 'yp'], ax=axs_top[axnum])
    axs_top[axnum].set_title('')
    axs_top[axnum].get_legend().remove()

    # plot ct delta
    axnum = 2
    h = Acc6D.coordinates['ct']
    v = Acc6D.coordinates['delta']
    sel = [a and b == min_pl['xp'] and c == min_pl['x'] and d == min_pl['y'] and e == min_pl['yp']
           for a, b, c, d, e in zip(
        Acc6D.survived,
        Acc6D.coordinates['xp'],
        Acc6D.coordinates['x'],
        Acc6D.coordinates['y'],
        Acc6D.coordinates['yp'])]

    axs_top[axnum] = plot_base(Acc6D, h, v, sel, ['ct', 'delta'], ax=axs_top[axnum])
    axs_top[axnum].set_title('')
    axs_top[axnum].get_legend().remove()

    # plot x xp
    axnum = 0
    h = Acc6D.coordinates['x']
    v = Acc6D.coordinates['xp']
    sel = [a and b == min_pl['delta'] and c == min_pl['ct'] and d == min_pl['y'] and e == min_pl['yp']
           for a, b, c, d, e in zip(
        Acc6D.survived,
        Acc6D.coordinates['delta'],
        Acc6D.coordinates['ct'],
        Acc6D.coordinates['y'],
        Acc6D.coordinates['yp'])]

    axs_bottom[axnum] = plot_base(Acc6D, h, v, sel, ['x', 'xp'], ax=axs_bottom[axnum])
    axs_bottom[axnum].set_title('')
    axs_bottom[axnum].get_legend().remove()

    # plot y yp
    axnum = 1
    h = Acc6D.coordinates['y']
    v = Acc6D.coordinates['yp']
    sel = [a and b == min_pl['delta'] and c == min_pl['ct'] and d == min_pl['x'] and e == min_pl['xp']
           for a, b, c, d, e in zip(
        Acc6D.survived,
        Acc6D.coordinates['delta'],
        Acc6D.coordinates['ct'],
        Acc6D.coordinates['x'],
        Acc6D.coordinates['xp'])]

    axs_bottom[axnum] = plot_base(Acc6D,h, v, sel, ['y', 'yp'], ax=axs_bottom[axnum])
    axs_bottom[axnum].set_title('')
    axs_bottom[axnum].get_legend().remove()

    # plot delta x
    axnum = 2
    h = Acc6D.coordinates['delta']
    v = Acc6D.coordinates['x']
    sel = [a and b == min_pl['y'] and c == min_pl['ct'] and d == min_pl['yp'] and e == min_pl['xp']
           for a, b, c, d, e in zip(
        Acc6D.survived,
        Acc6D.coordinates['y'],
        Acc6D.coordinates['ct'],
        Acc6D.coordinates['yp'],
        Acc6D.coordinates['xp'])]

    axs_bottom[axnum] = plot_base(Acc6D, h, v, sel, ['delta', 'x'], ax=axs_bottom[axnum])
    axs_bottom[axnum].set_title('')
    axs_bottom[axnum].get_legend().remove()

    plt.tight_layout()

    tit = axs_top[0].get_title().replace('\n', '')

    now = datetime.now()
    t_string = now.strftime("%d/%m/%Y %H:%M:%S")

    axs_top[0].set_title(t_string + '\n' + tit + '\n', fontdict={'fontsize': 12,
                                                   'verticalalignment': 'center',
                                                   'horizontalalignment': 'left'})

    if file_name_save:
        plt.savefig(file_name_save + '_6D', dpi=600)

    return axs_top, axs_bottom

def get_border(h, v, sel):
    """
    get border from output of Acc6D.compute()


    """

    h_s = []
    v_s = []

    fig, ax = plt.subplots()
    cs = ax.tricontour(h, v, sel, linewidths=2)

    if len(cs.allsegs) > 3:
        dat0 = cs.allsegs[-2][0]
        h_s = dat0[:, 0]
        v_s = dat0[:, 1]
    else:
        print('DA limit could not be computed (probably no closed contour)')

    plt.close(fig)

    return h_s, v_s


def plot_off_energy_dynamic_aperture(max_neg_x, deltaps, da, inject_from_inside = True, file_name_save = ''):
    """
    use output of dynamic_aperture.off_energy_dynamic_aperture to produce a figure
    """
    # make figure
    fig, ax = plt.subplots()
    fig.set_size_inches(8.0, 5.0)
    ll, bb, ww, hh = ax.get_position().bounds
    ax.set_position([ll+0.05, bb + 0.05, ww-0.05, hh - 0.05])
    ax.plot([dp*100 for dp in deltaps], [_mx*1e3 for _mx in max_neg_x])
    ax.set_xlabel('dpp [%]', fontsize=16)
    if inject_from_inside:
        ax.set_ylabel('max(x<0) [mm]', fontsize=16)
        ax.set_ylim([1.1*min(max_neg_x)*1e3, 0.0])
    else:
        ax.set_ylabel('max(x>0) [mm]', fontsize=16)
        ax.set_ylim([0.0, 1.1 * max(max_neg_x) * 1e3])
    plt.rcParams['font.size'] = '16'
    # Set tick font size
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(16)

    ax.set_title('{t} turns at {ll}'.format(t=da.number_of_turns, ll=da.ring[0].FamName))

    plt.grid()

    if file_name_save:
        plt.savefig(file_name_save, dpi=600)

    pass

def plot_momentum_acceptance(mom_acc, s, da=None, file_name_save=''):
    """
    use output of dynamic_aperture.momentum_acceptance to produce a figure

    :param mom_acc: 2xN array of maximum positive and negative momentum acceptance at N locations
    :param s: 1xN array, longitudinal location corresponding to the data in mom_acc
    :param da: instance of Acceptance6D class (in pyat/at/physics/dynamic_aperture.py)
    :param file_name_save: string, a file name to save the figure produced in png format

    """
    # make figure
    fig, ax = plt.subplots()
    fig.set_size_inches(8.0, 5.0)
    ll, bb, ww, hh = ax.get_position().bounds
    ax.set_position([ll+0.05, bb+0.05, ww-0.05, hh - 0.05])
    ax.plot(s, [dp * 100 for dp in mom_acc[0]])
    ax.plot(s, [dp * 100 for dp in mom_acc[1]])
    plt.rcParams['font.size'] = '16'
    # Set tick font size
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(16)
    ax.set_xlabel('s [m]', fontsize=16)
    ax.set_ylabel('delta p/p [%]', fontsize=16)

    if file_name_save:
        plt.savefig(file_name_save, dpi=600)

    pass

def plot_dynamic_aperture(h, v, da, search, file_name_save=''):
    """
    use output of dynamic_aperture.dynamic_aperture to produce a figure

    :param h: 1xN array of horizontal coordinates  (can be any coordinate, according to info stored in da input)
    :param v: 1xN array of vertical coordinates  (can be any coordinate, according to info stored in da input)
    :param da: instance of Acceptance6D class used for dynamic apertutre computation
               (in pyat/at/physics/dynamic_aperture.py)
    :param search: Boolean, output of dynamic_aperture.dynamic_aperture().
    :param file_name_save: string, a file name to save the figure produced in png format
    """

    if search:
        # fig, ax = plt.subplots()
        ax = plot(da)
        fig = plt.gcf()
        fig.set_size_inches(8.0, 5.0)
        ll, bb, ww, hh = ax.get_position().bounds
        ax.set_position([ll + 0.05, bb + 0.05, ww - 0.05, hh - 0.1])
        ax.plot([_h * da.dict_units['x'][0] for _h in h],
                [_v * da.dict_units['y'][0] for _v in v],
                color='darkgray', label='DA limit')
        ax.set_xlabel('x [' + da.dict_units['x'][1] + ']', fontsize=16)
        ax.set_ylabel('y [' + da.dict_units['y'][1] + ']', fontsize=16)
        plt.rcParams['font.size'] = '16'
        # Set tick font size
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(16)
        ax.set_title('{m} for {t} turns\n at {ll}, dp/p= {dpp}%, rad: {rad}\n'.format(
            m=da.mode, t=da.number_of_turns, ll=da.ring[0].FamName, dpp=da.dp * 100, rad=da.ring.radiation))
        ax.legend()
        if file_name_save:
            plt.savefig(file_name_save, dpi=600)
    else:
        ax = plot(da, file_name_save=file_name_save)

    pass