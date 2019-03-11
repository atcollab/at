"""AT generic plotting function"""
from itertools import chain, repeat
SLICES = 400

try:
    # noinspection PyPackageRequirements
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError) as exc:
    print(exc)
    print('Plotting is disabled')
    # noinspection PyUnusedLocal
    def baseplot(ring, plot_function, *args, **kwargs):
        raise ImportError('Matplotlib is not available: plotting is disabled')
else:
    def baseplot(ring, plot_function, *args, **kwargs):
        """
        baseplot divides the region of interest of ring into small elements,
        calls the specified function to get the plot data and calls matplotlib
        functions to generate the plot.
        By default it creates a new figure for the plot, but if provided with
        axes objects it can be used as part of a GUI

        PARAMETERS
            ring            Lattice object
            plot_function   specific data generating function to be called

            All other positional parameters are sent to the plotting function

            plot_function is called as:

            title, left, right = plot_function(ring, refpts, *args, **kwargs)

            and should return 3 output:

            title   plot title or None
            left    tuple returning the data for the main (left) axis
            right   tuple returning the data for the secondary (right) axis
                    or None if no secondary axis is needed

                left[0]   xdata: (N,) array
                left[1]   ydata: list of (N,) or (N,M) arrays. Lines from a
                          (N, M) array share the same style and label
                left[2]   y-axis label
                left[3]   labels: (optional) list of strings as long as ydata

        KEYWORDS
            s_range=None    plot range, defaults to the full ring
            axes=None       axes for plotting as (primary_axes, secondary_axes)
                            Default: create new axes
            slices=400      Number of slices
            legend=False    Show a legend on the plot

            Other keywords are sent to the plotting function

        RETURN
            axis1, axes2    primary and secondary plot axes
        """
        def plot1(ax, x, y, yaxis_label='', labels=()):
            lines = []
            for y1, prop, label in zip(y, props, chain(labels, repeat(None))):
                ll = ax.plot(x, y1, **prop)
                if label is not None:
                    ll[0].set_label(label)
                lines += ll
            ax.set_ylabel(yaxis_label)
            return lines

        def labeled(line):
            return not line.properties()['label'].startswith('_')

        s_range = kwargs.pop('s_range', None)
        slices = kwargs.pop('slices', SLICES)
        axes = kwargs.pop('axes', None)
        legend = kwargs.pop('legend', False)
        try:
            cycle_props = plt.rcParams['axes.prop_cycle']
        except KeyError:
            cycle_props = plt.rcParams['axes.color_cycle']

        # slice the ring
        rg = ring.slice(s_range=s_range, slices=slices)

        # get the data for the plot
        title, left, right = plot_function(rg, rg.i_range, *args, **kwargs)

        # prepare the axes
        if axes is None:
            # Create new axes
            fig = plt.figure()
            ax1 = fig.add_subplot(111, title=title, xlim=rg.s_range, xlabel='s [m]')
            ax2 = None if right is None else ax1.twinx()
        else:
            # Use existing axes
            ax1, ax2 = axes

        props = iter(cycle_props())

        # left plot
        lines1 = plot1(ax1, *left)
        # right plot
        lines2 = [] if (right is None or ax2 is None) else plot1(ax2, *right)
        if legend:
            if ax2 is None:
                ax1.legend(handles=[l for l in lines1 if labeled(l)])
            elif ax1.get_shared_x_axes().joined(ax1, ax2):
                ax1.legend(handles=[l for l in lines1+lines2 if labeled(l)])
            else:
                ax1.legend(handles=[l for l in lines1 if labeled(l)])
                ax2.legend(handles=[l for l in lines2 if labeled(l)])
        plt.show()
        return ax1, ax2
