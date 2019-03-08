"""AT generic plotting function"""

SLICES = 400

try:
    # noinspection PyPackageRequirements
    import matplotlib.pyplot as plt
except ImportError:
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
                left[1]   ydata: (N,) or (N,M) array, for plotting M curves
                left[2]   y-axis label
                left[3]   (optional) list of M labels associated with each curve
                left[4]   (optional) dictionary of keyword arguments for 'plot'

        KEYWORDS
            s_range=None    plot range, defaults to the full ring
            axes=None       axes for plotting as (primary_axes, secondary_axes)
                            Default: create new axes
            slices=400      Number of slices

            Other keywords are sent to the plotting function

        RETURN
            axis1, axes2    primary and secondary plot axes
        """
        def plot1(ax, x, y, yaxis_label='', labels=None, options=None):
            if labels is None:
                labels = []
            if options is None:
                options = {}
            lines = ax.plot(x, y, **options)
            ax.set_ylabel(yaxis_label)
            for line, label in zip(lines, labels):
                line.set_label(label)

        s_range = kwargs.pop('s_range', None)
        slices = kwargs.pop('slices', SLICES)
        axes = kwargs.pop('axes', None)

        # slice the ring
        rg = ring.slice(s_range=s_range, slices=slices)

        # get the data for the plot
        title, left, right = plot_function(rg, rg.i_range, *args, **kwargs)

        # prepare the axes
        if axes is None:
            # Create new axes
            fig, ax1 = plt.subplots()
            ax1.set_xlim(*rg.s_range)
            ax1.set_xlabel('s [m]')
            ax1.set_title(title)
            ax2 = None if right is None else ax1.twinx()
        else:
            # Use existing axes
            ax1, ax2 = axes

        # left plot
        plot1(ax1, *left)
        # right plot
        if not ((right is None) or (ax2 is None)):
            plot1(ax2, *right)
        plt.show()
        return ax1, ax2
