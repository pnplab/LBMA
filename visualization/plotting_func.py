import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math
from scipy import stats


# %%
#  Function to convert 1D vector to symmetric matrix

def vec_to_sym(vec):
    x = len(vec)
    n = int((math.sqrt(8 * x + 1) - 1) / 2)

    idx = np.triu_indices(n, k=0, m=n)
    matrix = np.zeros((n, n)).astype(int)
    matrix[idx] = vec
    matrix_sym = matrix + matrix.T - np.diag(np.diag(matrix))
    return matrix_sym


# %%

def plot_reg_analysis(data, contrast, X_type, Y_type, color_list, ax=None):
    max_axis = int(math.ceil(data.max().max() / 10.0)) * 10
    ax.set_xlim(-5, max_axis);
    ax.set_ylim(-5, max_axis);
    ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3");
    ax.set_aspect('equal');

    for i in contrast:
        X = data[i][X_type]
        Y = data[i][Y_type]
        sns.regplot(x=X, y=Y, scatter=True, ci=None, robust=True, color=color_list[i], x_jitter=0.2, y_jitter=0.2,
                    ax=ax);
    return


# %%# Plot graphs and calculate correlation

def plot_corr_analysis(data, X_contrast, X_type, Y_contrast, Y_type, significance_table=None, ax=None):
    X_reg = data[X_contrast][X_type]
    Y_reg = data[Y_contrast][Y_type]
    corr_pears, p_value = stats.pearsonr(X_reg, Y_reg)
    max_axis = int(math.ceil(data.max().max() / 10.0)) * 10

    if significance_table is None:
        style = None
        style_order = None
        markers = "o"
    else:
        X_significance = significance_table[X_contrast][X_type].replace([0, 1], ['X Non-significant', 'X Significant'])
        Y_significance = significance_table[Y_contrast][Y_type].replace([0, 1], ['Y Non-significant', 'Y Significant'])
        XY_significance = X_significance + '/' + Y_significance
        style = XY_significance
        style_order = ['X Significant/Y Significant', 'X Non-significant/Y Non-significant',
                       'X Non-significant/Y Significant', 'X Significant/Y Non-significant']
        markers = ["o", "P", "v", "^"]

    sns.scatterplot(x=X_reg, y=Y_reg, style=style,
                    style_order=style_order,
                    markers=markers, color='dimgray', s=100, ax=ax);
    ax.set_xlim(-5, max_axis);
    ax.set_ylim(-5, max_axis);
    ax.set_aspect('equal');
    ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3");
    sns.regplot(x=X_reg, y=Y_reg, scatter=None, color='dimgray', ci=None, ax=ax);

    if not (X_contrast == Y_contrast):
        ax.set_xlabel(X_contrast);
        ax.set_ylabel(Y_contrast);
    ax.text(0.10, 62, 'corr_coeff = @ \n p_value = %'.replace('@', str(round(corr_pears, 3))).replace('%', str(
        round(p_value, 3))), fontsize=12);
    return


# %% Plot pie chart with subgroups

def plot_double_pie(group_size, subgroup_size, group_names, subgroup_names, outer_colors, inner_colors, ax):
    # First Ring (outside)
    ax.axis('equal');
    mypie, _ = ax.pie(group_size, radius=1.3, labels=group_names, colors=outer_colors, textprops={'fontsize': 12},
                      counterclock=False);
    plt.setp(mypie, width=0.3, edgecolor='white');

    # Second Ring (Inside)
    mypie2, _ = ax.pie(subgroup_size, radius=1.3 - 0.3, labels=subgroup_names, labeldistance=0.45, colors=inner_colors,
                       textprops={'fontsize': 12}, counterclock=False);
    plt.setp(mypie2, width=0.4, edgecolor='white');
    plt.margins(0, 0);
    return


# %%

def plot_scatter3D(x, y, z, N, colors='black', markers='o', ax3d=None):
    ax3d.scatter(x, y, z, marker=markers, color=colors, s=100);

    ax3d.scatter(x, z, marker=markers, color=colors, zdir='y', zs=80, alpha=0.3)
    ax3d.scatter(y, z, marker=markers, color=colors, zdir='x', zs=0, alpha=0.3)
    ax3d.scatter(x, y, marker=markers, color=colors, zdir='z', zs=0, alpha=0.3)

    for i in range(N):
        ax3d.plot(np.array([0, x[i]]), np.array([y[i], y[i]]),
                  np.array([z[i], z[i]]), 'k--', alpha=0.1);
        ax3d.plot(np.array([x[i], x[i]]), np.array([70, y[i]]),
                  np.array([z[i], z[i]]), 'k--', alpha=0.1);
        ax3d.plot(np.array([x[i], x[i]]), np.array([y[i], y[i]]),
                  np.array([0, z[i]]), 'k--', alpha=0.1);

    ax3d.view_init(azim=-42, elev=20);
    return


# %%

def colorbar(palette, n_colors, color_min, color_max, ax):
    col_x = [0] * len(palette)  # Fixed x coordinate for the bars
    bar_y = np.linspace(color_min, color_max, n_colors)  # y coordinates for each of the n_colors bars

    bar_height = bar_y[1] - bar_y[0]
    ax.barh(
        y=bar_y,
        width=[5] * len(palette),  # Make bars 5 units wide
        left=col_x,  # Make bars start at 0
        height=bar_height,
        color=palette,
        linewidth=0
    )
    ax.set_xlim(1, 2)  # Bars are going from 0 to 5, so lets crop the plot somewhere in the middle
    ax.grid(False)  # Hide grid
    ax.set_facecolor('white')  # Make background white
    ax.set_xticks([])  # Remove horizontal ticks
    ax.set_yticks(np.linspace(min(bar_y), max(bar_y), 3))  # Show vertical ticks for min, middle and max
    ax.yaxis.tick_right()  # Show vertical ticks on the right
    return


# %%

# Heatmap with sized markers

def heatmap(x, y, ax=None, linewidths=None, linestyle=None, palette_list=None, **kwargs):
    if 'color' in kwargs:
        color = kwargs['color']
    else:
        color = [1] * len(x)

    if 'palette' in kwargs:
        palette = kwargs['palette']
        n_colors = len(palette)

    else:
        n_colors = 256  # Use 256 colors for the diverging color palette
        palette = sns.color_palette("Blues", n_colors)

    if 'color_range' in kwargs:
        color_min, color_max = kwargs['color_range']
    else:
        color_min, color_max = min(color), max(color)  # Range of values that will be mapped to the palette,
        # i.e. min and max possible correlation

    def value_to_color(val):
        if color_min == color_max:
            return palette[-1]
        else:
            val_position = float((val - color_min)) / (
                    color_max - color_min)  # position of value in the input range, relative to the length of the input range
            val_position = min(max(val_position, 0), 1)  # bound the position betwen 0 and 1
            ind = int(val_position * (n_colors - 1))  # target index in the color palette
            return palette[ind]

    if 'size' in kwargs:
        size = kwargs['size']
    else:
        size = [1] * len(x)

    if 'size_range' in kwargs:
        size_min, size_max = kwargs['size_range'][0], kwargs['size_range'][1]
    else:
        size_min, size_max = min(size), max(size)

    size_scale = kwargs.get('size_scale', 500)

    def value_to_size(val):
        if size_min == size_max:
            return 1 * size_scale
        else:
            val_position = (val - size_min) * 0.99 / (size_max - size_min) + 0.01
            val_position = min(max(val_position, 0), 1)  # bound the position betwen 0 and 1
            return val_position * size_scale

    if 'x_order' in kwargs:
        x_names = [t for t in kwargs['x_order']]
    else:
        x_names = [t for t in sorted(set([v for v in x]))]
    x_to_num = {p[1]: p[0] for p in enumerate(x_names)}

    if 'y_order' in kwargs:
        y_names = [t for t in kwargs['y_order']]
    else:
        y_names = [t for t in sorted(set([v for v in y]))]
    y_to_num = {p[1]: p[0] for p in enumerate(y_names)}

    bool = 0
    if ax is None:
        plot_grid = plt.GridSpec(1, 15, hspace=0.2, wspace=0.1)  # Setup a 1x10 grid
        ax = plt.subplot(plot_grid[:, :-1])  # Use the left 14/15ths of the grid for the main plot
        bool = True

    marker = kwargs.get('marker', 's')

    kwargs_pass_on = {k: v for k, v in kwargs.items() if k not in [
        'color', 'palette', 'color_range', 'size', 'size_range', 'size_scale', 'marker', 'x_order', 'y_order'
    ]}

    if isinstance(marker, pd.Series):
        for mark in marker.unique():
            ax.scatter(
                x=[x_to_num[v] for v in x.loc[(marker == mark)]],
                y=[y_to_num[v] for v in y.loc[(marker == mark)]],
                marker=mark,
                s=[value_to_size(v) for v in size.loc[(marker == mark)]],
                c=[value_to_color(v) for v in color.loc[(marker == mark)]],
                linewidths=linewidths.loc[(marker == mark)],
                linestyle=linestyle.loc[(marker == mark)],
                edgecolors='black',
                **kwargs_pass_on
            )

    elif palette_list is not None:
        n_colors = 256
        for colorCode in palette_list.unique():
            palette = sns.light_palette(colorCode, n_colors)
            ax.scatter(
                x=[x_to_num[v] for v in x.loc[(palette_list == colorCode)]],
                y=[y_to_num[v] for v in y.loc[(palette_list == colorCode)]],
                s=[value_to_size(v) for v in size.loc[(palette_list == colorCode)]],
                c=[value_to_color(v) for v in color.loc[(palette_list == colorCode)]],
                marker=marker,
                linewidths=linewidths.loc[(palette_list == colorCode)],
                linestyle=linestyle.loc[(palette_list == colorCode)],
                edgecolors='black',
                **kwargs_pass_on
            )
    else:
        ax.scatter(
            x=[x_to_num[v] for v in x],
            y=[y_to_num[v] for v in y],
            marker=marker,
            s=[value_to_size(v) for v in size],
            c=[value_to_color(v) for v in color],
            linewidths=linewidths,
            linestyle=linestyle,
            edgecolors='black',
            **kwargs_pass_on
        )
    ax.set_xticks([v for k, v in x_to_num.items()])
    ax.set_xticklabels([k for k in x_to_num], rotation=45, horizontalalignment='right')
    ax.set_yticks([v for k, v in y_to_num.items()])
    ax.set_yticklabels([k for k in y_to_num])

    ax.grid(False, 'major')
    ax.grid(True, 'minor')
    ax.set_xticks([t + 0.5 for t in ax.get_xticks()], minor=True)
    ax.set_yticks([t + 0.5 for t in ax.get_yticks()], minor=True)

    ax.set_xlim([-0.5, max([v for v in x_to_num.values()]) + 0.5])
    ax.set_ylim([-0.5, max([v for v in y_to_num.values()]) + 0.5])
    ax.set_facecolor('#F1F1F1')

    if bool:
        # Add color legend on the right side of the plot
        if color_min < color_max:
            ax = plt.subplot(plot_grid[:, -1])  # Use the rightmost column of the plot
            colorbar(palette, n_colors, color_min, color_max, ax=ax)
    return

# %%
