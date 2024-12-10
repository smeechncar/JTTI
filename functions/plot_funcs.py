"""
plot_funcs.py

Created by: Jared A. Lee (jaredlee@ucar.edu)

This file contains functions for non-map plots.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, InsetPosition

def performance_diagram(df, obj_labels, colors, markers, fname, figsize=(8, 8), dpi=300,
                        xlabel='Success Ratio (1-FAR)',
                        ylabel='Probability of Detection (POD)', ticks=np.arange(0, 1.1, 0.1),
                        csi_cmap='Blues', csi_label='Critical Success Index',
                        title='Performance Diagram',
                        legend_params=None, bootstrap_sets=None, ci=(2.5, 97.5),
                        label_fontsize=14, title_fontsize=16, tick_fontsize=12):

    """
    Borrow much of the code from:
    https://github.com/djgagne/hagelslag/blob/master/hagelslag/evaluation/MetricPlotter.py

    A performance diagram is a variation on the ROC curve in which the Probability of False Detection on the
    x-axis has been replaced with the Success Ratio (1-False Alarm Ratio or Precision). The diagram also shows
    the Critical Success Index (CSI or Threat Score) as a series of curved contours, and the frequency bias as
    angled diagonal lines. Points along the 1:1 diagonal are unbiased, and better performing models should appear
    in the upper right corner. The performance diagram is particularly useful for displaying verification for
    severe weather warnings as it displays all three commonly used statistics (POD, FAR, and CSI) simultaneously
    on the same chart.
    """

    if legend_params is None:
        legend_params = dict(loc=4, fontsize=12, framealpha=1, frameon=True)
    plt.figure(figsize=figsize)
    grid_ticks = np.arange(0, 1.01, 0.01)
    sr_g, pod_g = np.meshgrid(grid_ticks, grid_ticks)
    bias = pod_g / sr_g
    csi = 1.0 / (1.0 / sr_g + 1.0 / pod_g - 1.0)
    csi_contour = plt.contourf(sr_g, pod_g, csi, np.arange(0.1, 1.1, 0.1), extend='max', cmap=csi_cmap)
    b_contour = plt.contour(sr_g, pod_g, bias, [0.5, 1, 1.5, 2, 4], colors='k', linestyles='dashed')
    plt.clabel(b_contour, fmt='%1.1f', manual=[(0.2, 0.9), (0.4, 0.9), (0.6, 0.9), (0.7, 0.7)])

    # In hagelslag, DJ created "DistributedROC objects". Running these DistributedROC objects through the
    # performance_curve() function returns a pandas dataframe with POD, FAR, and thresholds as columns.
    # bootstrap_sets allows for confidence intervals to be drawn in both directions for points on the diagram.
    # TODO: Eventually create a similar object that holds multiple dataframes, but start with just one for now.
    #if bootstrap_sets is not None:
    #    for b, b_set in enumerate(bootstrap_sets):
    #        perf_curves = np.dstack([b_roc.performance_curve().values for b_roc in b_set])
    #        pod_range = np.nanpercentile(perf_curves[:, 0], ci, axis=1)
    #        sr_range = np.nanpercentile(1 - perf_curves[:, 1], ci, axis=1)
    #        pod_poly = np.concatenate((pod_range[1], pod_range[0, ::-1]))
    #        sr_poly = np.concatenate((sr_range[1], sr_range[0, ::-1]))
    #        pod_poly[np.isnan(pod_poly)] = 0
    #        sr_poly[np.isnan(sr_poly)] = 1
    #        plt.fill(sr_poly, pod_poly, alpha=0.5, color=colors[b])
    #for r, roc_obj in enumerate(roc_objs):
    #    perf_data = roc_obj.performance_curve()
    #    plt.plot(1 - perf_data["FAR"], perf_data["POD"], marker=markers[r], color=colors[r], label=obj_labels[r])

    #plt.plot(1 - df['FAR'], df['POD'], marker=markers, color=colors, label=obj_labels)
    plt.plot(df['SR'], df['POD'], marker=markers, color=colors, label=obj_labels)

    cbar = plt.colorbar(csi_contour)
    cbar.set_label(csi_label)
    plt.xlabel(xlabel, fontsize=label_fontsize)
    plt.ylabel(ylabel, fontsize=label_fontsize)
    plt.xticks(ticks, fontsize=tick_fontsize)
    plt.yticks(ticks, fontsize=tick_fontsize)
    plt.title(title, fontsize=title_fontsize)
    plt.legend(**legend_params)
    plt.savefig(fname, dpi=dpi, bbox_inches='tight')
    plt.close()
