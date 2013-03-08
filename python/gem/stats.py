#!/usr/bin/env python

import json
import locale

#import matplotlib stuff
try:
    from pylab import *
    from matplotlib.ticker import FuncFormatter
except Exception, e:
    # todo: trigger warning about missing imports
    raise


def write_general_stats(data, out_dir):
    def percent_format(x, pos=0):
        return '%1.0f%%' % (100 * x)

    def reads_format(x, pos=0):
        #return '{0:,.0f}'.format(x)
        return locale.format("%.0f", x, True)

    paired = True

    num_blocks = (float)(data["num_blocks"])
    num_split_maps = data["splits_profile"]["num_mapped_with_splitmaps"]
    num_mapped = data["num_mapped"]

    if paired:
        num_blocks = num_blocks / 2

    num_unmapped = num_blocks - num_mapped
    total = (float)(num_mapped + num_unmapped)

    fig, ax1 = plt.subplots(figsize=(10, 10))
    tick_params(top=False, bottom=False)
    xticks(visible=False)
    xlim([0, 3])
    colors = ['#3182bd', '#6baed6', '#9ecae1', '#c6dbef', '#e6550d', '#fd8d3c']

    plt.bar([1.5], num_mapped + num_unmapped, 1, color=colors[5], label="Unmapped")
    plt.bar([1.5], num_mapped, 1, color=colors[1], label="Mapped")
    plt.bar([1.5], num_split_maps, 1, color=colors[0], label="Split-Mapped")

    ax1.yaxis.set_major_formatter(FuncFormatter(reads_format))
    ax2 = ax1.twinx()
    ax2.set_ylabel("Percent")
    ax2.yaxis.set_major_formatter(FuncFormatter(percent_format))
    ax1.legend(loc="upper left", bbox_to_anchor=(1.1, 1., 1., .0))

    # add lines
    locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

    num_split_maps_p = (num_split_maps / total)
    ax1.axhline(y=num_split_maps, color=colors[0])
    text(0.1, num_split_maps_p, "Split-Maps %s (%.1f%%)" % (locale.format("%.0f", num_split_maps, True), (num_split_maps_p * 100.0)), verticalalignment='bottom')

    num_split_maps_p = (num_mapped / total)
    ax1.axhline(y=num_mapped, color=colors[1])
    text(0.1, num_split_maps_p + 0.01, "Mapped %s (%.1f%%)" % (locale.format("%.0f", num_mapped, True), (num_split_maps_p * 100.0)), verticalalignment='top')

    num_split_maps_p = (num_unmapped / total)
    ax1.axhline(y=total, color=colors[5])
    ax1.text(0.1, total, "Unmapped %s (%.1f%%)" % (locale.format("%.0f", num_unmapped, True), (num_split_maps_p * 100.0)), verticalalignment='top')

    savefig('%s/general.png' % (out_dir))

if __name__ == "__main__":
    f = "python/testdata/teststats.stats.json"
    # load stats
    with open(f) as of:
        data = json.load(of)

    out = "test_report"
    write_general_stats(data, out)
