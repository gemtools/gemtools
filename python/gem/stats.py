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

## set locale to get thousands separator easily in 2.6
locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

## default colors
__colors = ['#3182bd', '#6baed6', '#9ecae1', '#c6dbef', '#e6550d', '#fd8d3c']


# the tick formatter for percentage and reads
def __percent_format(x, pos=0):
    """Tick formatter to render as
    percentage
    """
    return '%1.0f%%' % (100 * x)


def __reads_format(x, pos=0):
    """Render number of reads with thousands
    separator
    """
    return locale.format("%.0f", x, True)


def write_general_stats(data, out_dir, paired=True):
    """Generate general stats plot and
    save it to the given out_dir

    data    -- the stats data
    out_dir -- the output directory
    paired  -- paired reads
    """

    num_blocks = (float)(data["num_blocks"])
    num_split_maps = data["splits_profile"]["num_mapped_with_splitmaps"]
    num_mapped = data["num_mapped"]

    if paired:
        num_blocks = num_blocks / 2

    num_unmapped = num_blocks - num_mapped
    total = (float)(num_mapped + num_unmapped)

    fig = plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(111)

    ax1.tick_params(top=False, bottom=False)
    plt.xticks(visible=False)
    plt.xlim([0, 3])

    ax1.bar([1.5], num_mapped + num_unmapped, 1, color=__colors[5], label="Unmapped")
    ax1.bar([1.5], num_mapped, 1, color=__colors[1], label="Mapped")
    ax1.bar([1.5], num_split_maps, 1, color=__colors[0], label="Split-Mapped")

    ax1.yaxis.set_major_formatter(FuncFormatter(__reads_format))
    ax1.set_ylabel("Reads")

    ax2 = ax1.twinx()
    ax2.set_ylabel("Percent")
    ax2.yaxis.set_major_formatter(FuncFormatter(__percent_format))

    # legend
    lgd = ax1.legend(loc='lower center', bbox_to_anchor=(0.5, -0.1), ncol=3)

    # add descriptor lines
    percent = (num_split_maps / total)
    ax1.axhline(y=num_split_maps, color=__colors[0])
    ax1.text(0.1, num_split_maps, "Split-Maps %s (%.1f%%)" % (locale.format("%.0f", num_split_maps, True), (percent * 100.0)), verticalalignment='bottom')

    percent = (num_mapped / total)
    ax1.axhline(y=num_mapped, color=__colors[1])
    ax1.text(0.1, num_mapped, "Mapped %s (%.1f%%)" % (locale.format("%.0f", num_mapped, True), (percent * 100.0)), verticalalignment='top')

    percent = (num_unmapped / total)
    ax1.axhline(y=total, color=__colors[5])
    ax1.text(0.1, total, "Unmapped %s (%.1f%%)" % (locale.format("%.0f", num_unmapped, True), (percent * 100.0)), verticalalignment='top')

    fig.savefig('%s/general.png' % (out_dir), bbox_extra_artists=(lgd,), bbox='tight')

if __name__ == "__main__":
    f = "python/testdata/teststats.stats.json"
    # load stats
    with open(f) as of:
        data = json.load(of)

    out = "test_report"
    write_general_stats(data, out)
