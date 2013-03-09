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






# # -*- coding: utf-8 -*-
# # <nbformat>3.0</nbformat>

# # <codecell>

# import json
# from pylab import *

# # <codecell>

# with open("C171CACXX_471D.stats.json") as f:
#     data = json.load(f)

# # <codecell>

# ## general
# from matplotlib.ticker import FuncFormatter

# def percent_format(x, pos=0):
#     return '%1.0f%%'%(100*x)

# def reads_format(x, pos=0):
#     return '{0:,.0f}'.format(x)

# paired = True
# pair_label = ""

# num_blocks = (float) (data["num_blocks"])
# num_split_maps = data["splits_profile"]["num_mapped_with_splitmaps"]
# num_mapped = data["num_mapped"]

# if paired:
#     pair_label = " (Pairs)"
#     num_blocks = num_blocks / 2

# num_unmapped = num_blocks - num_mapped
# total = (float)(num_mapped+num_unmapped)

# fig, ax1 = plt.subplots()
# tick_params(top=False, bottom=False)
# xticks(visible=False)
# xlim([0,3])
# colors = ['#3182bd', '#6baed6', '#9ecae1', '#c6dbef', '#e6550d', '#fd8d3c']

# plt.bar([1.5], num_mapped+num_unmapped, 1, color=colors[5], label="Unmapped")
# plt.bar([1.5], num_mapped, 1, color=colors[1], label="Mapped")
# plt.bar([1.5], num_split_maps, 1, color=colors[0], label="Split-Mapped")

# ax1.yaxis.set_major_formatter(FuncFormatter(reads_format))
# ax2 = ax1.twinx()
# ax2.set_ylabel("Percent")
# ax2.yaxis.set_major_formatter(FuncFormatter(percent_format))
# ax1.legend(loc="upper left", bbox_to_anchor=(1.1, 1., 1., .0))

# # add lines
# num_split_maps_p = (num_split_maps / total)
# ax1.axhline(y=num_split_maps, color=colors[0])
# text(0.1,num_split_maps_p, "Split-Maps {0:,.0f} ({1:.1f}%)".format(num_split_maps, (num_split_maps_p * 100.0)), verticalalignment='bottom')

# num_split_maps_p = (num_mapped / total)
# ax1.axhline(y=num_mapped, color=colors[1])
# text(0.1,num_split_maps_p-0.02, "Mapped {0:,.0f} ({1:.1f}%)".format(num_mapped, (num_split_maps_p * 100.0)), verticalalignment='top')

# num_split_maps_p = (num_unmapped / total)
# ax1.axhline(y=total, color=colors[5])
# ax1.text(0.1,total, "Unmapped {0:,.0f} ({1:.1f}%)".format(num_unmapped, (num_split_maps_p * 100.0)), verticalalignment='top')

# # <codecell>

# ## mmap ranges
# mmap_ranges_values = data["mmap"]
# mmap_ranges_labels = data["mmap_description"]
# alignments = (float) (data["num_alignments"])
# rest = [alignments - sum(mmap_ranges_values)]
# [rest.append(d) for d in mmap_ranges_values]
# rest = [(d/alignments)*100.0 for d in rest]
# grid(True)
# bar(xrange(9), rest, color="#3366FF", align="center", alpha=0.5)
# xticks(xrange(10), ("0", "1", "(1,5]", "(5,10]", "(10,50]", "(50,100]", "(100,500]", "(500,1000]", "(1000,inf]"), rotation=45, )
# title("Multi-Map Ranges")
# xlabel("Ranges")
# ylabel("% Alignments")



# # <codecell>

# ## unique ranges
# uniq_ranges_values = data["uniq"]
# rest = [uniq_ranges_values[-1]]
# [rest.append(d) for d in uniq_ranges_values[:7]]
# alignments = (float)(sum(rest))
# rest = [(d/alignments)*100.0 for d in rest]
# grid(True)
# bar(xrange(8), rest, color="#3366FF", align="center", alpha=0.5)
# xticks(xrange(8), ("X", "0", "1", "2", "3", "(3,10]", "(10,50]", "(50,inf]"), rotation=45, )
# title("Unique Ranges")
# xlabel("Ranges")
# ylabel("% Alignments")

# # <codecell>

# ## error profiles
# def plot_e_profile(da, _title, xlab="% Errors", ylab="% Alignments"):
#     labels = ("0%", "1%", "2%", "3%", "4%", "5%", "6%", "7%", "8%", "9%", "10%", "(10,20]%", "(20,50]%", "(50,100]%")
#     num_maps = (float) (data["num_maps"])
#     rest = da
#     rest = [(d/num_maps)*100.0 for d in rest]
#     grid(True)
#     bar(xrange(14), rest, color="#3366FF", alpha=0.5)
#     xticks(xrange(14), labels, rotation=45, )
#     title(_title)
#     xlabel(xlab)
#     ylabel(ylab)
#     ylim([0,100])

# figure(figsize=(10, 10))
# subplot(2,2,1)
# plot_e_profile(data["maps_profile"]["mismatches"], "Mismatch Profile")
# subplot(2,2,2)
# plot_e_profile(data["maps_profile"]["insertion_length"], "Insertion Profile")
# subplot(2,2,3)
# plot_e_profile(data["maps_profile"]["deletion_length"], "Deletion Profile")
# subplot(2,2,4)
# plot_e_profile(data["maps_profile"]["levenshtein"], "Levenshtein Profile")

# # <codecell>

# # error positions
# max_len = data["max_length"]
# error_events = (float)(data["maps_profile"]["total_errors_events"])
# da = data["maps_profile"]["error_position"]
# da = [(d/error_events)*100.0 for d in da]
# plot(da[:max_len],color="#3366FF")
# fill_between(range(max_len), da[:max_len], color="#3366FF", alpha=0.5)
# xlim([0,max_len])
# grid(True)
# title("Error events")
# xlabel("Position")
# ylabel("% Alignments")


# # <codecell>

# ## insert size ranges
# inss = data["maps_profile"]["inss"]
# labels = ["(-inf, 0)", "(-100, 0)", "(0, 100]", "(100, 200]", "(200, 300]", "(300, 400]", "(400, 500]", "(500, 600]", "(600, 700]", "(700, 800]", "(800, 900]", "(900, 1000]", "(1000, 2000]", "(2000, 5000]", "(5000, 10000]", "(1000, inf]"]
# num_maps = (float) (data["num_maps"])
# rest = [(d/num_maps)*100.0 for d in inss]

# grid(True)
# bar(xrange(16), rest, color="#3366FF", alpha=0.5)
# xticks(xrange(16), labels, rotation=45, )
# title("Insert sizes")
# xlabel("Ranges")
# ylabel("% Alignments")

# # <codecell>

# # quality mismatches
# offset = 33
# max_len = 41

# total = (float)(data["maps_profile"]["total_mismatches"])
# da = data["maps_profile"]["qual_score_misms"][offset:offset+max_len]
# da = [(d/total)*100.0 for d in da]

# plot(da,color="#3366FF")
# fill_between(range(max_len), da[:max_len], color="#3366FF", alpha=0.5)
# xlim([0,max_len])
# ylim([0,max(da)+3])
# grid(True)
# title("Quality value for mismatch")
# xlabel("Quality")
# ylabel("% Mismatches")

# # <codecell>

# offset = 33
# max_len = 41

# total = (float)(data["maps_profile"]["total_errors_events"])
# da = data["maps_profile"]["qual_score_errors"][offset:offset+max_len]
# da = [(d/total)*100.0 for d in da]

# plot(da,color="#3366FF")
# fill_between(range(max_len), da[:max_len], color="#3366FF", alpha=0.5)
# xlim([0,max_len])
# ylim([0,max(da)+3])
# grid(True)
# title("Quality value for errors")
# xlabel("Quality")
# ylabel("% Errors")

# # <codecell>

# offset = 33
# max_len = 41

# total = (float)(data["maps_profile"]["total_errors_events"])
# da = data["maps_profile"]["qual_score_errors"][offset:offset+max_len]
# da = [(d/total)*100.0 for d in da]

# plot(da,color="#FF5533", label="Errors")
# fill_between(range(max_len), da[:max_len], color="#FF5533", alpha=0.5)

# total = (float)(data["maps_profile"]["total_mismatches"])
# da = data["maps_profile"]["qual_score_misms"][offset:offset+max_len]
# da = [(d/total)*100.0 for d in da]

# plot(da,color="#3366FF", label="Mismatches")
# fill_between(range(max_len), da[:max_len], color="#3366FF", alpha=0.5)
# legend(loc="upper left")


# xlim([0,max_len])
# ylim([0,max(da)+3])
# grid(True)
# title("Quality value for errors")
# xlabel("Quality")
# ylabel("% Errors/Mismatches")

# # <codecell>

# # transisitons
# total = (float)(data["maps_profile"]["total_mismatches"])
# da = data["maps_profile"]["misms_transition"][:5*5]
# da = [(d/total)*100.0 for d in da]
# da = array([da[5*x:5*x+5] for x in range(5)])
# column_labels = list('ACGTN')
# row_labels = list('ACGTN')
# fig, ax = plt.subplots()
# pcolor(da, cmap=plt.cm.Blues, edgecolors="black")
# colorbar()
# # put the major ticks at the middle of each cell
# ax.set_xticks(np.arange(da.shape[0])+0.5, minor=False)
# ax.set_yticks(np.arange(da.shape[1])+0.5, minor=False)

# # want a more natural, table-like display
# ax.invert_yaxis()
# ax.xaxis.tick_top()

# ax.set_xticklabels(row_labels, minor=False, family='monospace')
# ax.set_yticklabels(column_labels, minor=False,family='monospace')
# tick_params(top=False, left=False, right=False)
# for x in range(5):
#     for y in range(5):
#         text(0.5+x,0.5+y, "%.2f%%" % (da[y,x]), horizontalalignment='center', verticalalignment='center')
# ylabel("Transitions")

# # <codecell>

# # transisitons
# total = (float)(data["maps_profile"]["total_mismatches"])
# da = data["maps_profile"]["misms_1context"][:(4*4*4)*5]
# da = [(d/total)*100.0 for d in da]
# da = array([da[5*x:5*x+5] for x in range(4*4*4)])
# row_labels = list('ACGTN')
# column_labels = list([a+b+c for b in "ACGT" for a in "ACGT" for c in "ACGT"])
# #figure(num=None, figsize=(3, 6), dpi=80, facecolor='w', edgecolor='k')
# fig, ax = plt.subplots(figsize=(5, 15))
# pcolor(da, cmap=plt.cm.Blues, edgecolors="black")
# colorbar()
# # put the major ticks at the middle of each cell
# ax.set_xticks(np.arange(da.shape[1])+0.5, minor=False)
# ax.set_yticks(np.arange(da.shape[0])+0.5, minor=False)

# # want a more natural, table-like display
# ax.invert_yaxis()
# ax.xaxis.tick_top()
# tick_params(top=False, left=False, right=False)

# ax.set_xticklabels(row_labels, minor=False, family='monospace')
# ax.set_yticklabels(column_labels, minor=False, family='monospace')
# for x in range(5):
#     for y in range((4*4*4)):
#         text(0.5+x,0.5+y, "%.2f%%" % (da[y,x]), horizontalalignment='center', verticalalignment='center')
# ylabel("Transitions")

# # <codecell>

# ## junctions
# ## mmap ranges
# sp = data["splits_profile"]
# max_len = data["max_length"]
# total_junctions = (float)(sp["total_junctions"])

# figure(figsize=(15,10))

# subplot2grid((2,3), (0,0))
# da = [(d/total_junctions)*100.0 for d in sp["num_junctions"]]
# labels = ["[1]", "[2]", "[3]", "(3, inf)"]
# pie(da, labels=labels, autopct="%1.1f%%", shadow=False, colors=['#3182bd', '#6baed6', '#9ecae1', '#c6dbef', '#e6550d', '#fd8d3c'])
# title("Number of Junctions")

# subplot2grid((2,3), (0,1))
# da = [(d/total_junctions)*100.0 for d in sp["length_junctions"]]
# labels = ["[0,100]", "(100, 1000]", "(1000, 5000]", "(5000, 10000]", "(10000, 50000]", "(50000, inf)"]
# pie(da, labels=labels, autopct="%1.1f%%", shadow=False, colors=['#3182bd', '#6baed6', '#9ecae1', '#c6dbef', '#e6550d', '#fd8d3c'])
# title("Junction Lengths")

# subplot2grid((2,3), (0,2))
# pe = [sp["pe_rm_rm"], sp["pe_sm_rm"], sp["pe_sm_sm"]]
# sum_pe = (float) (sum(pe))
# da = [(d/sum_pe)*100.0 for d in pe]
# labels = ["RM+RM", "SM+RM", "SM+SM"]
# pie(da, labels=labels, autopct="%1.1f%%", shadow=False, colors=['#3182bd', '#6baed6', '#9ecae1', '#c6dbef', '#e6550d', '#fd8d3c'])
# title("Pair combinations")


# subplot2grid((2,3), (1,0), colspan=3)
# da = sp["junction_position"]
# da = [(d/total_junctions)*100.0 for d in da]
# plot(da[:max_len],color="#3366FF")
# fill_between(range(max_len), da[:max_len], color="#3366FF", alpha=0.5)
# xlim([0,max_len])
# grid(True)
# title("Junction Positions")
# xlabel("Position")
# ylabel("% Junctions")



# # <codecell>


