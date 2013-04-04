#!/usr/bin/env python

from string import Template
import errno
import json
import locale
import os
import re
import shutil
import zipfile

#import matplotlib stuff
__plotlib_avail = False
try:
    from pylab import *
    from matplotlib.ticker import FuncFormatter
    __plotlib_avail = True
except Exception, e:
    pass

## set locale to get thousands separator easily in 2.6
locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

## default colors
__colors = ['#3182bd', '#6baed6', '#9ecae1', '#c6dbef', '#e6550d', '#fd8d3c']
__default_color = __colors[1]

__default_template = '''
<html>
<head>
    <link rel="stylesheet" type="text/css" href="style.css">
    <style type="text/css">
        /*RESET DEFAULTS*/
        html, body, div, span, applet, object, iframe,
        h1, h2, h3, h4, h5, h6, p, blockquote, pre,
        a, abbr, acronym, address, big, cite, code,
        del, dfn, em, font, ins, kbd, q, s, samp,
        small, strike, strong, sub, sup, tt, var,
        dl, dt, dd, ol, ul, li,
        fieldset, form, label, legend,
        table, caption, tbody, tfoot, thead, tr, th, td {
            border: 0;
            font-family: inherit;
            font-size: 100%;
            font-style: inherit;
            font-weight: inherit;
            margin: 0;
            padding: 0;
            vertical-align: baseline;
        }

        :focus { /* remember to define focus styles! */
            outline: 0;
        }

        body, input, textarea {
            color: #373737;
            font: 15px "Helvetica Neue", Helvetica, Arial, sans-serif;
            font-weight: 300;
            line-height: 1.625;
        }
        body {
            background: #FFF;
        }
        /* Alignment */
        .alignleft {
            display: inline;
            float: left;
            margin-right: 2em;
        }
        .alignright {
            display: inline;
            float: right;
        }
        .aligncenter {
            clear: both;
            display: block;
            margin-left: auto;
            margin-right: auto;
        }

        strong {
            font-weight: bold;
        }
        a {
            color: #1C231C;
            text-decoration: none;
        }
        a:focus,
        a:active,
        a:hover {
            text-decoration: underline;
        }

        .page{
            margin: 2em;
        }
        /*HEADLINES*/
        .header{
            border-bottom: 1px solid #DDD;
        }
        .header h1{
            font-size: 2.5em;
            text-align: center;
        }

        h1{
            font-size: 1.6em;
            text-decoration: underline;
        }

        .data{
            border-bottom: 1px solid #DDD;
        }

        /*PLOTS*/
        .plot{
            width: 60em;
        }

        .general_plot{
            width: 40em;
            margin: -35px 0 0 0;
        }

        .transitions_plot{
            width: 30em;
        }

        .transitions_1context_plot{
            width: 35em;
        }
    </style>
</head>
<body>
    <div class="page">
        <div class="header">
            <h1>${name} mapping stats</h1>
        </div>
        <div class="general data">
            <h1>General Stats</h1>
            <img src="general.png" class="general_plot plot"/>
            <div class="alignleft">
                <table>
                    <tr>
                        <td>#Reads</td>
                        <td>${reads}</td>
                    </tr>
                    <tr>
                        <td>Reads length (min, avg, max)</td>
                        <td>${min}, ${avg}, ${max}</td>
                    </tr>
                    <tr>
                        <td>Reads Mapped</td>
                        <td>${mapped} (${mapped_p})</td>
                    </tr>
                    <tr>
                        <td>#Alignments</td>
                        <td>${alignments}</td>
                    </tr>
                    <tr>
                        <td>#Maps</td>
                        <td>${maps} (${maps_p} map/alg)</td>
                    </tr>
                </table>
            </div>
        </div>
        <div class="errorprofile data">
            <h1>Error Profile</h1>
            <img src="error_profile.png" class="errors_plot plot"/>
        </div>
        <div class="ranges data">
            <h1>Ranges</h1>
            <img src="ranges.png" class="ranges_plot plot"/>
        </div>
        <div class="transitions data">
            <h1>Transitions</h1>
            <div>
                <img src="transitions.png" class="transitions_plot plot"/>
            </div>
            <div>
                <img src="transitions_1context.png" class="transitions_1context_plot plot"/>
            </div>
        </div>
        <div class="junctions data">
            <h1>Junctions Profile</h1>
            <img src="junctions_profile.png" class="junctions_plot plot"/>
        </div>
    </div>
</body>
</html>
'''


def __zipfolder(foldername, target_dir):
    zipobj = zipfile.ZipFile(foldername + '.zip', 'w', zipfile.ZIP_DEFLATED)
    rootlen = len(target_dir) + 1
    bd = os.path.basename(target_dir)
    for base, dirs, files in os.walk(target_dir):
        for file in files:
            fn = os.path.join(base, file)
            zipobj.write(fn, os.path.join(bd, fn[rootlen:]))


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

    fig.savefig('%s/general.png' % (out_dir), bbox='tight')


def write_mmaps_and_uniq_ranges(data, out_dir):
    fig = plt.figure(figsize=(20, 15))

    # get data, add to list and transform to percent
    mmap_ranges_values = data["mmap"]
    alignments = (float)(data["num_alignments"])
    rest = [alignments - sum(mmap_ranges_values)]
    [rest.append(d) for d in mmap_ranges_values]
    rest = [(d / alignments) * 100.0 for d in rest]

    subplots_adjust(hspace=0.3)
    ## mmap ranges plot
    subplot2grid((3, 2), (0, 0))
    grid(True)
    bar(xrange(9), rest, color=__default_color, align="center")
    plt.xticks(xrange(10), ("0", "1", "(1,5]", "(5,10]", "(10,50]", "(50,100]", "(100,500]", "(500,1000]", "(1000,inf]"), rotation=45, )
    ylim([0, 100])
    title("Multi-Map Ranges")
    xlabel("Ranges")
    ylabel("% Alignments")

    # get data and transform to percentage
    uniq_ranges_values = data["uniq"]
    rest = [uniq_ranges_values[-1]]
    [rest.append(d) for d in uniq_ranges_values[:7]]
    alignments = (float)(sum(rest))
    rest = [(d / alignments) * 100.0 for d in rest]
    # plot
    ## unique ranges
    subplot2grid((3, 2), (0, 1))
    grid(True)
    ylim([0, 100])
    bar(xrange(8), rest, color=__default_color, align="center")
    xticks(xrange(8), ("X", "0", "1", "2", "3", "(3,10]", "(10,50]", "(50,inf]"), rotation=45, )
    title("Unique Ranges")
    xlabel("Ranges")
    ylabel("% Alignments")

    subplot2grid((3, 2), (1, 0), colspan=2)
    inss = data["maps_profile"]["inss"]
    labels = ["(-inf, 0)", "(-100, 0)", "(0, 100]", "(100, 200]", "(200, 300]", "(300, 400]", "(400, 500]", "(500, 600]", "(600, 700]", "(700, 800]", "(800, 900]", "(900, 1000]", "(1000, 2000]", "(2000, 5000]", "(5000, 10000]", "(10000, inf]"]
    num_maps = (float)(data["num_maps"])
    rest = [(d / num_maps) * 100.0 for d in inss]
    grid(True)
    ylim([0, 100])
    bar(xrange(16), rest, color=__default_color)
    xticks(xrange(16), labels, rotation=45, )
    title("Insert sizes")
    xlabel("Ranges")
    ylabel("% Alignments")

    subplot2grid((3, 2), (2, 0), colspan=2)
    inss = data["maps_profile"]["inss_fine_grain"]
    num_maps = (float)(data["num_maps"])
    rest = [(d / num_maps) * 100.0 for d in inss]
    grid(True)
    xlim([-1000, 10000])
    bar(xrange(-1000, 10000, 10), rest, color=__default_color)

    # labels = ["<1000"]
    # ticks = [0]
    # for x in range(-900, 10000, 100):
    #     labels.append(str(x))
    #     ticks.append(x + 901)
    # xticks(ticks, labels, rotation=0, )
    title("Insert sizes (fine grained)")
    xlabel("Ranges")
    ylabel("% Alignments")

    fig.savefig('%s/ranges.png' % (out_dir), bbox_inches='tight')


def write_error_profiles(data, out_dir, offset=33):
    # ## error profiles
    def plot_e_profile(da, _title, xlab="% Errors", ylab="% Alignments"):
        labels = ("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "(10,20]", "(20,50]", "(50,100]")
        num_maps = (float)(data["num_maps"])
        rest = da
        rest = [(d / num_maps) * 100.0 for d in rest]
        grid(True)
        bar(xrange(14), rest, color=__default_color)
        xticks(xrange(14), labels[:len(rest)], rotation=45, )
        title(_title)
        xlabel(xlab)
        ylabel(ylab)
        ylim([0, 100])
        xlim([0, len(rest)])

    figure(figsize=(15, 20))
    subplots_adjust(hspace=0.5)
    subplot2grid((4, 2), (0, 0))
    plot_e_profile(data["maps_profile"]["mismatches"], "Mismatch Profile")
    subplot2grid((4, 2), (0, 1))
    plot_e_profile(data["maps_profile"]["insertion_length"], "Insertion Profile")
    subplot2grid((4, 2), (1, 0))
    plot_e_profile(data["maps_profile"]["deletion_length"], "Deletion Profile")
    subplot2grid((4, 2), (1, 1))
    plot_e_profile(data["maps_profile"]["levenshtein"], "Levenshtein Profile")

    # plot errors and mismatches
    subplot2grid((4, 2), (2, 0), colspan=2)
    max_len = 41
    total = (float)(data["maps_profile"]["total_errors_events"])
    da = data["maps_profile"]["qual_score_errors"][offset:offset + max_len]
    da = [(d / total) * 100.0 for d in da]
    plot(da, color="#FF5533", label="Errors")
    fill_between(range(max_len), da[:max_len], color="#FF5533", alpha=0.5)
    total = (float)(data["maps_profile"]["total_mismatches"])
    da = data["maps_profile"]["qual_score_misms"][offset:offset + max_len]
    da = [(d / total) * 100.0 for d in da]
    plot(da, color=__default_color, label="Mismatches")
    fill_between(range(max_len), da[:max_len], color=__default_color, alpha=0.5)
    ylim(bottom=0)
    title("Quality Errors/Mismatches Profile")
    xlabel("Quality Score")
    ylabel("Errors/Mismatches")
    legend(loc="upper left")

    subplot2grid((4, 2), (3, 0), colspan=2)
    max_len = data["max_length"]
    error_events = (float)(data["maps_profile"]["total_errors_events"])
    da = data["maps_profile"]["error_position"]
    da = [(d / error_events) * 100.0 for d in da]
    plot(da[:max_len], color=__default_color)
    fill_between(range(max_len), da[:max_len], color=__default_color, alpha=0.5)
    ylim(bottom=0)
    xlim([0, max_len])
    grid(True)
    title("Error events")
    xlabel("Position")
    ylabel("% Alignments")

    savefig('%s/error_profile.png' % (out_dir), bbox_inches='tight')


def __exclude_zero(data, labels, delta=0.1):
    """Exclude values <= delta from data and labels
    """
    d = []
    l = []
    for i, x in enumerate(data):
        if x > delta:
            d.append(x)
            l.append(labels[i])
    return (d, l)


def write_junctions_profile(data, out_dir):
    # ## junctions
    sp = data["splits_profile"]
    max_len = data["max_length"]
    total_junctions = (float)(sp["total_junctions"])
    if total_junctions == 0:
        total_junctions = 1.0

    figure(figsize=(20, 12))
    #subplots_adjust( hspace=0, wspace=0 )

    subplot2grid((2, 3), (0, 0))
    da = [(d / total_junctions) * 100.0 for d in sp["num_junctions"]]
    labels = ["[1]", "[2]", "[3]", "(3, inf)"]
    da, labels = __exclude_zero(da, labels)
    pie(da, labels=labels, autopct="%1.1f%%", shadow=False, colors=__colors)
    title("Number of Junctions")

    subplot2grid((2, 3), (0, 1))
    da = [(d / total_junctions) * 100.0 for d in sp["length_junctions"]]
    labels = ["[0,100]", "(100, 1000]", "(1000, 5000]", "(5000, 10000]", "(10000, 50000]", "(50000, inf)"]
    da, labels = __exclude_zero(da, labels)
    pie(da, labels=labels, autopct="%1.1f%%", shadow=False, colors=__colors)
    title("Junction Lengths")

    subplot2grid((2, 3), (0, 2))
    pe = [sp["pe_rm_rm"], sp["pe_sm_rm"], sp["pe_sm_sm"]]
    sum_pe = (float)(sum(pe))
    if sum_pe == 0:
        sum_pe = 1.0

    da = [(d / sum_pe) * 100.0 for d in pe]
    labels = ["RM+RM", "SM+RM", "SM+SM"]
    da, labels = __exclude_zero(da, labels)
    pie(da, labels=labels, autopct="%1.1f%%", shadow=False, colors=__colors)
    title("Pair combinations")

    subplot2grid((2, 3), (1, 0), colspan=3)
    da = sp["junction_position"]
    da = [(d / total_junctions) * 100.0 for d in da]
    max_len = min(max_len, len(da))
    plot(da[:max_len], color="black")
    fill_between(range(max_len), da[:max_len], color=__default_color)
    xlim([0, max_len])
    ylim([0, max(da) + 1])
    grid(True)
    title("Junction Positions")
    xlabel("Position")
    ylabel("% Junctions")
    savefig('%s/junctions_profile.png' % (out_dir), bbox_inches='tight')


def write_transitions(data, out_dir):

    figure(figsize=(10, 10))
    # transisitons
    total = (float)(data["maps_profile"]["total_mismatches"])
    da = data["maps_profile"]["misms_transition"][:5 * 5]
    da = [(d / total) * 100.0 for d in da]
    da = array([da[5 * x:5 * x + 5] for x in range(5)])
    column_labels = list('ACGTN')
    row_labels = list('ACGTN')

    subplot2grid((2, 1), (0, 0))
    ax = gca()
    pcolor(da, cmap=plt.cm.Blues, edgecolors="black")
    colorbar()
    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(da.shape[0]) + 0.5, minor=False)
    ax.set_yticks(np.arange(da.shape[1]) + 0.5, minor=False)

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    ax.set_xticklabels(row_labels, minor=False, family='monospace')
    ax.set_yticklabels(column_labels, minor=False, family='monospace')
    tick_params(top=False, left=False, right=False)
    for x in range(5):
        for y in range(5):
            text(0.5 + x, 0.5 + y, "%.2f%%" % (da[y, x]), horizontalalignment='center', verticalalignment='center')
    ylabel("Transitions")
    savefig('%s/transitions.png' % (out_dir), bbox_inches='tight')

    # context 1 transitions
    raw = data["maps_profile"]["misms_1context"]
    raw = [(d / total) * 100.0 for d in raw]

    da = []

    def __get_index(a, b, c, i):
        return ((((a * 5 + b) * 5) + c) * 5 + i)

    for b in range(4):
        for a in range(4):
            for c in range(4):
                for i in range(5):
                    da.append(raw[__get_index(a, b, c, i)])

    da = array([da[5 * x:5 * x + 5] for x in range(4 * 4 * 4)])
    row_labels = list('ACGTN')
    column_labels = list([a + b + c for b in "ACGT" for a in "ACGT" for c in "ACGT"])

    figure(figsize=(10, 30))
    pcolor(da, cmap=plt.cm.Blues, edgecolors="black")
    colorbar()
    ax = gca()
    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(da.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(da.shape[0]) + 0.5, minor=False)

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    tick_params(top=False, left=False, right=False)

    ax.set_xticklabels(row_labels, minor=False, family='monospace')
    ax.set_yticklabels(column_labels, minor=False, family='monospace')
    for x in range(5):
        for y in range((4 * 4 * 4)):
            text(0.5 + x, 0.5 + y, "%.2f%%" % (da[y, x]), horizontalalignment='center', verticalalignment='center')
    ylabel("Transitions")

    savefig('%s/transitions_1context.png' % (out_dir), bbox_inches='tight')


def write_template(data, out, paired=True, name=None):
    if paired:
        avg_length = data["total_bases_aligned"] / float(data["num_mapped"] * 2)
    else:
        avg_length = data["total_bases_aligned"] / float(data["num_mapped"])

    num_blocks = (float)(data["num_blocks"])
    num_split_maps = data["splits_profile"]["num_mapped_with_splitmaps"]
    num_mapped = data["num_mapped"]

    if paired:
        num_blocks = num_blocks / 2

    num_unmapped = num_blocks - num_mapped
    total = (float)(num_mapped + num_unmapped)

    if name is None:
        name = os.path.basename(out)

    tmpl = Template(__default_template).safe_substitute({
        "name": name,
        "reads": data["num_blocks"],
        "min": data["mapped_min_length"],
        "max": data["mapped_max_length"],
        "avg": "%.0f" % avg_length,
        "mapped": data["num_mapped"],
        "mapped_p": "%.2f%%" % ((num_mapped / total) * 100.0),
        "alignments": data["num_blocks"],
        "maps": data["num_maps"],
        "maps_p": "%.3f" % (data["num_maps"] / float(data["num_mapped"])),
    })
    with open("%s/index.html" % (out), 'w') as f:
        f.write(tmpl)


def create_report(input_file, output_name, paired=True, extract=False, name=None):
    """Create a stats report from a json stats file and store it in a zip
    file using the given output name.

    Parameters
    ----------
    input_file: string or file handle
        The input file either as a string pointing to the json report file
        or as an open readable stream.
    output_name: stirng
        The output prefix is used to create a directory that hosts the
        html report
    paired: bool
        Set this to false if the input is single end
    extract: bool
        Set this to true to keep the directory next to the zip file
    name:
        Name of the dataset
    """
    if not __plotlib_avail:
        raise Exception("""
Matplotlib could not be imported. We need the matplotlib library
to render the stats report! 

Please install matplotlib and its dependencies. For example:

pip install numpy
pip install matplotlib
""")
    if input_file is None:
        raise ValueError("No input file specified")
    if output_name is None:
        raise ValueError("No output name specified")

    # load input data
    of = None
    if isinstance(input_file, basestring):
        of = open(input_file, 'r')

    data = json.load(of)
    of.close()

    # guess name
    if name is None:
        m = re.match("(.*)(\.stats\.(all|best)\.json$)", input_file)
        if m:
            name = m.group(1)
        else:
            idx = input_file.rfind(".")
            if idx > 0:
                name = input_file[:idx]

    # create output directory
    try:
        os.makedirs(output_name)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(output_name):
            pass
        else:
            raise

    # create plots
    write_general_stats(data, output_name, paired=paired)
    write_mmaps_and_uniq_ranges(data, output_name)
    write_error_profiles(data, output_name)
    write_junctions_profile(data, output_name)
    write_transitions(data, output_name)

    # print the data to the folder
    with open("%s/stats.json" % output_name, 'w') as of:
        json.dump(data, of, indent=2)

    # write the html template
    write_template(data, output_name, paired=paired, name=name)

    # zip the folder
    __zipfolder(output_name, output_name)

    # remove folder
    if not extract:
        shutil.rmtree(output_name)
