import statistics

from collections import namedtuple
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


PROCESSES = [
    'GET_ORFS',
    'PARSE_SEQUENCE',
    'ANTIFAM_HMMER_RUNNER', 'ANTIFAM_HMMER_PARSER',
    'CDD_RUNNER', 'CDD_PARSER', 'CDD_POSTPROCESS',
    'COILS_RUNNER', 'COILS_PARSER',
    'FUNFAM_HMMER_RUNNER', 'FUNFAM_HMMER_PARSER', 'FUNFAM_CATH_RESOLVE_HITS', 'FUNFAM_FILTER_MATCHES',
    'GENE3D_HMMER_RUNNER', 'GENE3D_HMMER_PARSER', 'GENE3D_CATH_RESOLVE_HITS', 'GENE3D_ADD_CATH_SUPERFAMILIES',  'GENE3D_FILTER_MATCHES',
    'HAMAP_HMMER_RUNNER', 'HAMAP_HMMER_PARSER', 'HAMAP_POST_PROCESSER', 'HAMAP_FILTER_MATCHES',
    'MOBIDB_RUNNER', 'MOBIDB_PARSER',
    'NCBIFAM_HMMER_RUNNER', 'NCBIFAM_HMMER_PARSER',
    'PANTHER_HMMER_RUNNER', 'PANTHER_HMMER_PARSER', 'PANTHER_POST_PROCESSER', 'PANTHER_FILTER_MATCHES',
    'PFAM_HMMER_RUNNER', 'PFAM_HMMER_PARSER', 'PFAM_FILTER_MATCHES',
    'PHOBUIS_RUNNER', 'PHOBIUS_PARSER',
    'PIRSF_HMMER_RUNNER', 'PIRSF_HMMER_PARSER', 'PIRSF_FILTER_MATCHES',
    'PIRSR_HMMER_RUNNER', 'PIRSR_HMMER_PARSER', 'PIRSR_FILTER_MATCHES',
    'PRINTS_RUNNER', 'PRINTS_PARSER',
    'PROSITE_PATTERNS_RUNNER', 'PROSITE_PATTERNS_PARSER',
    'PROSITE_PROFILES_RUNNER', 'PROSITE_PROFILES_PARSER',
    'SFLD_HMMER_RUNNER', 'SFLD_HMMER_PARSER', 'SFLD_POST_PROCESSER', 'SFLD_FILTER_MATCHES',
    'SIGNALP_RUNNER', 'SIGNALP_PARSER',
    'SMART_HMMER2_RUNNER', 'HMMER2_PARSER', 'SMART_FILTER_MATCHES',
    'SUPERFAMILY_HMMER_RUNNER',  'SUPERFAMILY_POST_PROCESSER', 'SUPERFAMILY_FILTER_MATCHES',
    'AGGREGATE_RESULTS', 'AGGREGATE_PARSED_SEQS',
    'XREFS:ENTRIES', 'XREFS:PAINT_ANNOTATIONS', 'XREFS:GOTERMS', 'XREFS:PATHWAYS',
    'WRITE_RESULTS',
]


def plot_total_runtime(
    df: pd.DataFrame,
    group: str|None,
    outdir: Path,
    fig_formats: set,
    group_order: list
):
    """Plot the total run time with one series per group"""
    total_run_times = []
    grp_col = group if group else "Groups"
    for grp in set(df[grp_col]):
        for run in set(df['Run']):
            rows = df[(df[grp_col] == grp) & (df['Run'] == run)]
            total_run_times.append([grp, sum(rows['raw_realtime'])])

    total_rt_df = pd.DataFrame(total_run_times, columns=[grp_col, 'Real Runtime (s)'])
    total_rt_df = total_rt_df.sort_values(grp_col)

    fig, ax = plt.subplots()
    g = sns.stripplot(
        data=total_rt_df,
        x=grp_col,
        y='Real Runtime (s)',
        color='black',
        size=4,
        ax=ax,
        hue_order=group_order,
        dodge=True,
    )
    g = sns.boxplot(
        data=total_rt_df,
        x=grp_col,
        y='Real Runtime (s)',
        fliersize=0,
        ax=ax,
        hue_order=group_order,
        order=group_order,
    )
    plt.xticks(rotation=90)
    g.set_ylabel(ylabel='Real Runtime (s)')
    g.set_xlabel(xlabel=grp_col)
    g.set(title='Total InterProScan Runtime')

    for fig_format in fig_formats:
        plt.savefig((outdir / f"total_runtime.{fig_format}"), bbox_inches='tight', format=fig_format)

    plt.clf()


def plot_process_runtime(
    df: pd.DataFrame,
    group: str|None,
    outdir: Path,
    fig_formats: set,
    group_order: list
):
    """Plot the runtimes, broken down by process. This will help to identify
    slower running (or bottleneck) processes."""
    x = "raw_process"
    y = "raw_realtime"
    grp_col = group if group else "Groups"
    df = df.sort_values(grp_col)
    hue = grp_col if len(set(df[grp_col])) > 1 else None
    process_order = [_ for _ in PROCESSES if _ in set(df["raw_process"])]

    # If there is only one group then colour code by process
    if len(set(df[grp_col])) > 1:
        hue_order = group_order
    else:
        hue_order = process_order

    # Write the processes out in the correct order
    df.sort_values(
        by="raw_process",
        key=lambda column: column.map(lambda e: process_order.index(e)),
        inplace=True
    )

    if len(process_order) > 30:
        fig, ax = plt.subplots(figsize=(20, 5))
    else:
        fig, ax = plt.subplots()

    sns.set(font_scale=0.7)

    if hue:
        g = sns.stripplot(
            data=df,
            x=x,
            y=y,
            hue=hue,
            hue_order=hue_order,
            size=3,
            ax=ax,
            linewidth=1,
            dodge=True
        )
    else:
        g = sns.stripplot(
            data=df,
            x=x,
            y=y,
            color='black',
            size=3,
            ax=ax,
            linewidth=1,
            dodge=True
        )
    g = sns.boxplot(
        data=df,
        x=x,
        y=y,
        hue=hue,
        fliersize=0,
        ax=ax,
        hue_order=hue_order,
        width=0.75,
        linewidth=0.75,
    )
    plt.xticks(rotation=90)
    g.set_ylabel(ylabel='Real Runtime (s)')
    g.set_xlabel(xlabel='Process')
    g.set_title('Runtime Per Process')

    for fig_format in fig_formats:
        plt.savefig((outdir / f"process_runtime.{fig_format}"), bbox_inches='tight', format=fig_format)

    plt.clf()


def plot_process_runtime_piechart(
    df: pd.DataFrame,
    group: str|None,
    outdir: Path,
    fig_formats: set
):
    """Plot the runtimes, broken down by process. This will help to identify
    slower running (or bottleneck) processes."""
    processes = []
    Process = namedtuple('Process', ['name', 'average', 'sd'])
    for process_name in set(df['raw_process']):
        p_rows = df[df['raw_process'] == process_name]
        p_mean = round(statistics.mean(p_rows['raw_realtime']), 1)
        p_sd = round(statistics.stdev(p_rows['raw_realtime']), 1) if len(p_rows['raw_realtime']) >= 2 else "na"
        processes.append(
            Process(
                process_name,
                p_mean,
                p_sd
            )
        )

    segment_values = [_.average for _ in processes]
    labels = [f"{_.name} ({_.sd} SD)" for _ in processes]

    plt.pie(
        segment_values,
        labels=labels,
        autopct='%.0f%%',
        pctdistance=0.8,
        labeldistance=1.2,
        textprops={'fontsize': 5}
    )

    for fig_format in fig_formats:
        plt.savefig(
            (outdir / f"process_runtime_piechart.{fig_format}"),
            bbox_inches='tight',
            format=fig_format
        )

    pie_df = pd.DataFrame([[p.name, p.average, p.sd] for p in processes], columns=['Process', 'Mean Runtime', 'SD'])
    pie_df.to_csv((outdir / "pie_chart_values.csv"))

    plt.clf()


def plot_overall_summary(
    df: pd.DataFrame,
    group: str|None,
    outdir: Path,
    fig_formats: set,
    fig_name: str,
    x: str,
    x_label: str,
    y: str,
    y_label: str,
    title: str,
    hue_order: str,
    group_order: list
):
    """Plot all data across all processes grouped together,
    only separate by the user defined 'groups'."""
    process_order = [_ for _ in PROCESSES if _ in set(df["raw_process"])]
    grp_col = group if group else "Groups"
    hue = grp_col if len(set(df[grp_col])) > 1 else None

    if hue_order == 'groups':  # do not break up the data by process
        hue_order = group_order
        x_axis_order = group_order
        legend_title = grp_col
        bbox_to_anchor = (.5, -0.5)

    else:  # break up the data by process, then sub-group by user defined group
        if len(group_order) > 1:
            hue_order = group_order
        else:
            hue_order = process_order
        x_axis_order = process_order
        df.sort_values(
            by="raw_process",
            key=lambda column: column.map(lambda e: process_order.index(e)),
            inplace=True
        )
        legend_title = "Process"
        bbox_to_anchor = (.5, -1.5)

    if len(x_axis_order) > 30:
        if hue:
            if len(df[hue]) > 1:
                fig, ax = plt.subplots(figsize=(30, 5))
        else:
            fig, ax = plt.subplots(figsize=(12, 5))
    else:
        fig, ax = plt.subplots()

    df[x] = df[x].astype(str)

    sns.set(font_scale=0.75)

    if hue:
        g = sns.stripplot(
            data=df,
            x=x,
            y=y,
            hue=hue,
            hue_order=hue_order,
            size=3,
            ax=ax,
            linewidth=1,
            dodge=True
        )
    else:
        g = sns.stripplot(
            data=df,
            x=x,
            y=y,
            color='black',
            size=3,
            ax=ax,
            linewidth=1,
            dodge=True
        )
    g = sns.boxplot(
        data=df,
        x=x,
        y=y,
        fliersize=0,
        ax=ax,
        hue=hue,
        hue_order=hue_order,
        order=x_axis_order,
        linewidth=0.75
    )
    plt.xticks(rotation=90)
    g.set_ylabel(ylabel=y_label)
    g.set_xlabel(xlabel=x_label)
    g.set(title=title)

    for fig_format in fig_formats:
        plt.savefig((outdir / f"{fig_name}.{fig_format}"), bbox_inches='tight', format=fig_format)

    plt.clf()
