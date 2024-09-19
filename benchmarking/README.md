# Benchmarking IPS6

Nextflow provides some built in options for assessing the operation of `IPS6`, including generating a HTML report. However, these reports are limited to presenting the resource usage from only a single run, and can only be generated if a run is successful. Consequently, we have packaged a simple benchmarking script into IPS6 to enable assessing the task duration and resource usage across multiple runs, and customised grouping of the data. For example, you may wish to clearly see differences in performance with altering the batch size, the number of CPUs or amount of memory allocated. 

## Using `Nextflow`

### Report

Run `IPS6` with the `-with-report [file name]` Nextflow parameter to generate a summary html report for the current run. You can find more in the [Nextflow docs](https://www.nextflow.io/docs/latest/tracing.html#execution-report).

### Timeline

Run `IPS6` with the `-with-timeline [file name]` Nextflow parameter to generate a visual representation of the operational timeline for the current run. You can find more in the [Nextflow docs](https://www.nextflow.io/docs/latest/tracing.html#timeline-report).

## Using IPS6 benchmarking scripts

1. **Setup**

Install all necessary third party packages listed in `benchmarking/requirements.txt`.

2. **Build a trace file as part of your InterProScan runs**

Include the following code in the `nextflow.config` file and rename the resulting trace file
`ips6.trace.txt` after each run.

```groovy
trace {
    enabled = true
    fields = 'task_id,process,realtime,duration,status,cpus,%cpu,time,memory,%mem,rss,peak_rss,submit,start,complete,queue'
    file = 'ips6.trace.txt'
    overwrite = true   
}
```

You can find more information the columns that can be included in the trace in the [Nextflow documentation](https://www.nextflow.io/docs/latest/tracing.html#trace-report). 

At a minimum, the `IPS6` benchmarking scripts require the following columns:
* process
* realtime
* rss
* peak_rss

3. **Update the benchmarking config file**

Update the `tracefile.json` file (or create a new JSON file) keyed by the name of each group of analyses,
e.g. the number of CPU assigned to the group or the batch size used, and valued by a list of string 
representations of paths to `IPS6` trace files.

For example, if you wanted to assess the impact of altering the batch size on performance:

```json
{
    "500": [
        "benchmarking/24.08.27.batch.500.report.3.tsv",
        "benchmarking/24.08.26.batch.500.report.2.tsv",
        "benchmarking/24.08.26.batch.500.report.1.tsv"
    ],
    "1000": [
        "benchmarking/24.08.27.batch.1000.report.3.tsv",
        "benchmarking/24.08.26.batch.1000.report.2.tsv",
        "benchmarking/24.08.26.batch.1000.report.1.tsv"
    ],
    "5000": [
        "benchmarking/24.08.27.batch.5000.report.3.tsv",
        "benchmarking/24.08.26.batch.5000.report.2.tsv",
        "benchmarking/24.08.26.batch.5000.report.1.tsv"
    ]
}
```

4. **Run the benchmarking**

The onl required argument for running the benchmarking is the path to the JSON file listing the 
paths to the trace files.

```bash
# running from the root of the IPS6 project dir
# and using the benchmarking/tracefiles.json file
python3 benchmarking/benchmark_ips6.py benchmarking/tracefiles.json
```

You can print a help message to prin out the argument options:
```bash
python3 benchmarking/benchmark_ips6.py --help
```

By default, the benchmarking will label the groupings as 'Groups' on the resulting plot axes and 
legends. You can name the groupings using the `--group_name` flag and providing the name you 
wish to be assigned to the axes and legends, e.g. `--group_name "Batch Sizes"`, or `--group_name "Number of CPU"`.

By default, the resulting figures are only written out in `PDF` format. Use the `--format` flag to 
list the desired file outputs. Accepted outputs: png, pdf, and svg. For example to generate svg and 
png files use `--format png,svg`.

By default the trace file writes the in human readable format, but can be configured to write the raw
values. If this is the case, include the `--raw` flag in the `benchmark_ips6.py` command.

By default, the output figures will be written to the current working directory. To write the files 
to a desired output directory use the `--outdir` flag and provide the path for the output dir. The 
scripts will build all necessary parent directories for the output dir.

If you wish to perform further analyses on the data, use the `--save_data` flag to configure 
`benchmark_ips6.py` to write out the dataframe it generates to a CSV file in the output dir.

## Output:

Each run of `benchmark_ips6.py` will produce the following figures (note, references to 'group' refers to the keys in the input JSON file, each key represents a different 'group'):

1. `total_runtime.*` - Shows the total run time of IPS6 per group in the input JSON file
2. `process_runtime.*` - Shows the total run time per process in IPS6
3. `process_runtime_piechart.*` - Shows the percentage of the total runtime contributed by each process
4. `pie_chart_values.csv` - Contains the data used to build the `process_runtime_piechart.*` figure. If many processes are included the legends in the pie chart can often overlap. Use this CSV file to plot the pie chart (or alternative chart).
5. `overall_memory_usage.*` - Plots the overall memory usage per group in the input JSON file
6. `overall_max_memory_usage.*` - Plots the overall maximum memory used per group in the input JSON file
7. `memory_per_process.*` - Plots the memory usage per process (and per group if multiple groups are defined in the input JSON file)
8. `max_memory_per_process.*` - Plots the maximum memory usage per process (and per group if multiple groups are defined in the input JSON file)

Each box and whisker plot is overlaid by a strip plot with each point of the strip plot representing the value
from a single run.
