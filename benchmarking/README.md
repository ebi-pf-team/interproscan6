# Benchmarking IPS6

Nextflow provides some built in options for assessing the operation of `IPS6`, including generating a HTML report. However, these reports are limited to presenting the resource utilisation from only single run, and can only be generated if a run is successful. Consequently, we have packaged a simple benchmarking script into IPS6 to allow for assessing the task duration and resource usage across multiple runs, and customised grouping of the data. For example, you may wish to clearly see differences in performance with altering the batch size. 

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
