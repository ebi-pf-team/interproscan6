# Benchmarking IPS6

Nextflow provides some built in options for assessing the operation of `IPS6`, including generating a HTML report. However, these reports are limited to presenting the resource usage from only a single run, and can only be generated if a run is successful. Consequently, we have packaged a simple benchmarking script into IPS6 to enable assessing the task duration and resource usage across multiple runs, and customised grouping of the data. For example, you may wish to clearly see differences in performance with altering the batch size, the number of CPUs or amount of memory allocated. 

## Using IPS6 benchmarking scripts

1. Install all necessary third party packages listed in `benchmarking/requirements.txt`.

2. Build a trace file as part of your InterProScan runs.

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

3. Update or create the benchmarking JSON config file

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

4. Run the benchmarking

```bash
# running from the root of the IPS6 project dir
# and using the benchmarking/tracefiles.json file
python3 benchmarking/benchmark_ips6.py benchmarking/tracefiles.json
```

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
