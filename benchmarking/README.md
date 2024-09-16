# Benchmarking IPS6

Nextflow provides some built in options for assessing the operation of `IPS6`, including generating a HTML report. However, these reports are limited to presenting the resource utilisation from only single run, and can only be generated if a run is successful. Consequently, we have packaged a simple benchmarking script into IPS6 to allow for assessing the task duration and resource usage across multiple runs, and customised grouping of the data. For example, you may wish to clearly see differences in performance with altering the batch size. 

## Using `Nextflow`

### Report

Run `IPS6` with the `-with-report [file name]` Nextflow parameter to generate a summary html report for the current run. You can find more in the [Nextflow docs](https://www.nextflow.io/docs/latest/tracing.html#execution-report).

### Timeline

Run `IPS6` with the `-with-timeline [file name]` Nextflow parameter to generate a visual representation of the operational timeline for the current run. You can find more in the [Nextflow docs](https://www.nextflow.io/docs/latest/tracing.html#timeline-report).

## Using IPS6 benchmarking scripts

1. **Setup**

2. **Build a trace file as part of your InterProScan runs**

3. **Update the benchmarking config file**

4. **Run the benchmarking**

