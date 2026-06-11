# Licensed applications

Some InterProScan analyses depend on third-party applications that must be licensed and installed separately.

You do not need a separate license to run InterProScan itself or to analyse your own sequences with the bundled analyses.

In InterProScan 6, the licensed applications are:

- `DeepTMHMM`
- `Phobius`
- `SignalP-Euk`
- `SignalP-Prok`

These analyses remain unavailable until you install the third-party software and configure InterProScan to use the extracted installation directory.

!!! tip

    You only need to download and configure the licensed applications you plan to run.

## Activation steps

To enable a licensed application in InterProScan:

1. Obtain a license from the software provider.
2. Download and extract the application archive.
3. Create a Nextflow config file, for example `licensed.conf`.
4. Add the absolute path to the extracted directory under `params.appsConfig`.
5. Run InterProScan with `-c licensed.conf`.

## Obtaining the software

### DeepTMHMM 1.0

Request a standalone copy of DeepTMHMM 1.0 by emailing <licensing@biolib.com>.

After receiving the archive:

```bash
unzip -q DeepTMHMM-v1.0.zip
```

### Phobius 1.01

Download Phobius 1.01 from [Erik Sonnhammer's website](https://software.sbc.su.se/phobius.html), then extract it:

```bash
tar -zxf phobius101_linux.tgz
```

!!! info

    Phobius skips sequences containing the non-standard or ambiguous residues `O`, `B`, `Z`, or `J`. Other InterProScan analyses continue to run normally for those sequences.

### SignalP 6.0

SignalP 6.0 is available in two model variants:

- Full (`slow`)
- Distilled (`fast`)

Use the full model when you need the most accurate region boundaries and have sufficient time and memory available.

Licenses and downloads are available from the [DTU SignalP 6.0 page](https://services.healthtech.dtu.dk/services/SignalP-6.0/).

After downloading the archive:

```bash
tar -zxf signalp-6.0i.fast.tar.gz
```

## Configuration

Store the application paths in a config file such as `licensed.conf`:

```groovy
params {
  appsConfig {
    deeptmhmm {
      dir = "/full/path/to/DeepTMHMM"
    }
    phobius {
      dir = "/full/path/to/phobius"
    }
    signalp_euk {
      dir = "/full/path/to/signal6p_fast"
      mode = "fast"
    }
    signalp_prok {
      dir = "/full/path/to/signal6p_fast"
      mode = "fast"
    }
  }
}
```

!!! danger

    For SignalP, set `mode = "fast"` for the distilled model, and `mode = "slow"` or `mode = "slow-sequential"` for the full model. `slow` runs the full model in parallel and requires substantially more memory. `slow-sequential` runs the same full model sequentially, uses roughly the same memory as `fast`, and is about six times slower.

Run InterProScan with the default analyses plus all licensed deep-learning analyses:

```bash
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.1 \
  -profile docker \
  -c licensed.conf \
  --input /path/to/sequences.faa \
  --run-ml  # (1)!
```

1. DeepTMHMM, SignalP-Euk, and SignalP-Prok use deep-learning models that are computationally intensive. They must be enabled explicitly with `--run-ml`, or selected explicitly with `--applications`.

Run only the licensed applications explicitly:

```bash
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.1 \
  -profile docker \
  -c licensed.conf \
  --input /path/to/sequences.faa \
  --applications deeptmhmm,phobius,signalp_euk,signalp_prok
```

!!! info

    Running both `signalp_euk` and `signalp_prok` executes SignalP twice, once with the eukaryotic model and once with the prokaryotic model. When possible, run only the model that matches your dataset.
