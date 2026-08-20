# Matches API

InterProScan can retrieve precomputed matches from the InterPro Matches API instead of recomputing every analysis locally.

This can make runs much faster when your input proteins correspond to sequences in UniParc that have already been annotated by the InterPro production pipeline and exposed through the Matches API.

## How it works

InterPro annotates UniParc protein sequences with InterProScan as part of its production pipeline.

For each InterPro release, typically every eight weeks, those production annotations are published through the Matches API.

When the Matches API is enabled, InterProScan queries the API for precomputed annotations for the input sequences and reuses them when available.

Sequences that are not available through the API are still executed locally.

Use `--no-matches-api` to disable this behaviour and force all sequences to be annotated locally.

!!! info

    There are currently no plans to add DeepTMHMM annotations to the Matches API. TMbed is the selected predictor of transmembrane proteins in the InterPro production pipeline.

## Deploying your own Matches API instance

InterProScan can query any compatible Matches API server through `--matches-api-url`.

The reference server implementation is available on [GitHub](https://github.com/ProteinsWebTeam/interpro-matches-api).

To deploy your own instance, download the published Matches API data, start the reference server, and point InterProScan at the server's base URL.

### Download the API data

Fetch the latest archive and checksum:

```bash
curl -OJ https://ftp.ebi.ac.uk/pub/databases/interpro/releases/latest/matches-api-data.tar.gz
curl -OJ https://ftp.ebi.ac.uk/pub/databases/interpro/releases/latest/matches-api-data.tar.gz.md5
```

Verify and unpack it:

```bash
md5sum -c matches-api-data.tar.gz.md5
mkdir matches-api-data
tar -zxvf matches-api-data.tar.gz \
  --strip-components=1 \
  -C matches-api-data
```

### Start the API server

With Docker:

```bash
docker run --rm \
  -v $PWD/matches-api-data:/data \
  -e MATCHES_API_PATH=/data \
  -p 8000:8000 \
  interpro/matches-api:0.5.0
```

With Singularity:

```bash
singularity run \
  -B $PWD/matches-api-data:/data \
  --env "MATCHES_API_PATH=/data" \
  docker://interpro/matches-api:0.5.0
```

The server listens on port `8000`.

### Point InterProScan at your server

Set the Matches API base URL when launching InterProScan:

```bash
interproscan.sh --matches-api-url http://your-host:8000
```

InterProScan will query the server's `/matches` route automatically.

You can also test the API directly by POSTing up to 100 MD5 hashes to `/matches`:

```bash
curl -X POST "http://your-host:8000/matches" \
  -H 'Content-Type: application/json' \
  -d '{
    "md5": [
      "020DA9322D8466E699BDD584593749FC",
      "0204A4724F5991CB9B7E1013CBDA3367",
      "00C5D66EC4232D18A9639ECF3FC4BDDB"
    ]
  }'
```

## Configuration

The Matches API is enabled by default.

Relevant options:

- `--matches-api-url <URL>`: use an alternative Matches API endpoint
- `--no-matches-api`: disable the Matches API and run all analyses locally

See [Options / Parameters](options.md) for the full option reference.
