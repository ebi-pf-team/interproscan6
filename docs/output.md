# `InterProScan` Output

## `InterProScan-5` Output:

For each signature:
* `signature`: dict summarising the signature
	* `accession`: signature accession
	* `name`: shortname from InterPro
	* `description`: long name from InterPro
	* `signatureLibraryRelease`: dict {`library`: application name, `version`: release version number}
	* `entry`: ???
* `locations` : location of domains:
    * List of dicts
    * Each dict is a domain hit
        * `start` start location of domain
        * `end` end location of domain
        * `hmmStart` start location of hmm
        * `hmmEnd` end location of hmm
        * `evalue`: domain evalue
        * `score`: bitscore
        * `envelopesStart`: start of envelop (the envelope is the region or range of the protein sequence were the domain may be located) / `envelopeEnd`
        * `location-fragments` -- list of dicts, one dict per fragment --> start, end and dc-status
        * `sites`: list of dicts, one dict per site (rememeber a site can have multiple residues/locations)
            * `description`: str
            * `numLocations`: basically the number of residues, most locations are a single residue
            * `label`: ???
            * `group`: ????
            * `hmmStart`: ???? -- surely not relevant here?
            * `hmmEnd`: ???? -- surely not relevant here?      
            * `siteLocations`: list, one dict per location:
                * `start`: int
                * `end`: int
                * `residue`: char
* `evalue`: overall, full sequence evalue
* `score`: overall, full sequence bit-score
* `model-ac`: signature accession --> but this is already stored under the `signature`:`accession` keys.

Example output file: `testsfld.json`
example: A0A0H3E4R3_BACA1