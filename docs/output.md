<style>
  :root {
    color: red;
  }

  /* .md-typeset table:not([class]) {
      display: table;
  } */

  .md-typeset table td:first-child,
  .md-typeset table th:first-child{
    width: 1%;
    min-width: 0;
    white-space: nowrap;
  }

  .md-typeset table td:last-child,
  .md-typeset table th:last-child{
    white-space: nowrap;
  }
</style>

# Output

InterProScan can produce `gff3`, `json`, `jsonl`, `tsv`, and `xml` outputs.

## GFF3

The GFF3 output follows the standard GFF3 format with InterProScan-specific headers and a FASTA tail section.

### Headers

```
##gff-version 3.1.26
##interpro-version 108.0
##interproscan-version 6.0.1
```

### Features

Features follow the 9 standard GFF3 columns:

1. Sequence identifier
2. Source of the annotation
3. Feature type
4. Feature start position (1-based, inclusive)
5. Feature end position (1-based, inclusive)
6. Score for the feature, or `.`
7. Strand, typically `+` or `-`, or `.` when not applicable
8. Phase: reading frame phase for CDS features, or `.` when not applicable
9. Attributes: semicolon-separated tag-value pairs

Common attributes include: `Name`, `Alias`, `Parent` (for ORF/nucleic mode), `Dbxref`, `Ontology_term`, `type`, and `representative`.

!!! ifno

    When [representative domains and families](representative_locations.md) are selected, all matches are still reported; the `representative` attribute is `true` only for selected representative features and `false` for the others. With `--skip-repr-locations`, it is `false` for all features.


### Sequences

`##FASTA` section with input sequences.

### Example

```gff3
##gff-version 3.1.26
##interpro-version 108.0
##interproscan-version 6.0.1
##sequence-region tr|A0A086JQP8|A0A086JQP8_TOXGO 1 847
tr|A0A086JQP8|A0A086JQP8_TOXGO	COILS	coiled_coil	485	526	.	.	.	Name=Coil;type=Region;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	PROSITE patterns	polypeptide_motif	106	115	.	.	.	Name=HSP90;Alias=PS00298;Dbxref=InterPro:IPR019805;type=Conserved_site;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	MobiDB-lite	polypeptide_region	51	88	.	.	.	Name=disorder_prediction;Alias=mobidb-lite;type=Region;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	MobiDB-lite	polypeptide_region	297	348	.	.	.	Name=disorder_prediction;Alias=mobidb-lite;type=Region;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	MobiDB-lite	polypeptide_region	299	309	.	.	.	Name=disorder_prediction;Alias=mobidb-lite;type=Region;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	MobiDB-lite	polypeptide_region	310	336	.	.	.	Name=disorder_prediction;Alias=mobidb-lite;type=Region;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	MobiDB-lite	polypeptide_region	799	847	.	.	.	Name=disorder_prediction;Alias=mobidb-lite;type=Region;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	MobiDB-lite	polypeptide_region	809	836	.	.	.	Name=disorder_prediction;Alias=mobidb-lite;type=Region;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	MobiDB-lite	polypeptide_region	837	847	.	.	.	Name=disorder_prediction;Alias=mobidb-lite;type=Region;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	CDD	polypeptide_domain	95	283	3.31697E-97	.	.	Name=HATPase_Hsp90-like;Alias=cd16927;Dbxref=InterPro:IPR020575;type=Domain;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	PROSITE profiles	polypeptide_domain	1	25	6.0	.	.	Name=PROKAR_LIPOPROTEIN;Alias=PS51257;type=Domain;representative=true
tr|A0A086JQP8|A0A086JQP8_TOXGO	SMART	polypeptide_domain	108	263	1.6E-4	.	.	Name=HATPase_c;Alias=SM00387;Dbxref=InterPro:IPR003594;type=Domain;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	PRINTS	polypeptide_motif	86	106	3.0E-66	.	.	Name=HEATSHOCK90;Alias=PR00775;Dbxref=InterPro:IPR020575;type=Family;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	PRINTS	polypeptide_motif	107	129	3.0E-66	.	.	Name=HEATSHOCK90;Alias=PR00775;Dbxref=InterPro:IPR020575;type=Family;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	PRINTS	polypeptide_motif	156	173	3.0E-66	.	.	Name=HEATSHOCK90;Alias=PR00775;Dbxref=InterPro:IPR020575;type=Family;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	PRINTS	polypeptide_motif	174	191	3.0E-66	.	.	Name=HEATSHOCK90;Alias=PR00775;Dbxref=InterPro:IPR020575;type=Family;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	PRINTS	polypeptide_motif	199	221	3.0E-66	.	.	Name=HEATSHOCK90;Alias=PR00775;Dbxref=InterPro:IPR020575;type=Family;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	PRINTS	polypeptide_motif	251	268	3.0E-66	.	.	Name=HEATSHOCK90;Alias=PR00775;Dbxref=InterPro:IPR020575;type=Family;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	PRINTS	polypeptide_motif	269	287	3.0E-66	.	.	Name=HEATSHOCK90;Alias=PR00775;Dbxref=InterPro:IPR020575;type=Family;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	PIRSR	polypeptide_region	76	495	5.3E-181	.	.	Name=PIRSR002583-1;Alias=PIRSR002583-1;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	HAMAP	polypeptide_region	83	790	25.80159	.	.	Name=HSP90;Alias=MF_00505;Dbxref=InterPro:IPR001404;type=Family;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	SUPERFAMILY	polypeptide_domain	356	641	2.88E-82	.	.	Name=SSF54211;Dbxref=InterPro:IPR020568;type=Homologous_superfamily;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	SUPERFAMILY	polypeptide_domain	85	292	7.11E-68	.	.	Name=SSF55874;Dbxref=InterPro:IPR036890;type=Homologous_superfamily;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	SUPERFAMILY	polypeptide_domain	665	788	4.32E-28	.	.	Name=SSF110942;Dbxref=InterPro:IPR037196;type=Homologous_superfamily;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	PIRSF	polypeptide_region	76	828	4.5E-240	.	.	Name=Hsp90;Alias=PIRSF002583;Dbxref=InterPro:IPR001404;type=Family;representative=true
tr|A0A086JQP8|A0A086JQP8_TOXGO	Pfam	polypeptide_region	265	815	5.5E-176	.	.	Name=HSP90;Alias=PF00183;Dbxref=InterPro:IPR001404;type=Family;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	Pfam	polypeptide_domain	110	225	3.6E-6	.	.	Name=HATPase_c_3;Alias=PF13589;type=Domain;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	NCBIFAM	polypeptide_region	86	786	1.6E-187	.	.	Name=PRK05218.1;Alias=NF003555;Dbxref=InterPro:IPR001404;type=Family;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	PANTHER	polypeptide_region	74	811	5.8E-253	.	.	Name=PTHR11528;Dbxref=InterPro:IPR001404;type=Family;representative=false
tr|A0A086JQP8|A0A086JQP8_TOXGO	CATH-Gene3D	polypeptide_domain	71	309	1.7E-86	.	.	Name=G3DSA:3.30.565.10;Dbxref=InterPro:IPR036890;type=Homologous_superfamily;representative=true
tr|A0A086JQP8|A0A086JQP8_TOXGO	CATH-Gene3D	polypeptide_domain	354	553	3.2E-71	.	.	Name=G3DSA:3.30.230.80;type=Homologous_superfamily;representative=true
tr|A0A086JQP8|A0A086JQP8_TOXGO	CATH-Gene3D	polypeptide_domain	554	642	1.5E-26	.	.	Name=G3DSA:3.40.50.11260;type=Homologous_superfamily;representative=true
tr|A0A086JQP8|A0A086JQP8_TOXGO	CATH-Gene3D	polypeptide_domain	643	793	1.4E-50	.	.	Name=G3DSA:1.20.120.790;Dbxref=InterPro:IPR037196;type=Homologous_superfamily;representative=true
tr|A0A086JQP8|A0A086JQP8_TOXGO	CATH-FunFam	polypeptide_domain	76	309	4.9E-102	.	.	Name=G3DSA:3.30.565.10:FF:000005;type=Region;representative=false
##FASTA
>tr|A0A086JQP8|A0A086JQP8_TOXGO
MSPAGRRTPKKLAFAALLLGVSVACTSSFFSASVSPSALWVAATETDAAEPLTAEEAPRS
LPIDESEKAAAPLTAEEQEAVQKSQESHQYQTEVSRLMDIIINSLYTQREVFLRELISNA
VDALEKVRFTALSHPEVLEPKKNLDIRIEFDADAKTLSIIDSGIGMTKQDLINNLGTVAK
SGTSNFLEAMAQGNDVNLIGQFGVGFYSAFLVADKVTVVSKNVEDDQHIWESSADAKFHV
AKDPRGNTLGRGTCVTLHLKEDATEFLNEWKLKDLTTRFSQFMSYPIYVRTSRTVTEEVP
IEDEEAETKDEDKDKDEDKDKDDVEVTEGDKDEKKDKPKTKKVEKKKDEWEQVNTQKAIW
LRPKEEIEEKEYNEFYKSVSKDWSDPLAHIHFSAEGEVEFKALLYIPKRAPSDIYSNYFD
KQTSVKVYVRRVLVADQFDDLLPKYLHFVKGVVDSDDLPLNVSREQLQQHKILNVISKKL
VRKTLDTMRKLSVDALKEREEMEKELEQEEDEAKKKELQKKLKEKSVYERFYDEFSRNLK
LGCYEDDTNRNKLLKLLRFHTSKSGPERSVTLESFVAKLPENQPNIYYAAGESAEQLMKA
PEMQIFLKKDIEVLFLLEAMDEPCIQRVMDFEGKKFVSIQKGDVQLDQTEEEKKTEKRLK
KAFEPLLSWWKKLLGEKVTKVEVSKRLVEAPCAVVASEWGYSAQMEKIMKTQTFADPRHV
RMMAGQKVFEINPHHRMIQYLLAQVQKEGDNVGSKEIEMARLLFEVAKLASGFEVEDPKD
VAASLYKAVAADLTLPTDEPMIAEYELPREEEDEKVGDEDAKDEEKNEEGEADEPEEKEH
TEKHDEL
```

## TSV

One line per reported match location. Columns are tab-separated:

1. Sequence identifier
2. MD5 checksum of the protein sequence
3. Sequence length
4. Analysis or member database name (e.g. `Pfam`, `SMART`)
5. Signature or model accession
6. Signature description, or `-`
7. Match start position (1-based, inclusive)
8. Match end position (1-based, inclusive)
9. Score or E-value reported by the analysis, or `-`
10. Match status, always `T`
11. Date (format: `dd-MM-yyyy`)
12. InterPro entry accession, or `-`
13. InterPro entry description, or `-`
14. GO terms and their source under parentheses, separated by pipes (`|`), e.g. `GO:0015934(InterPro)|GO:0022625(PANTHER)`
15. Patway cross-references, separated by pipes (`|`), or `0`

!!! info

    For nucleic-acid input (`--nucleic`), the sequence identifier is 
    reported as `<nucleotide_id>_<orf_id>`.

### Example

```txt
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	COILS	Coil	Coil	485	526	-	COILS	17-04-2026	-	-	-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	PROSITE patterns	PS00298	Heat shock hsp90 proteins family signature	106	115	-	PROSITE patterns	17-04-2026	IPR019805	Heat shock protein Hsp90, conserved site	-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	MobiDB-lite	mobidb-lite	Consensus disorder prediction	51	88	-	MobiDB-lite	17-04-2026	-	--	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	MobiDB-lite	mobidb-lite	Consensus disorder prediction	297	348	-	MobiDB-lite	17-04-2026	-	--	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	MobiDB-lite	mobidb-lite	Consensus disorder prediction	299	309	-	MobiDB-lite	17-04-2026	-	--	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	MobiDB-lite	mobidb-lite	Consensus disorder prediction	310	336	-	MobiDB-lite	17-04-2026	-	--	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	MobiDB-lite	mobidb-lite	Consensus disorder prediction	799	847	-	MobiDB-lite	17-04-2026	-	--	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	MobiDB-lite	mobidb-lite	Consensus disorder prediction	809	836	-	MobiDB-lite	17-04-2026	-	--	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	MobiDB-lite	mobidb-lite	Consensus disorder prediction	837	847	-	MobiDB-lite	17-04-2026	-	--	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	CDD	cd16927	Histidine kinase-like ATPase domain of human cytosolic Hsp90 and its homologs including Escherichia coli HtpG, and related domains	95	283	3.31697E-97	CDD	17-04-2026	IPR020575	Heat shock protein Hsp90, N-terminal	-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	PROSITE profiles	PS51257	Prokaryotic membrane lipoprotein lipid attachment site profile	1	25	6.0	PROSITE profiles	17-04-2026	-	-	-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	SMART	SM00387	Histidine kinase-like ATPases	108	263	1.6E-4	SMART	17-04-2026	IPR003594	Histidine kinase/HSP90-like ATPase domain	-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	PRINTS	PR00775	-	86	106	3.0E-66	PRINTS	17-04-2026	IPR020575	Heat shock protein Hsp90, N-terminal-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	PRINTS	PR00775	-	107	129	3.0E-66	PRINTS	17-04-2026	IPR020575	Heat shock protein Hsp90, N-terminal-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	PRINTS	PR00775	-	156	173	3.0E-66	PRINTS	17-04-2026	IPR020575	Heat shock protein Hsp90, N-terminal-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	PRINTS	PR00775	-	174	191	3.0E-66	PRINTS	17-04-2026	IPR020575	Heat shock protein Hsp90, N-terminal-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	PRINTS	PR00775	-	199	221	3.0E-66	PRINTS	17-04-2026	IPR020575	Heat shock protein Hsp90, N-terminal-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	PRINTS	PR00775	-	251	268	3.0E-66	PRINTS	17-04-2026	IPR020575	Heat shock protein Hsp90, N-terminal-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	PRINTS	PR00775	-	269	287	3.0E-66	PRINTS	17-04-2026	IPR020575	Heat shock protein Hsp90, N-terminal-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	PIRSR	PIRSR002583-1	-	76	495	5.3E-181	PIRSR	17-04-2026	-	-	-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	HAMAP	MF_00505	Chaperone protein HtpG [htpG]	83	790	25.80159	HAMAP	17-04-2026	IPR001404	Heat shock protein Hsp90 family	-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	SUPERFAMILY	SSF54211	Ribosomal protein S5 domain 2-like	356	641	2.88E-82	SUPERFAMILY	17-04-2026	IPR020568	Ribosomal protein uS5 domain 2-type superfamily	-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	SUPERFAMILY	SSF55874	ATPase domain of HSP90 chaperone/DNA topoisomerase II/histidine kinase	85	292	7.11E-68	SUPERFAMILY	17-04-2026	IPR036890	Histidine kinase/HSP90-like ATPase superfamily	-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	SUPERFAMILY	SSF110942	HSP90 C-terminal domain	665	788	4.32E-28	SUPERFAMILY	17-04-2026	IPR037196	HSP90, C-terminal domain	-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	PIRSF	PIRSF002583	Heat shock protein, HSP90/HTPG types	76	828	4.5E-240	PIRSF	17-04-2026	IPR001404	Heat shock protein Hsp90 family	-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	Pfam	PF00183	Hsp90 protein	265	815	5.5E-176	Pfam	17-04-2026	IPR001404	Heat shock protein Hsp90 family	-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	Pfam	PF13589	Histidine kinase-, DNA gyrase B-, and HSP90-like ATPase	110	225	3.6E-6	Pfam	17-04-2026	-	--	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	NCBIFAM	NF003555	molecular chaperone HtpG	86	786	1.6E-187	NCBIFAM	17-04-2026	IPR001404	Heat shock protein Hsp90 family	-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	PANTHER	PTHR11528	HEAT SHOCK PROTEIN 90 FAMILY MEMBER	74	811	3.9E-253	PANTHER	17-04-2026	IPR001404	Heat shock protein Hsp90 family	-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	CATH-Gene3D	G3DSA:3.30.565.10	Histidine kinase-like ATPase, C-terminal domain	71	309	1.7E-86	CATH-Gene3D17-04-2026	IPR036890	Histidine kinase/HSP90-like ATPase superfamily	-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	CATH-Gene3D	G3DSA:3.30.230.80	-	354	553	3.2E-71	CATH-Gene3D	17-04-2026	-	-	-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	CATH-Gene3D	G3DSA:3.40.50.11260	-	554	642	1.5E-26	CATH-Gene3D	17-04-2026	-	-	-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	CATH-Gene3D	G3DSA:1.20.120.790	Heat shock protein 90, C-terminal domain	643	793	1.4E-50	CATH-Gene3D17-04-2026	IPR037196	HSP90, C-terminal domain	-	-
tr|A0A086JQP8|A0A086JQP8_TOXGO	F647E1911ECD589A3D9A1B46899E62B7	847	CATH-FunFam	G3DSA:3.30.565.10:FF:000005	Heat shock protein 90	76	309	4.9E-102	CATH-FunFam	17-04-2026	-
```


## JSON

The JSON output provides a structured representation of InterProScan results, grouping each analysed sequence with its associated matches, signatures, locations, and InterPro annotations. This is the recommended output format for downstream processing because it is easier to parse reliably than the tab-delimited or XML formats.

The schema is aligned with the response returned by the [Matches API](matches_api.md), which makes it easier to move between local InterProScan runs and programmatic access through the API using the same overall data model.

### Example: protein sequence

```json
{
  "interproscan-version": "6.0.1",
  "interpro-version": "108.0",
  "results": [
    {
      "sequence": "MKAKEIREM...",  // (1)!
      "md5": "04A129E51F351B91B6373645AFDBCC58", // (2)!
      "matches": [
        {  // (3)!
          "signature": {  // (4)!
            "accession": "PF00831",
            "name": "Ribosomal_L29",
            "description": "Ribosomal L29 protein",
            "type": "Family",
            "signatureLibraryRelease": {
              "library": "Pfam",
              "version": "38.1"
            },
            "entry": {  // (5)!
              "accession": "IPR001854",
              "name": "Ribosomal_uL29",
              "description": "Large ribosomal subunit protein uL29",
              "type": "Family",
              "goXRefs": [
                {
                  "name": "structural constituent of ribosome",
                  "databaseName": "GO",
                  "category": "MOLECULAR_FUNCTION",
                  "id": "GO:0003735"
                },
                {
                  "name": "ribosome",
                  "databaseName": "GO",
                  "category": "CELLULAR_COMPONENT",
                  "id": "GO:0005840"
                },
                {
                  "name": "translation",
                  "databaseName": "GO",
                  "category": "BIOLOGICAL_PROCESS",
                  "id": "GO:0006412"
                }
              ],
              "pathwayXRefs": []
            }
          },
          "model-ac": "PF00831",
          "evalue": 2.5e-18,
          "score": 77.1,
          "source": "Pfam",
          "locations": [  // (6)!
            {
              "start": 4,
              "end": 59,
              "representative": true, // (7)!
              "hmmStart": 1,
              "hmmEnd": 56,
              "hmmLength": 57,
              "hmmBounds": "N_TERMINAL_COMPLETE",
              "evalue": 2.9e-18,
              "score": 76.9,
              "envelopeStart": 4,
              "envelopeEnd": 59,
              "location-fragments": [  // (8)!
                {
                  "start": 4,
                  "end": 59,
                  "dc-status": "CONTINUOUS"
                }
              ]
            }
          ]
        }
      ],
      "xref": [  // (9)!
        {
          "name": "sp|B8FES7|RL29_DESAL Large ribosomal subunit protein uL29 OS=Desulfatibacillum aliphaticivorans OX=218208 GN=rpmC PE=3 SV=1",
          "id": "sp|B8FES7|RL29_DESAL"
        }
      ]
    }
  ]
}
```

1. Amino-acid sequence of the input protein.
2. MD5 checksum of the protein sequence. This is a stable identifier for the exact sequence content and is used to detect identical sequences.
3. Match object. Each match represents one hit from a member database signature to the input protein.
4. Signature metadata provided by the member database.
5. InterPro entry to which the signature is integrated. If the value is `null`, the signature has not yet been integrated.
6. Match locations, that is, the sequence coordinates where the signature matched the protein. A single signature can produce more than one location on the same protein.
7. Boolean flag indicating whether this location was selected as the representative domain or family.
8. Domain fragments. CATH-Gene3D, CATH-FunFam, Pfam, and SUPERFAMILY can identify discontinuous domains, where one biological domain is split across separate sequence regions. In these cases, the coordinates of the individual fragments are reported in `location-fragments`.
9. Sequence cross-references and metadata. If multiple input records contain exactly the same protein sequence, InterProScan reports the matches once and lists the corresponding input identifiers in `xref`.

!!! info

    When [representative domains and families](representative_locations.md) are selected, all matches are still reported; the `representative` field is `true` only for selected representative locations and `false` for the others. With `--skip-repr-locations`, it is `false` for all locations.

### Example: nucleotide sequence

```json
{
  "interproscan-version": "6.0.1",
  "interpro-version": "108.0",
  "results": [
    {
      "sequence": "ATGAGGGATTCGCCCGATGAAGTCAGCGTCGACGAGCTGGTGAACATGGCCGTGGCCGGTGGCATCGACGAAGGAACGGCGTTGGACGCCTTACAGGGTAAGCTGGACCCGTACAAGGTAATGCGGGCTGCACACGAGGCCCGACTTAAGATCGTCGGTGAACACGTCACGTTCGTGGTGAACAGGAACATCAACTTCACCAACGTGTGCATTAACAGATGTCGATTCTGTGCGTTCCGGAGGGATCCGGACGACCCGGATGCTTACCGTATGACGCCGGAGGAGGTGGGCGAGCGGGCAGCGGAAGCCCGTGACGCTGGAGCTACGGAAGTATGTCTTCAGGGCGGACTGCATCCCGAGGCGACGTTTGAGTACTACCTGGAAATGTTGGACGAGATCAAGTCCCAAGCCCCGGACATCCACGTGCACGGGTACTCACCGATGGAGGTGAAGTACTGCGCCAAGCTGGCGGGAGAGGACATCGAAGACGTACTACGAGAGCTGAAGCGAGCCGGTCTCGATTCGATGCCCGGAACGGCCGCGGAGATATTCTCCCCTGAGGTGAGGAAGCGGCTATGTCCTGATAAGTTGGAAGCCGATGAGTGGGAACATATCATCAGGATCGCGCACGAGTTGGGAATTCCCACCACTTGTACTATGATGTACGGTCACATCGACTCACCGAGGGACTGGATCGACCACATGAAGCGGCTTCGAGGGATCCAAGAGGACACGGGAGGCTTCACGGAGTTCGTGCCGCTCTCCTTCGTACATTCGAACGCACCGATTTACCGACGAGGAGGGGCGCGACCCGGAGTATCGGGTATGACGGACGTACTCGTGCACGCTGTGGCCCGATTGTACTTCGGACCGTTGATTCCGAACATACAGGCTTCCTGGGTGAAGCTCGGAGTGAAGCTGGCTCAGATGACGCTGCACGCCGGGGCGAACGATCTAGGTGGCACCCTCATGGAAGAGAACATCTCCCGGGAGGCCGGAGCGACCGAGGGCGAGCAGCTCGAGCCCGAGGAGATAGTGGAGATCATTCGGGAGGCGGGCTTCACCCCCGTGCAGCGCACCACGCTCTACGAGCCGGTGAAGGTGTACTAA",
      "md5": "0A4E057B197C1182FD211A4BFF5271CE",
      "crossReferences": [
        {
          "id": "ENA|AAM02110|AAM02110.1 Methanopyrus kandleri AV19 Predicted enzyme related to thiamine biosynthesis enzyme ThiH",
          "name": "ENA|AAM02110|AAM02110.1"
        }
      ],
      "openReadingFrames": [
        {
          "start": 397,
          "end": 326,
          "strand": "ANTISENSE",
          "protein": {
            "sequence": "SRPTFPGSTQTSPRDAVRPEDILP",
            "md5": "3B52A60935909B9B722B11EBF380610C",
            "matches": [
              {
                "signature": {
                  "accession": "mobidb-lite",
                  "name": "disorder_prediction",
                  "description": "Consensus disorder prediction",
                  "type": "Region",
                  "signatureLibraryRelease": {
                    "library": "MobiDB-lite",
                    "version": "4.0"
                  },
                  "entry": null
                },
                "model-ac": "mobidb-lite",
                "source": "MobiDB-lite",
                "locations": [
                  {
                    "start": 1,
                    "end": 24,
                    "representative": false,
                    "location-fragments": [
                      {
                        "start": 1,
                        "end": 24,
                        "dc-status": "CONTINUOUS"
                      }
                    ],
                    "sequence-feature": null
                  },
                  {
                    "start": 1,
                    "end": 13,
                    "representative": false,
                    "location-fragments": [
                      {
                        "start": 1,
                        "end": 13,
                        "dc-status": "CONTINUOUS"
                      }
                    ],
                    "sequence-feature": "Polar"
                  },
                  {
                    "start": 15,
                    "end": 24,
                    "representative": false,
                    "location-fragments": [
                      {
                        "start": 15,
                        "end": 24,
                        "dc-status": "CONTINUOUS"
                      }
                    ],
                    "sequence-feature": "Polyampholyte"
                  }
                ]
              }
            ],
            "xref": [
              {
                "name": "orf1260 source=ENA|AAM02110|AAM02110.1 coords=397..326 length=24 frame=6 desc=Methanopyrus kandleri AV19 Predicted enzyme related to thiamine biosynthesis enzyme ThiH",
                "id": "orf1260"
              }
            ]
          }
        }
      ]
    }
  ]
}
```

## JSON Lines

JSON Lines (`.jsonl`) output provides the same data model as the standard JSON output, but writes one JSON object per line. Each line contains exactly one sequence result, so the `results` array always contains a single item.


## XML

The XML output provides a structured representation of InterProScan results, grouping each analysed sequence with its associated matches, signatures, locations, and InterPro annotations. The root element is `<results>` with version attributes, and the document contains `<protein>` results for protein input or `<nucleotide-sequence>` results with nested ORFs and translated proteins for nucleic input.

Compared with the InterProScan 5 XML schema, the InterProScan 6 XML format is simpler and closer to the JSON data model. It uses a `<results>` root, generic `<match>` and `<location>` elements, and a more consistent structure across analyses instead of the more heavily specialised InterProScan 5 schema.

### Example: protein sequence

```xml
<?xml version='1.0' encoding='UTF-8'?>
<results interproscan-version="6.0.1" interpro-version="108.0">
  <protein>
    <sequence md5="04A129E51F351B91B6373645AFDBCC58">MKAKEIREMGADEIRRKIDDSTQEMFNLRFQHATGQLENTARLNKTKKEVARLKTILKEVEQ</sequence>
    <xref id="sp|B8FES7|RL29_DESAL" name="sp|B8FES7|RL29_DESAL Large ribosomal subunit protein uL29 OS=Desulfatibacillum aliphaticivorans OX=218208 GN=rpmC PE=3 SV=1"/>
    <matches>
      <match evalue="2.5E-18" score="77.1" source="Pfam">
        <signature ac="PF00831" name="Ribosomal_L29" type="Family">
          <signature-library-release library="Pfam" version="38.1"/>
          <entry ac="IPR001854" desc="Large ribosomal subunit protein uL29" name="Ribosomal_uL29" type="Family">
            <go-xref category="MOLECULAR_FUNCTION" db="GO" id="GO:0003735" name="structural constituent of ribosome"/>
            <go-xref category="CELLULAR_COMPONENT" db="GO" id="GO:0005840" name="ribosome"/>
            <go-xref category="BIOLOGICAL_PROCESS" db="GO" id="GO:0006412" name="translation"/>
          </entry>
        </signature>
        <model-ac>PF00831</model-ac>
        <locations>
          <location start="4" end="59" representative="true" hmm-start="1" hmm-end="56" hmm-length="57" hmm-bounds="N_TERMINAL_COMPLETE" evalue="2.9E-18" score="76.9" env-start="4" env-end="59">
            <location-fragments>
              <fragment start="4" end="59" dc-status="CONTINUOUS"/>
            </location-fragments>
          </location>
        </locations>
      </match>
    </matches>
  </protein>
</results>
```

!!! info

    When [representative domains and families](representative_locations.md) are selected, all matches are still reported; the `representative` attribute is `true` only for selected representative locations and `false` for the others. With `--skip-repr-locations`, it is `false` for all locations.

### Example: nucleotide sequence

```xml
<?xml version='1.0' encoding='UTF-8'?>
<results interproscan-version="6.0.1" interpro-version="108.0">
  <nucleotide-sequence>
    <sequence md5="0A4E057B197C1182FD211A4BFF5271CE">ATGAGGGATTCGCCCGATGAAGTCAGCGTCGACGAGCTGGTGAACATGGCCGTGGCCGGTGGCATCGACGAAGGAACGGCGTTGGACGCCTTACAGGGTAAGCTGGACCCGTACAAGGTAATGCGGGCTGCACACGAGGCCCGACTTAAGATCGTCGGTGAACACGTCACGTTCGTGGTGAACAGGAACATCAACTTCACCAACGTGTGCATTAACAGATGTCGATTCTGTGCGTTCCGGAGGGATCCGGACGACCCGGATGCTTACCGTATGACGCCGGAGGAGGTGGGCGAGCGGGCAGCGGAAGCCCGTGACGCTGGAGCTACGGAAGTATGTCTTCAGGGCGGACTGCATCCCGAGGCGACGTTTGAGTACTACCTGGAAATGTTGGACGAGATCAAGTCCCAAGCCCCGGACATCCACGTGCACGGGTACTCACCGATGGAGGTGAAGTACTGCGCCAAGCTGGCGGGAGAGGACATCGAAGACGTACTACGAGAGCTGAAGCGAGCCGGTCTCGATTCGATGCCCGGAACGGCCGCGGAGATATTCTCCCCTGAGGTGAGGAAGCGGCTATGTCCTGATAAGTTGGAAGCCGATGAGTGGGAACATATCATCAGGATCGCGCACGAGTTGGGAATTCCCACCACTTGTACTATGATGTACGGTCACATCGACTCACCGAGGGACTGGATCGACCACATGAAGCGGCTTCGAGGGATCCAAGAGGACACGGGAGGCTTCACGGAGTTCGTGCCGCTCTCCTTCGTACATTCGAACGCACCGATTTACCGACGAGGAGGGGCGCGACCCGGAGTATCGGGTATGACGGACGTACTCGTGCACGCTGTGGCCCGATTGTACTTCGGACCGTTGATTCCGAACATACAGGCTTCCTGGGTGAAGCTCGGAGTGAAGCTGGCTCAGATGACGCTGCACGCCGGGGCGAACGATCTAGGTGGCACCCTCATGGAAGAGAACATCTCCCGGGAGGCCGGAGCGACCGAGGGCGAGCAGCTCGAGCCCGAGGAGATAGTGGAGATCATTCGGGAGGCGGGCTTCACCCCCGTGCAGCGCACCACGCTCTACGAGCCGGTGAAGGTGTACTAA</sequence>
    <xref id="ENA|AAM02110|AAM02110.1" name="ENA|AAM02110|AAM02110.1 Methanopyrus kandleri AV19 Predicted enzyme related to thiamine biosynthesis enzyme ThiH"/>
    <orf start="397" end="326" strand="ANTISENSE">
      <protein>
        <sequence md5="3B52A60935909B9B722B11EBF380610C">SRPTFPGSTQTSPRDAVRPEDILP</sequence>
        <xref id="orf1260" name="orf1260 source=ENA|AAM02110|AAM02110.1 coords=397..326 length=24 frame=6 desc=Methanopyrus kandleri AV19 Predicted enzyme related to thiamine biosynthesis enzyme ThiH"/>
        <matches>
          <match source="MobiDB-lite">
            <signature ac="mobidb-lite" name="disorder_prediction" type="Region">
              <signature-library-release library="MobiDB-lite" version="4.0"/>
            </signature>
            <model-ac>mobidb-lite</model-ac>
            <locations>
              <location start="1" end="24" representative="false">
                <location-fragments>
                  <fragment start="1" end="24" dc-status="CONTINUOUS"/>
                </location-fragments>
              </location>
              <location start="1" end="13" representative="false" sequence-feature="Polar">
                <location-fragments>
                  <fragment start="1" end="13" dc-status="CONTINUOUS"/>
                </location-fragments>
              </location>
              <location start="15" end="24" representative="false" sequence-feature="Polyampholyte">
                <location-fragments>
                  <fragment start="15" end="24" dc-status="CONTINUOUS"/>
                </location-fragments>
              </location>
            </locations>
          </match>
        </matches>
      </protein>
    </orf>
  </nucleotide-sequence>
</results>
```
