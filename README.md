# interproscan6

!! WORK IN PROGRESS !! 

Before to run you need to download:
- members bin (you can find in interproscan5 on path "interproscan/core/jms-implementation/target/interproscan-5-dist/bin")
- members data (curl ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.62-94.0/alt/interproscan-data-5.62-94.0.tar.gz --output interproscan-data-5.62-94.0.tar.gz)
- xref (entries, goterms and pathways) files (Use this script: interproscan6/files_test/get_data_to_i6.py)

IMPORTANT: Change the input params in input.yaml file if you want to test different flows (see in main.nf)

remember to build docker interproscan6 image (necessary to hmmer process):

    docker build -t interproscan6 .

How to run:

     nextflow run main.nf [-resume]

The results will apear in result/ folder.
You can debug seen work/ dir.

PS: batchsize parameter on nextflow.config define the quantity of threads of parallelism
