# interproscan6


!! WORK IN PROGRESS !! 
PS: The code still need to be refactored (the focus now is in obtain the correct results for all members in all flows)

Before run you need to upload:
- members bin 
- members data
- xref files (entries, goterms and pathways) 

Change the input params in input_opt.yaml file if want to test different flows (see in main.nf)

How to run (example):

     nextflow run main.nf --input files_test/test_all_appl.fasta -params-file input_opt.yaml

Temporary (scripts that still not in nextflow flow to make easy the test):
    CDD:
        If want to run cdd.py build docker and run:
        
            docker build -t interproscan6 .
            docker run -v ./results:/opt/interproscan6/results interproscan6