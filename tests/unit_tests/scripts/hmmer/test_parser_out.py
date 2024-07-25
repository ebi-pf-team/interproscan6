"""Test the python script scripts/hmmer/test_parser_out.py

These test are intend to be run from the root of the repository using:
pytest -v
"""

import os
import sys
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../../interproscan/scripts/hmmer')))
from interproscan.scripts.hmmer.parser_out import parse, get_accession_regex

hmmer_multiple_sequences = """
Query:       PALP  [M=295]
Accession:   PF00291.30
Description: Pyridoxal-phosphate dependent enzyme
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence       Description
    ------- ------ -----    ------- ------ -----   ---- --  --------       -----------
    9.3e-39  145.1   0.0    1.1e-38  144.8   0.0    1.0  1  OUTROSEMLOOKUP  
    5.5e-13   60.5   0.1    6.6e-13   60.2   0.1    1.0  1  TESTESEMLOOKUP  


Domain annotation for each sequence (and alignments):
>> OUTROSEMLOOKUP  
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  144.8   0.0   3.7e-46   1.1e-38       4     287 ..      12     315 ..       9     320 .. 0.91

  Alignments for each domain:
  == domain 1  score: 144.8 bits;  conditional E-value: 3.7e-46
                     CSSS--EEEECCTCCCTTCE..EEEEEGGGST...TSBTTHHHHHHHHHH.HHTTTTTSEEEEEBSSHHHHHHHHHHHHHT-EEEEEEETTS-.. CS
            PALP   4 gigpTPlvrlprlskelgve..vylKlEslnp...tgSfKdRgalnllar.lkegkggktvveassGNhGaalAaaaarlGlkvtivvpekas.. 90 
                       gp P+ +l+rls++lg++  +y+K+E++n    +g++K R++++l+ + +++g+++ + + + ++N  + +Aa+aa+lG+k+++v +++++  
  OUTROSEMLOOKUP  12 TFGPSPITPLKRLSQHLGGKveLYAKREDCNSglaFGGNKTRKLEYLIPEaIEQGCDTLVSIGGIQSNQTRQVAAVAAHLGMKCVLVQENWVNys 106
                     579****************************8799***************6777799999999999***************************88 PP

                     ......HHHHHHHHHTT-EEEEECCH....HHHHHHH.HHHHHHHSTTCEE--TTT..SHHHHHHHHTHHHHHH...HHHTTTEEEEEEE-SSSH CS
            PALP  91 ......peklaliralGaevvlvggd....ydeavel.akelaeegegayyinqyd..npaniegyktiglEil...eqlggkpdavvvpvGgGg 169
                            +++++ r +Ga+v+l  ++    ++ + e+ + +++e+g +++ i+ +   +p++  g+   + E+    ++lg k+d++vv+  +G+
  OUTROSEMLOOKUP 107 davydrVGNIEMSRIMGADVRLDAAGfdigIRPSWEKaMSDVVEQGGKPFPIPAGCseHPYGGLGFVGFAEEVRqqeKELGFKFDYIVVCSVTGS 201
                     7777766****************66623333333333355899999********9999***************97777777************** PP

                     HHHHHHHHHHHHSTTSEEEEEEECTTCGGTTCCS--S-SS--B-SSS-CCSTCSTCGTTCCHHHHSEEEEEEEEEEHHHHHHHHHHHHHHHSB-B CS
            PALP 170 liaGiarglkelgpevrvigvepegapalaksleagrpvkvksadtiadglgvgpepgelalelldeyvdevvtvsdeealeamrllarregilv 264
                     + aG++ g+++ g + +vig++++  p  +k+        +  a+++a+ ++ g e +e+++ l  +   + +++++e +lea+rl  + eg+l+
  OUTROSEMLOOKUP 202 TQAGMVVGFAADGRSKNVIGIDASAKPEQTKAQ------ILRIARHTAELVELGREITEEDVVLDTRFAYPEYGLPNEGTLEAIRLCGSLEGVLT 290
                     ***********88888*******9888888888......88889*************************************************** PP

                     -H.HHHHHHHHH.HHHHHHCCTTCE CS
            PALP 265 ep.ssaaalaal.klreagelkegd 287
                     +p + +++++++ +++++ge++eg+
  OUTROSEMLOOKUP 291 DPvYEGKSMHGMiEMVRRGEFPEGS 315
                     **********************994 PP

>> TESTESEMLOOKUP  
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   60.2   0.1   2.1e-20   6.6e-13      14     104 ..      12     116 ..      10     118 .. 0.87

  Alignments for each domain:
  == domain 1  score: 60.2 bits;  conditional E-value: 2.1e-20
                     CCTCCCTTCE..EEEEEGGGST...TSBTTHHHHHHHHHHHHTTTTTSEEEEEBS..SHHHHHHHHHHHHHT-EEEEEEETTS-........HHH CS
            PALP  14 prlskelgve..vylKlEslnp...tgSfKdRgalnllarlkegkggktvveass..GNhGaalAaaaarlGlkvtivvpekas........pek 93 
                     +rls++lg++  +y+K+E++n    +g++K R++++++ ++    g++t+v  ++  +Nh + +Aa+aa+ G+k+++v +++++         ++
  TESTESEMLOOKUP  12 SRLSAHLGGKvdLYAKREDCNSglaFGGNKLRKLEYIVPDAIAS-GADTLVSIGGvqSNHTRMVAAVAAKIGMKCRLVQEAWVPhedavydrVGN 105
                     5899******9**********8799***************7666.666666666555*********************999999666565555** PP

                     HHHHHHTT-EE CS
            PALP  94 laliralGaev 104
                     + l r +Ga+v
  TESTESEMLOOKUP 106 IMLSRIMGADV 116
                     *********98 PP
//
"""
expected_result_multiple_sequences = {
  "OUTROSEMLOOKUP": {
    "PF00291": {
      "accession": "PF00291",
      "name": "PALP",
      "description": "Pyridoxal-phosphate dependent enzyme",
      "evalue": 9.3e-39,
      "score": 145.1,
      "qlen": 295,
      "bias": "0.0",
      "member_db": "pfam",
      "version": "37.0",
      "model-ac": "PF00291",
      "locations": [
        {
          "start": 12,
          "end": 315,
          "representative": "",
          "hmmStart": 4,
          "hmmEnd": 287,
          "hmmLength": 295,
          "rawHmmBounds": "..",
          "hmmBounds": "INCOMPLETE",
          "evalue": 1.1e-38,
          "score": 144.8,
          "envelopeStart": "9",
          "envelopeEnd": "320",
          "alignment": "TFGPSPITPLKRLSQHLGGKveLYAKREDCNSglaFGGNKTRKLEYLIPEaIEQGCDTLVSIGGIQSNQTRQVAAVAAHLGMKCVLVQENWVNysdavydrVGNIEMSRIMGADVRLDAAGfdigIRPSWEKaMSDVVEQGGKPFPIPAGCseHPYGGLGFVGFAEEVRqqeKELGFKFDYIVVCSVTGSTQAGMVVGFAADGRSKNVIGIDASAKPEQTKAQ------ILRIARHTAELVELGREITEEDVVLDTRFAYPEYGLPNEGTLEAIRLCGSLEGVLTDPvYEGKSMHGMiEMVRRGEFPEGS",
          "cigar_alignment": "20M2I10M3I15M1I42M8I20M4I7M1I18M2I16M3I51M6D58M1I9M1I12M"
        }
      ]
    }
  },
  "TESTESEMLOOKUP": {
    "PF00291": {
      "accession": "PF00291",
      "name": "PALP",
      "description": "Pyridoxal-phosphate dependent enzyme",
      "evalue": 5.5e-13,
      "score": 60.5,
      "qlen": 295,
      "bias": "0.1",
      "member_db": "pfam",
      "version": "37.0",
      "model-ac": "PF00291",
      "locations": [
        {
          "start": 12,
          "end": 116,
          "representative": "",
          "hmmStart": 14,
          "hmmEnd": 104,
          "hmmLength": 295,
          "rawHmmBounds": "..",
          "hmmBounds": "INCOMPLETE",
          "evalue": 6.6e-13,
          "score": 60.2,
          "envelopeStart": "10",
          "envelopeEnd": "118",
          "alignment": "SRLSAHLGGKvdLYAKREDCNSglaFGGNKLRKLEYIVPDAIAS-GADTLVSIGGvqSNHTRMVAAVAAKIGMKCRLVQEAWVPhedavydrVGNIMLSRIMGADV",
          "cigar_alignment": "10M2I10M3I19M1D10M2I27M8I14M"
        }
      ]
    }
  }
}

hmmer_multiple_domains = """
Query:       3jt0B00-i2  [M=121]
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence       Description
    ------- ------ -----    ------- ------ -----   ---- --  --------       -----------
    0.00016   23.9   0.0     0.0005   22.3   0.0    1.8  2  WP_338726824.1  extracellular exonuclease ExeM [Shewanella ba


Domain annotation for each sequence (and alignments):
>> WP_338726824.1  extracellular exonuclease ExeM [Shewanella baltica]
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   22.3   0.0   7.7e-09    0.0005       4      88 ..      25     113 ..      22     123 .. 0.74
   2 ?   -3.3   0.0      0.65   4.2e+04      19      37 ..     776     796 ..     776     806 .. 0.77

  Alignments for each domain:
  == domain 1  score: 22.3 bits;  conditional E-value: 7.7e-09
                     xxxxxxx...xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx..xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx RF
      3jt0B00-i2   4 nveIaev...dpdgkfveLeNtsdkdvdLggWrlkrkvggqevs..ykfppsfvLkpgqtvtvwsagagathsppsdlvwksqdswgtgd 88 
                     +v I+e    ++++k +eL+N++d ++dL g++l+r ++g +v+  +   ++  ++p++t  + +++a  t ++  d +  s  s++ gd
  WP_338726824.1  25 DVMITEYvegSSNNKAIELYNSGDTAIDLAGYKLVRYKDGATVAsdMVALDGQSIAPKTTKVILNSSAVITLDQGVDSYSGSL-SFNGGD 113
                     56664443337899************************996655347778999999999999998888776666666665555.555444 PP

  == domain 2  score: -3.3 bits;  conditional E-value: 0.65
                     xxxxx..xxxxxxxxxxxxxx RF
      3jt0B00-i2  19 LeNts..dkdvdLggWrlkrk 37 
                     L+N s  ++ vd+++W + + 
  WP_338726824.1 776 LANASlaKRLVDIEDWHINSV 796
                     56777767789*****99765 PP

//
"""
expected_result_multiple_domains = {
    "WP_338726824.1": {
        "3jt0B00-i2": {
          "accession": "3jt0B00-i2",
          "name": "3jt0B00-i2",
          "description": "",
          "evalue": 0.00016,
          "score": 23.9,
          "qlen": 121,
          "bias": "0.0",
          "member_db": "gene3d",
          "version": "4.3.0",
          "model-ac": "3jt0B00-i2",
          "locations": [
            {
              "start": 22,
              "end": 123,
              "representative": "",
              "hmmStart": 4,
              "hmmEnd": 88,
              "hmmLength": 121,
              "rawHmmBounds": "..",
              "hmmBounds": "INCOMPLETE",
              "evalue": 0.0005,
              "score": 22.3,
              "envelopeStart": "22",
              "envelopeEnd": "123",
              "alignment": "DVMITEYvegSSNNKAIELYNSGDTAIDLAGYKLVRYKDGATVAsdMVALDGQSIAPKTTKVILNSSAVITLDQGVDSYSGSL-SFNGGD",
              "cigar_alignment": "7M3I34M2I37M1D6M"
            },
            {
              "start": 776,
              "end": 806,
              "representative": "",
              "hmmStart": 19,
              "hmmEnd": 37,
              "hmmLength": 121,
              "rawHmmBounds": "..",
              "hmmBounds": "INCOMPLETE",
              "evalue": 42000.0,
              "score": -3.3,
              "envelopeStart": "776",
              "envelopeEnd": "806",
              "alignment": "LANASlaKRLVDIEDWHINSV",
              "cigar_alignment": "5M2I14M"
            }
          ]
        },
        "3g91A00-i2": {
          "accession": "3g91A00-i2",
          "name": "3g91A00-i2",
          "description": "",
          "evalue": 6.1e-09,
          "score": 37.9,
          "qlen": 265,
          "bias": "0.0",
          "member_db": "gene3d",
          "version": "4.3.0",
          "model-ac": "3g91A00-i2",
          "locations": [
            {
              "start": 463,
              "end": 832,
              "representative": "",
              "hmmStart": 4,
              "hmmEnd": 256,
              "hmmLength": 265,
              "rawHmmBounds": "..",
              "hmmBounds": "INCOMPLETE",
              "evalue": 1.3e-08,
              "score": 36.8,
              "envelopeStart": "463",
              "envelopeEnd": "832",
              "alignment": "LRVASFNVlNFFNDVdggdtnpsgsnrgaltleemvlQRTKIVSAITAMNADIVGLMEIANNGfgeksaiKNLVDALNEKqtpENAYSFVEITDadkyDGKYFGtdaitVGMLYRggkvtlagaaqaidT----------PeqhasagsvtrtkdgktetnPGNDAYQRHSLAQTfkihdESLTVVVNHLKSKGSGClEdwanfeesvdpadqqgkCNAFRVSAAKVLGETL---KDVKGDLLIIGDMNaygmedpiRVLTDFDASKSDRDIMTASWTTLDGKVFERQGSKiekgyGLINLNTKAHGAgTYSYSYNGELGN---------LDHALANASLAKRLVDIEDwhinsvesnlfeygkkfsgdlaksENAFSASDHDPVIVALS",
              "cigar_alignment": "8M1I6M22I26M7I10M3I11M4I6M5I6M14I1M10D1M20I14M5I17M1I1M17I16M3D14M8I34M5I12M1I12M9D19M24I17M"
            }
          ]
        },
        "1akoA00-i2": {
          "accession": "1akoA00-i2",
          "name": "1akoA00-i2",
          "description": "",
          "evalue": 5.3e-07,
          "score": 31.6,
          "qlen": 268,
          "bias": "0.2",
          "member_db": "gene3d",
          "version": "4.3.0",
          "model-ac": "1akoA00-i2",
          "locations": [
            {
              "start": 464,
              "end": 831,
              "representative": "",
              "hmmStart": 1,
              "hmmEnd": 266,
              "hmmLength": 268,
              "rawHmmBounds": "[.",
              "hmmBounds": "N_TERMINAL_COMPLETE",
              "evalue": 1.3e-06,
              "score": 30.3,
              "envelopeStart": "464",
              "envelopeEnd": "831",
              "alignment": "LRVASFNVlNFFNDvdggdtnpsgsnrgaltleemvlQRTKIVSAITAMNADIVGLMEIANNGFGeksaiknlvdaLNEKQTPENAYSFVEiTdAdkyDGKYFGtdaitVGMLYRGgkvTLA-------GaaqaidtpeqhasagsvtrtkdgktetnPGNDAyQRHSLAQTFKIHDESLTVVVNHLKSKGSGCledwanfeesvdpadqQGKCNAFRVSAAKVLGETLK---DVKGDLLIIGDMNAYGMEDPIRVLTDFDASKSDRdiMTASWTTLDGKVFERQGskiekGYGLINLNTKAHGAGTYSYSYNGELGN---------LDHALANASLAKRLVDIEdwhINSVESnlfeygkkfsgdlakseNAFSASDHDPVIVAL",
              "cigar_alignment": "8M1I5M23I28M11I15M1I1M1I1M3I6M5I7M3I3M7D1M28I5M1I30M16I20M3D34M2I17M5I27M9D18M3I6M17I15M"
            }
          ]
        },
        "2j63A00-i2": {
          "accession": "2j63A00-i2",
          "name": "2j63A00-i2",
          "description": "",
          "evalue": 1e-05,
          "score": 27.1,
          "qlen": 358,
          "bias": "0.0",
          "member_db": "gene3d",
          "version": "4.3.0",
          "model-ac": "2j63A00-i2",
          "locations": [
            {
              "start": 445,
              "end": 831,
              "representative": "",
              "hmmStart": 33,
              "hmmEnd": 355,
              "hmmLength": 358,
              "rawHmmBounds": "..",
              "hmmBounds": "INCOMPLETE",
              "evalue": 1.9e-05,
              "score": 26.2,
              "envelopeStart": "445",
              "envelopeEnd": "831",
              "alignment": "PSVATKGDLRVASFNVLNFFNDVDGgdtnpsgsnrgaltleemvlQRTKIVSAITAMNADIVGLMEIAnngfgeksaiknlvdaLNEKQTPENAysfvEITDADKYDGKYFGTDaitvgmlyrgGKVTLAGAAQAIDTPEQHA-SAGSVTRTKDGKTETNPGNDAYQRhslAQTFKIHDESLTVVVNHLKSKG--SGCLEDWANFEESVDPADQQGKCNAFRV-SAAKVLGETLKDVKGD---LLIIGD---MNAYGMEDPIRVLTDFDASKSDRDIMTA-SWTTLDGKV-F---E-RQGSKIEKGYGLINLNTKAHGAGTYSYSYNGELGN---------LDHALANASLAKRLVDIedwhinsvesnlfeygkkfsgdlaKSENAFSASDHDPVIVAL",
              "cigar_alignment": "25M20I23M16I10M4I16M10I19M1D24M3I22M2D28M1D16M3D6M3D28M1D9M1D1M3D1M1D35M9D17M24I18M"
            }
          ]
        },
        "4fpvB00-i2": {
          "accession": "4fpvB00-i2",
          "name": "4fpvB00-i2",
          "description": "",
          "evalue": 1.4e-12,
          "score": 49.9,
          "qlen": 257,
          "bias": "0.0",
          "member_db": "gene3d",
          "version": "4.3.0",
          "model-ac": "4fpvB00-i2",
          "locations": [
            {
              "start": 458,
              "end": 830,
              "representative": "",
              "hmmStart": 9,
              "hmmEnd": 253,
              "hmmLength": 257,
              "rawHmmBounds": "..",
              "hmmBounds": "INCOMPLETE",
              "evalue": 2.3e-12,
              "score": 49.2,
              "envelopeStart": "458",
              "envelopeEnd": "830",
              "alignment": "GDLRVASFNVLnFFNDVDggdtnpsgsnrgaltleeMVLQRTKIVSAITAMNADIVGLMEIANNgFGEksAIKNLvdalnekqtpENAYSFVeitdadkydgkYFGTDAITVGMLYRGGKVTLAGAaQAIDTPEQhasagsvtrtkdgktetnpgndAYQRHSLAQTFKIHDESLTVVVNHLKSKgsgcledwanfeesvdpadqqGKCNAFRVSAAKVLGETLKDVKG--DLLIIGDMNAygmedpirvltdfdaS-KsdrdimtaswttldgkvfERQGsKIEKGYGLINLNTKAHGAG----TYSYSYN-------GELGNLDHALANASLAKRLVDIEDWHINSVESNlfeygkkfsgdlaksENAFSASDHDPVIV",
              "cigar_alignment": "11M1I6M18I28M1I3M2I5M10I7M11I23M1I8M22I28M21I23M2D10M15I1M1D1M18I4M1I19M4D7M7D33M15I14M"
            }
          ]
        }
      }
    }

hmmer_multipleseq_multipledom = """
    Query:       MF_00807  [M=338]
    Scores for complete sequences (score includes all domains):
       --- full sequence ---   --- best 1 domain ---    -#dom-
        E-value  score  bias    E-value  score  bias    exp  N  Sequence       Description
        ------- ------ -----    ------- ------ -----   ---- --  --------       -----------
       6.6e-220  716.4   0.7   1.3e-203  662.8   0.0    2.0  2  OUTROSEMLOOKUP
        9.5e-68  216.2   0.3    1.1e-67  215.9   0.3    1.1  1  TESTESEMLOOKUP


    Domain annotation for each sequence (and alignments):
    >> OUTROSEMLOOKUP
       #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
     ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
       1 !  662.8   0.0  8.8e-204  1.3e-203       1     316 [.       1     316 [.       1     317 [. 1.00
       2 !   54.4   0.2     1e-18   1.5e-18     302     338 .]     316     352 .]     315     352 .] 0.98

      Alignments for each domain:
      == domain 1  score: 662.8 bits;  conditional E-value: 8.8e-204
            MF_00807   1 mnlekferypltfgptpieklkrlsehlggkvelyakredcnsglafggnklrkleylipdaiasgadtlvsiggiqsnqtrqvaavaaklglkc 95
                         mnl++ferypltfgp+pi++lkrls+hlggkvelyakredcnsglafggnk+rkleylip+ai++g+dtlvsiggiqsnqtrqvaavaa+lg+kc
      OUTROSEMLOOKUP   1 MNLNRFERYPLTFGPSPITPLKRLSQHLGGKVELYAKREDCNSGLAFGGNKTRKLEYLIPEAIEQGCDTLVSIGGIQSNQTRQVAAVAAHLGMKC 95
                         89********************************************************************************************* PP

            MF_00807  96 vlvqenwvnyedavydrvgnielsrilgadvrlvedgfdigirkswekaledvkerggkpyaipagasehpygglgfvgfaeevrkqekelgfkf 190
                         vlvqenwvny+davydrvgnie+sri+gadvrl+++gfdigir+sweka++dv+e+ggkp++ipag+sehpygglgfvgfaeevr+qekelgfkf
      OUTROSEMLOOKUP  96 VLVQENWVNYSDAVYDRVGNIEMSRIMGADVRLDAAGFDIGIRPSWEKAMSDVVEQGGKPFPIPAGCSEHPYGGLGFVGFAEEVRQQEKELGFKF 190
                         *********************************************************************************************** PP

            MF_00807 191 dyivvcsvtgstqagmvvgfaadgrakrvigidasakpeqtkeqvlriakntaelvelgreiteddvvlderfaypeyglpneetleairlaarl 285
                         dyivvcsvtgstqagmvvgfaadgr+k+vigidasakpeqtk+q+lria++taelvelgreite+dvvld+rfaypeyglpne+tleairl+++l
      OUTROSEMLOOKUP 191 DYIVVCSVTGSTQAGMVVGFAADGRSKNVIGIDASAKPEQTKAQILRIARHTAELVELGREITEEDVVLDTRFAYPEYGLPNEGTLEAIRLCGSL 285
                         *********************************************************************************************** PP

            MF_00807 286 egvltdpvyegksmqglielvrkgefpegsk 316
                         egvltdpvyegksm+g+ie+vr+gefpegs
      OUTROSEMLOOKUP 286 EGVLTDPVYEGKSMHGMIEMVRRGEFPEGSM 316
                         *****************************96 PP

      == domain 2  score: 54.4 bits;  conditional E-value: 1e-18
            MF_00807 302 lielvrkgefpegskvlyahlggapalnaysflfrdg 338
                         +i+lvrkg+fp+g+kvlyahlggapaln+y+++fr+g
      OUTROSEMLOOKUP 316 MIDLVRKGFFPAGAKVLYAHLGGAPALNGYAYAFRNG 352
                         9***********************************8 PP

    >> TESTESEMLOOKUP
       #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
     ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
       1 !  215.9   0.3   7.7e-68   1.1e-67      21     126 ..      11     116 ..       2     119 .] 0.97

      Alignments for each domain:
      == domain 1  score: 215.9 bits;  conditional E-value: 7.7e-68
            MF_00807  21 lkrlsehlggkvelyakredcnsglafggnklrkleylipdaiasgadtlvsiggiqsnqtrqvaavaaklglkcvlvqenwvnyedavydrvgn 115
                          +rls+hlggkv+lyakredcnsglafggnklrkley++pdaiasgadtlvsigg+qsn+tr+vaavaak+g+kc+lvqe+wv++edavydrvgn
      TESTESEMLOOKUP  11 QSRLSAHLGGKVDLYAKREDCNSGLAFGGNKLRKLEYIVPDAIASGADTLVSIGGVQSNHTRMVAAVAAKIGMKCRLVQEAWVPHEDAVYDRVGN 105
                         57********************************************************************************************* PP

            MF_00807 116 ielsrilgadv 126
                         i+lsri+gadv
      TESTESEMLOOKUP 106 IMLSRIMGADV 116
                         **********8 PP



    Internal pipeline statistics summary:
    -------------------------------------
    Query model(s):                            1  (338 nodes)
    Target sequences:                          3  (1341 residues searched)
    Passed MSV filter:                         2  (0.666667); expected 0.1 (0.02)
    Passed bias filter:                        2  (0.666667); expected 0.1 (0.02)
    Passed Vit filter:                         2  (0.666667); expected 0.0 (0.001)
    Passed Fwd filter:                         2  (0.666667); expected 0.0 (1e-05)
    Initial search space (Z):                  3  [actual number of targets]
    Domain search space  (domZ):               2  [number of targets reported over threshold]
    # CPU time: 0.00u 0.00s 00:00:00.00 Elapsed: 00:00:00.01
    # Mc/sec: 25.49
    //
    """
expected_result_multipleseq_multipledom = {
  "OUTROSEMLOOKUP": {
    "": {
      "accession": "",
      "name": "MF_00807",
      "description": "",
      "evalue": 6.6e-220,
      "score": 716.4,
      "qlen": 338,
      "bias": "0.7",
      "member_db": "panther",
      "version": "18.0",
      "model-ac": "",
      "locations": [
        {
          "start": 1,
          "end": 316,
          "representative": "",
          "hmmStart": 1,
          "hmmEnd": 316,
          "hmmLength": 338,
          "rawHmmBounds": "[.",
          "hmmBounds": "N_TERMINAL_COMPLETE",
          "evalue": 1.3e-203,
          "score": 662.8,
          "envelopeStart": "1",
          "envelopeEnd": "317",
          "alignment": "MNLNRFERYPLTFGPSPITPLKRLSQHLGGKVELYAKREDCNSGLAFGGNKTRKLEYLIPEAIEQGCDTLVSIGGIQSNQTRQVAAVAAHLGMKCVLVQENWVNYSDAVYDRVGNIEMSRIMGADVRLDAAGFDIGIRPSWEKAMSDVVEQGGKPFPIPAGCSEHPYGGLGFVGFAEEVRQQEKELGFKFDYIVVCSVTGSTQAGMVVGFAADGRSKNVIGIDASAKPEQTKAQILRIARHTAELVELGREITEEDVVLDTRFAYPEYGLPNEGTLEAIRLCGSLEGVLTDPVYEGKSMHGMIEMVRRGEFPEGSM",
          "cigar_alignment": "316M"
        },
        {
          "start": 316,
          "end": 352,
          "representative": "",
          "hmmStart": 302,
          "hmmEnd": 338,
          "hmmLength": 338,
          "rawHmmBounds": ".]",
          "hmmBounds": "C_TERMINAL_COMPLETE",
          "evalue": 1.5e-18,
          "score": 54.4,
          "envelopeStart": "315",
          "envelopeEnd": "352",
          "alignment": "MIDLVRKGFFPAGAKVLYAHLGGAPALNGYAYAFRNG",
          "cigar_alignment": "37M"
        }
      ]
    }
  },
  "TESTESEMLOOKUP": {
    "": {
      "accession": "",
      "name": "MF_00807",
      "description": "",
      "evalue": 9.5e-68,
      "score": 216.2,
      "qlen": 338,
      "bias": "0.3",
      "member_db": "panther",
      "version": "18.0",
      "model-ac": "",
      "locations": [
        {
          "start": 11,
          "end": 116,
          "representative": "",
          "hmmStart": 21,
          "hmmEnd": 126,
          "hmmLength": 338,
          "rawHmmBounds": "..",
          "hmmBounds": "INCOMPLETE",
          "evalue": 1.1e-67,
          "score": 215.9,
          "envelopeStart": "2",
          "envelopeEnd": "119",
          "alignment": "QSRLSAHLGGKVDLYAKREDCNSGLAFGGNKLRKLEYIVPDAIASGADTLVSIGGVQSNHTRMVAAVAAKIGMKCRLVQEAWVPHEDAVYDRVGNIMLSRIMGADV",
          "cigar_alignment": "106M"
        }
      ]
    }
  }
}


parsed_result_multiple_sequences = ""
file_path = os.path.abspath("0.1._.membername1._.out")
with open(file_path2, 'w') as f:
    f.write(hmmer_multiple_sequences)
    parsed_result_multiple_sequences = parse(file_path)
    os.remove(file_path)

parsed_result_multiple_domains = ""
file_path2 = os.path.abspath("0.2._.membername2._.out")
with open(file_path2, 'w') as f:
    f.write(hmmer_multiple_domains)
    parsed_result_multiple_domains = parse(file_path2)
    os.remove(file_path2)

parsed_result_multipleseq_multipledom = ""
file_path2 = os.path.abspath("0.3._.membername3._.out")
with open(file_path3, 'w') as f:
    f.write(hmmer_multiple_domains)
    parsed_result_multiple_domains = parse(file_path3)
    os.remove(file_path3)


@pytest.mark.parametrize((parsed_result, expected_result))
def test_match_multiple_domains_and_sequences(input_text, expected_output):
    assert parsed_output == expected_output, f"Expected {expected_output} but got {parsed_output}"


@pytest.mark.parametrize((parsed_result))
def test_result_types(result):
    # Make sure about correct types is extremely important for some members filtering
    assert isinstance(parsed_result["OUTROSEMLOOKUP"][""]["evalue"], float)
    assert isinstance(parsed_result["OUTROSEMLOOKUP"][""]["score"], float)
    assert isinstance(parsed_result["OUTROSEMLOOKUP"][""]["qlen"], int)
    assert isinstance(parsed_result["OUTROSEMLOOKUP"][""]["bias"], float)

    location = parsed_result["OUTROSEMLOOKUP"][""]["locations"][0]
    assert isinstance(location["start"], int)
    assert isinstance(location["end"], int)
    assert isinstance(location["hmmStart"], int)
    assert isinstance(location["hmmEnd"], int)
    assert isinstance(location["hmmLength"], int)
    assert isinstance(location["evalue"], float)
    assert isinstance(location["score"], float)
    assert isinstance(location["envelopeStart"], int)
    assert isinstance(location["envelopeEnd"], int)


# @pytest.mark.parametrize("member_db,query_line,expected_query_name,expected_qlen", [
#     ("PFAM", "Query: exonuc_ExeM-GG [M=884]", "exonuc_ExeM-GG", "884"),
#     ("PANTHER", "Query: PTHR31901.orig.30.pir [M=558]", "PTHR31901.orig.30.pir", "558"),
# ])
# def test_parse_query(member_db, query_line, expected_query_name, expected_qlen):
#     member_accession = get_accession_regex(member_db)
#     match = member_accession.match(query_line)
#     assert match is not None, f"Failed to match query line for {member_db}"
#     query_name = match.group(2)
#     qlen = query_line.split("[M=")[1].split("]")[0]
#     assert query_name == expected_query_name, f"Expected {expected_query_name} but got {query_name}"
#     assert qlen == expected_qlen, f"Expected {expected_qlen} but got {qlen}"


# @pytest.mark.parametrize("appl,query_line,expected_accession", [
#     ("NCBIFAM", "Accession: NF033680.1", "NF033680.1"),
#     ("ANTIFAM", "Accession: ANF00001", "ANF00001"),
# ])
# def test_parse_accession():
#     member_accession = get_accession_regex(appl)
#     match = member_accession.match(query_line)
#     assert match is not None, f"Failed to match accession line for {appl}"
#     accession = match.group(2)
#     assert accession == expected_accession, f"Expected {expected_accession} but got {accession}"

# def test_parse_filename():
