"""Test the python script scripts/hmmer/test_parser_out.py

These test are intend to be run from the root of the repository using:
pytest -v
"""

import os
import sys
import pytest
import tempfile
import textwrap

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../../interproscan/scripts/hmmer')))
from interproscan.scripts.hmmer.parser_out import parse, get_accession_regex

HMMER_OUT_MULTIPLE_SEQUENCES = textwrap.dedent("""
Query:       PALP  [M=295]
Accession:   PF00291.30
Description: Pyridoxal-phosphate dependent enzyme
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence                  Description
    ------- ------ -----    ------- ------ -----   ---- --  --------                  -----------
    8.2e-43  158.4   0.3    9.7e-43  158.1   0.3    1.0  1  Protein_with_panther_hits  
    8.2e-43  158.4   0.3    9.7e-43  158.1   0.3    1.0  1  sp|A2SLW2|1A1D_METPP       1-aminocyclopropane-1-carboxylate 
    9.3e-39  145.1   0.0    1.1e-38  144.8   0.0    1.0  1  OUTROSEMLOOKUP             
    5.5e-13   60.5   0.1    6.6e-13   60.2   0.1    1.0  1  TESTESEMLOOKUP             


Domain annotation for each sequence (and alignments):
>> Protein_with_panther_hits  
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  158.1   0.3   6.4e-50   9.7e-43       4     293 ..      12     321 ..      10     323 .. 0.90

  Alignments for each domain:
  == domain 1  score: 158.1 bits;  conditional E-value: 6.4e-50
                                CSSS--EEEECCTCCCTTCE..EEEEEGGGST...TSBTTHHHHHHHHHH.HHTTTTTSEEEEEBSSHHHHHHHHHHHHHT-EE CS
                       PALP   4 gigpTPlvrlprlskelgve..vylKlEslnp...tgSfKdRgalnllar.lkegkggktvveassGNhGaalAaaaarlGlkv 81 
                                  gpTP+++l+rls++lg++  +y+K+E++n    +g++K R++++l+ + l++g+++ + + + ++N  + +Aa+aa+lG+k+
  Protein_with_panther_hits  12 TFGPTPIQPLKRLSAHLGGQveLYAKREDCNSglaFGGNKTRKLEYLIPEaLAQGCDTLVSIGGIQSNQTRQVAAVAAHLGMKC 95 
                                579****************************8799***************99999***************************** PP

                                EEEEETTS-........HHHHHHHHHTT-EEEEECCHHHHHHHH.....HHHHHHHSTTCEE--TTT..SHHHHHHHHTHHHHH CS
                       PALP  82 tivvpekas........peklaliralGaevvlvggdydeavel.....akelaeegegayyinqyd..npaniegyktiglEi 150
                                ++v +++++         +++++ r lGa+v+l  +++d   +      ++++ + g +++ i+ +   + ++  g+   + E+
  Protein_with_panther_hits  96 VLVQENWVNysdavydrVGNIEMSRILGADVRLDAAGFDIGIRPsweqaMADVRAAGGKPFPIPAGCseHRLGGLGFVGFAEEV 179
                                *********887777766****************666443333322111445666689999999999556668888999999** PP

                                H...HHHTTTEEEEEEE-SSSHHHHHHHHHHHHHSTTSEEEEEEECTTCGGTTCCS--S-SS--B-SSS-CCSTCSTCGTTCCH CS
                       PALP 151 l...eqlggkpdavvvpvGgGgliaGiarglkelgpevrvigvepegapalaksleagrpvkvksadtiadglgvgpepgelal 231
                                    ++lg k+d++vv+  +G++ aG++ g+++ g + rvig++++  p  +++        +  a+ +a+ +g g e +++++
  Protein_with_panther_hits 180 RaqeAELGFKFDYIVVCSVTGSTQAGMVVGFAADGRAERVIGIDASAKPEQTHAQ------ILRIAQNTAELVGLGREITAQDV 257
                                97888888*************************99999********888888888......88888****************** PP

                                HHHSEEEEEEEEEEHHHHHHHHHHHHHHHSB-B-H.HHHHHHHHH.HHHHHHCCTTCEEEEEEE CS
                       PALP 232 elldeyvdevvtvsdeealeamrllarregilvep.ssaaalaal.klreagelkegdrvvvvl 293
                                 l  +y+++ +++++e +lea+rl ar eg+l++p + +++++++ + +++ge+++g+rv+  +
  Protein_with_panther_hits 258 VLDTRYGGPEYGLPSEGTLEAIRLCARQEGMLTDPvYEGKSMHGMiDKVKRGEFPAGSRVLYAH 321
                                ************************************************************9876 PP

>> sp|A2SLW2|1A1D_METPP  1-aminocyclopropane-1-carboxylate deaminase OS=Methylibium petroleiphilum (strain ATCC BAA-1232
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  158.1   0.3   6.4e-50   9.7e-43       4     293 ..      12     321 ..      10     323 .. 0.90

  Alignments for each domain:
  == domain 1  score: 158.1 bits;  conditional E-value: 6.4e-50
                           CSSS--EEEECCTCCCTTCE..EEEEEGGGST...TSBTTHHHHHHHHHH.HHTTTTTSEEEEEBSSHHHHHHHHHHHHHT-EEEEEEE CS
                  PALP   4 gigpTPlvrlprlskelgve..vylKlEslnp...tgSfKdRgalnllar.lkegkggktvveassGNhGaalAaaaarlGlkvtivvp 86 
                             gpTP+++l+rls++lg++  +y+K+E++n    +g++K R++++l+ + l++g+++ + + + ++N  + +Aa+aa+lG+k+++v +
  sp|A2SLW2|1A1D_METPP  12 TFGPTPIQPLKRLSAHLGGQveLYAKREDCNSglaFGGNKTRKLEYLIPEaLAQGCDTLVSIGGIQSNQTRQVAAVAAHLGMKCVLVQE 100
                           579****************************8799***************99999********************************** PP

                           TTS-........HHHHHHHHHTT-EEEEECCHHHHHHHH.....HHHHHHHSTTCEE--TTT..SHHHHHHHHTHHHHHH...HHHTTT CS
                  PALP  87 ekas........peklaliralGaevvlvggdydeavel.....akelaeegegayyinqyd..npaniegyktiglEil...eqlggk 157
                           ++++         +++++ r lGa+v+l  +++d   +      ++++ + g +++ i+ +   + ++  g+   + E+    ++lg k
  sp|A2SLW2|1A1D_METPP 101 NWVNysdavydrVGNIEMSRILGADVRLDAAGFDIGIRPsweqaMADVRAAGGKPFPIPAGCseHRLGGLGFVGFAEEVRaqeAELGFK 189
                           ****887777766****************666443333322111445666689999999999556668888999999**97888888** PP

                           EEEEEEE-SSSHHHHHHHHHHHHHSTTSEEEEEEECTTCGGTTCCS--S-SS--B-SSS-CCSTCSTCGTTCCHHHHSEEEEEEEEEEH CS
                  PALP 158 pdavvvpvGgGgliaGiarglkelgpevrvigvepegapalaksleagrpvkvksadtiadglgvgpepgelalelldeyvdevvtvsd 246
                           +d++vv+  +G++ aG++ g+++ g + rvig++++  p  +++        +  a+ +a+ +g g e +++++ l  +y+++ +++++
  sp|A2SLW2|1A1D_METPP 190 FDYIVVCSVTGSTQAGMVVGFAADGRAERVIGIDASAKPEQTHAQ------ILRIAQNTAELVGLGREITAQDVVLDTRYGGPEYGLPS 272
                           ***********************99999********888888888......88888********************************* PP

                           HHHHHHHHHHHHHHSB-B-H.HHHHHHHHH.HHHHHHCCTTCEEEEEEE CS
                  PALP 247 eealeamrllarregilvep.ssaaalaal.klreagelkegdrvvvvl 293
                           e +lea+rl ar eg+l++p + +++++++ + +++ge+++g+rv+  +
  sp|A2SLW2|1A1D_METPP 273 EGTLEAIRLCARQEGMLTDPvYEGKSMHGMiDKVKRGEFPAGSRVLYAH 321
                           *********************************************9876 PP

>> OUTROSEMLOOKUP  
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  144.8   0.0   7.5e-46   1.1e-38       4     287 ..      12     315 ..       9     320 .. 0.91

  Alignments for each domain:
  == domain 1  score: 144.8 bits;  conditional E-value: 7.5e-46
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
   1 !   60.2   0.1   4.3e-20   6.6e-13      14     104 ..      12     116 ..      10     118 .. 0.87

  Alignments for each domain:
  == domain 1  score: 60.2 bits;  conditional E-value: 4.3e-20
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



Internal pipeline statistics summary:
-------------------------------------
Query model(s):                            1  (295 nodes)
Target sequences:                         14  (5487 residues searched)
Passed MSV filter:                         5  (0.357143); expected 0.3 (0.02)
Passed bias filter:                        4  (0.285714); expected 0.3 (0.02)
Passed Vit filter:                         4  (0.285714); expected 0.0 (0.001)
Passed Fwd filter:                         4  (0.285714); expected 0.0 (1e-05)
Initial search space (Z):           61295632  [as set by --Z on cmdline]
Domain search space  (domZ):               4  [number of targets reported over threshold]
# CPU time: 0.00u 0.00s 00:00:00.00 Elapsed: 00:00:00.00
# Mc/sec: 235.36
//
""")
EXPECTED_RESULT_MULTIPLE_SEQUENCES = {
   "Protein_with_panther_hits":{
      "PF00291":{
         "accession":"PF00291",
         "name":"PALP",
         "description":"Pyridoxal-phosphate dependent enzyme",
         "evalue":8.2e-43,
         "score":158.4,
         "qlen":295,
         "bias":"0.3",
         "member_db":"pfam",
         "version":"0.0",
         "model-ac":"PF00291",
         "locations":[
            {
               "start":12,
               "end":321,
               "representative":"",
               "hmmStart":4,
               "hmmEnd":293,
               "hmmLength":295,
               "rawHmmBounds":"..",
               "hmmBounds":"INCOMPLETE",
               "evalue":9.7e-43,
               "score":158.1,
               "envelopeStart":"10",
               "envelopeEnd":"323",
               "alignment":"TFGPTPIQPLKRLSAHLGGQveLYAKREDCNSglaFGGNKTRKLEYLIPEaLAQGCDTLVSIGGIQSNQTRQVAAVAAHLGMKCVLVQENWVNysdavydrVGNIEMSRILGADVRLDAAGFDIGIRPsweqaMADVRAAGGKPFPIPAGCseHRLGGLGFVGFAEEVRaqeAELGFKFDYIVVCSVTGSTQAGMVVGFAADGRAERVIGIDASAKPEQTHAQ------ILRIAQNTAELVGLGREITAQDVVLDTRYGGPEYGLPSEGTLEAIRLCARQEGMLTDPvYEGKSMHGMiDKVKRGEFPAGSRVLYAH",
               "cigar_alignment":"20M2I10M3I15M1I42M8I27M5I18M2I16M3I51M6D58M1I9M1I18M"
            }
         ]
      }
   },
   "sp|A2SLW2|1A1D_METPP":{
      "PF00291":{
         "accession":"PF00291",
         "name":"PALP",
         "description":"Pyridoxal-phosphate dependent enzyme",
         "evalue":8.2e-43,
         "score":158.4,
         "qlen":295,
         "bias":"0.3",
         "member_db":"pfam",
         "version":"0.0",
         "model-ac":"PF00291",
         "locations":[
            {
               "start":12,
               "end":321,
               "representative":"",
               "hmmStart":4,
               "hmmEnd":293,
               "hmmLength":295,
               "rawHmmBounds":"..",
               "hmmBounds":"INCOMPLETE",
               "evalue":9.7e-43,
               "score":158.1,
               "envelopeStart":"10",
               "envelopeEnd":"323",
               "alignment":"TFGPTPIQPLKRLSAHLGGQveLYAKREDCNSglaFGGNKTRKLEYLIPEaLAQGCDTLVSIGGIQSNQTRQVAAVAAHLGMKCVLVQENWVNysdavydrVGNIEMSRILGADVRLDAAGFDIGIRPsweqaMADVRAAGGKPFPIPAGCseHRLGGLGFVGFAEEVRaqeAELGFKFDYIVVCSVTGSTQAGMVVGFAADGRAERVIGIDASAKPEQTHAQ------ILRIAQNTAELVGLGREITAQDVVLDTRYGGPEYGLPSEGTLEAIRLCARQEGMLTDPvYEGKSMHGMiDKVKRGEFPAGSRVLYAH",
               "cigar_alignment":"20M2I10M3I15M1I42M8I27M5I18M2I16M3I51M6D58M1I9M1I18M"
            }
         ]
      }
   },
   "OUTROSEMLOOKUP":{
      "PF00291":{
         "accession":"PF00291",
         "name":"PALP",
         "description":"Pyridoxal-phosphate dependent enzyme",
         "evalue":9.3e-39,
         "score":145.1,
         "qlen":295,
         "bias":"0.0",
         "member_db":"pfam",
         "version":"0.0",
         "model-ac":"PF00291",
         "locations":[
            {
               "start":12,
               "end":315,
               "representative":"",
               "hmmStart":4,
               "hmmEnd":287,
               "hmmLength":295,
               "rawHmmBounds":"..",
               "hmmBounds":"INCOMPLETE",
               "evalue":1.1e-38,
               "score":144.8,
               "envelopeStart":"9",
               "envelopeEnd":"320",
               "alignment":"TFGPSPITPLKRLSQHLGGKveLYAKREDCNSglaFGGNKTRKLEYLIPEaIEQGCDTLVSIGGIQSNQTRQVAAVAAHLGMKCVLVQENWVNysdavydrVGNIEMSRIMGADVRLDAAGfdigIRPSWEKaMSDVVEQGGKPFPIPAGCseHPYGGLGFVGFAEEVRqqeKELGFKFDYIVVCSVTGSTQAGMVVGFAADGRSKNVIGIDASAKPEQTKAQ------ILRIARHTAELVELGREITEEDVVLDTRFAYPEYGLPNEGTLEAIRLCGSLEGVLTDPvYEGKSMHGMiEMVRRGEFPEGS",
               "cigar_alignment":"20M2I10M3I15M1I42M8I20M4I7M1I18M2I16M3I51M6D58M1I9M1I12M"
            }
         ]
      }
   },
   "TESTESEMLOOKUP":{
      "PF00291":{
         "accession":"PF00291",
         "name":"PALP",
         "description":"Pyridoxal-phosphate dependent enzyme",
         "evalue":5.5e-13,
         "score":60.5,
         "qlen":295,
         "bias":"0.1",
         "member_db":"pfam",
         "version":"0.0",
         "model-ac":"PF00291",
         "locations":[
            {
               "start":12,
               "end":116,
               "representative":"",
               "hmmStart":14,
               "hmmEnd":104,
               "hmmLength":295,
               "rawHmmBounds":"..",
               "hmmBounds":"INCOMPLETE",
               "evalue":6.6e-13,
               "score":60.2,
               "envelopeStart":"10",
               "envelopeEnd":"118",
               "alignment":"SRLSAHLGGKvdLYAKREDCNSglaFGGNKLRKLEYIVPDAIAS-GADTLVSIGGvqSNHTRMVAAVAAKIGMKCRLVQEAWVPhedavydrVGNIMLSRIMGADV",
               "cigar_alignment":"10M2I10M3I19M1D10M2I27M8I14M"
            }
         ]
      }
   }
}

HMMER_OUT_MULTIPLE_DOMAINS = textwrap.dedent("""
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
  """)
EXPECTED_RESULT_MULTIPLE_DOMAINS = {
   "WP_338726824.1":{
      "3jt0B00-i2":{
         "accession":"3jt0B00-i2",
         "name":"3jt0B00-i2",
         "description":"",
         "evalue":0.00016,
         "score":23.9,
         "qlen":121,
         "bias":"0.0",
         "member_db":"gene3d",
         "version":"0.0",
         "model-ac":"3jt0B00-i2",
         "locations":[
            {
               "start":22,
               "end":123,
               "representative":"",
               "hmmStart":4,
               "hmmEnd":88,
               "hmmLength":121,
               "rawHmmBounds":"..",
               "hmmBounds":"INCOMPLETE",
               "evalue":0.0005,
               "score":22.3,
               "envelopeStart":"22",
               "envelopeEnd":"123",
               "alignment":"DVMITEYvegSSNNKAIELYNSGDTAIDLAGYKLVRYKDGATVAsdMVALDGQSIAPKTTKVILNSSAVITLDQGVDSYSGSL-SFNGGD",
               "cigar_alignment":"7M3I34M2I37M1D6M"
            },
            {
               "start":776,
               "end":806,
               "representative":"",
               "hmmStart":19,
               "hmmEnd":37,
               "hmmLength":121,
               "rawHmmBounds":"..",
               "hmmBounds":"INCOMPLETE",
               "evalue":42000.0,
               "score":-3.3,
               "envelopeStart":"776",
               "envelopeEnd":"806",
               "alignment":"LANASlaKRLVDIEDWHINSV",
               "cigar_alignment":"5M2I14M"
            }
         ]
      }
   }
}

HMMER_OUT_MULTIPLESEQ_MULTIPLEDOM = textwrap.dedent("""
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
      """)
EXPECTED_RESULT_MULTIPLESEQ_MULTIPLEDOM = {
   "OUTROSEMLOOKUP":{
      "":{
         "accession":"",
         "name":"MF_00807",
         "description":"",
         "evalue":6.6e-220,
         "score":716.4,
         "qlen":338,
         "bias":"0.7",
         "member_db":"panther",
         "version":"0.0",
         "model-ac":"",
         "locations":[
            {
               "start":1,
               "end":316,
               "representative":"",
               "hmmStart":1,
               "hmmEnd":316,
               "hmmLength":338,
               "rawHmmBounds":"[.",
               "hmmBounds":"N_TERMINAL_COMPLETE",
               "evalue":1.3e-203,
               "score":662.8,
               "envelopeStart":"1",
               "envelopeEnd":"317",
               "alignment":"MNLNRFERYPLTFGPSPITPLKRLSQHLGGKVELYAKREDCNSGLAFGGNKTRKLEYLIPEAIEQGCDTLVSIGGIQSNQTRQVAAVAAHLGMKCVLVQENWVNYSDAVYDRVGNIEMSRIMGADVRLDAAGFDIGIRPSWEKAMSDVVEQGGKPFPIPAGCSEHPYGGLGFVGFAEEVRQQEKELGFKFDYIVVCSVTGSTQAGMVVGFAADGRSKNVIGIDASAKPEQTKAQILRIARHTAELVELGREITEEDVVLDTRFAYPEYGLPNEGTLEAIRLCGSLEGVLTDPVYEGKSMHGMIEMVRRGEFPEGSM",
               "cigar_alignment":"316M"
            },
            {
               "start":316,
               "end":352,
               "representative":"",
               "hmmStart":302,
               "hmmEnd":338,
               "hmmLength":338,
               "rawHmmBounds":".]",
               "hmmBounds":"C_TERMINAL_COMPLETE",
               "evalue":1.5e-18,
               "score":54.4,
               "envelopeStart":"315",
               "envelopeEnd":"352",
               "alignment":"MIDLVRKGFFPAGAKVLYAHLGGAPALNGYAYAFRNG",
               "cigar_alignment":"37M"
            }
         ]
      }
   },
   "TESTESEMLOOKUP":{
      "":{
         "accession":"",
         "name":"MF_00807",
         "description":"",
         "evalue":9.5e-68,
         "score":216.2,
         "qlen":338,
         "bias":"0.3",
         "member_db":"panther",
         "version":"0.0",
         "model-ac":"",
         "locations":[
            {
               "start":11,
               "end":116,
               "representative":"",
               "hmmStart":21,
               "hmmEnd":126,
               "hmmLength":338,
               "rawHmmBounds":"..",
               "hmmBounds":"INCOMPLETE",
               "evalue":1.1e-67,
               "score":215.9,
               "envelopeStart":"2",
               "envelopeEnd":"119",
               "alignment":"QSRLSAHLGGKVDLYAKREDCNSGLAFGGNKLRKLEYIVPDAIASGADTLVSIGGVQSNHTRMVAAVAAKIGMKCRLVQEAWVPHEDAVYDRVGNIMLSRIMGADV",
               "cigar_alignment":"106M"
            }
         ]
      }
   }
}


@pytest.fixture
def mock_hmmer_out_multiple_sequences(tmpdir):
    input_file_path = os.path.join(tmpdir, "0.0._.pfam._.out")
    with open(input_file_path, 'w') as f:
        f.write(HMMER_OUT_MULTIPLE_SEQUENCES)
    return input_file_path


@pytest.fixture
def mock_hmmer_out_multiple_domains(tmpdir):
    input_file_path = os.path.join(tmpdir, "0.0._.gene3d._.out")
    with open(input_file_path, 'w') as f:
        f.write(HMMER_OUT_MULTIPLE_DOMAINS)
    return input_file_path


@pytest.fixture
def mock_hmmer_out_multiple_sequences_multiple_domains(tmpdir):
    input_file_path = os.path.join(tmpdir, "0.0._.panther._.out")
    with open(input_file_path, 'w') as f:
        f.write(HMMER_OUT_MULTIPLESEQ_MULTIPLEDOM)
    return input_file_path


def test_check_output_types(mock_hmmer_out_multiple_sequences):
    result = parse(mock_hmmer_out_multiple_sequences)

    sequence_data = result["sp|A2SLW2|1A1D_METPP"]["PF00291"]

    # Make sure about correct types is extremely important for some members filtering
    assert isinstance(sequence_data["evalue"], float)
    assert isinstance(sequence_data["score"], float)
    assert isinstance(sequence_data["qlen"], int)

    location = sequence_data["locations"][0]
    assert isinstance(location["start"], int)
    assert isinstance(location["end"], int)
    assert isinstance(location["hmmStart"], int)
    assert isinstance(location["hmmEnd"], int)
    assert isinstance(location["hmmLength"], int)
    assert isinstance(location["evalue"], float)
    assert isinstance(location["score"], float)


@pytest.mark.parametrize("expected_result", [EXPECTED_RESULT_MULTIPLE_SEQUENCES])
def test_parse_multiple_sequences(mock_hmmer_out_multiple_sequences, expected_result):
    result = parse(mock_hmmer_out_multiple_sequences)
    assert result == expected_result


@pytest.mark.parametrize("expected_result", [EXPECTED_RESULT_MULTIPLE_DOMAINS])
def test_parse_multiple_domains(mock_hmmer_out_multiple_domains, expected_result):
    result = parse(mock_hmmer_out_multiple_domains)
    assert result == expected_result


@pytest.mark.parametrize("expected_result", [EXPECTED_RESULT_MULTIPLESEQ_MULTIPLEDOM])
def test_parse_multiple_sequences_multiple_domains(mock_hmmer_out_multiple_sequences_multiple_domains, expected_result):
    result = parse(mock_hmmer_out_multiple_sequences_multiple_domains)
    assert result == expected_result
