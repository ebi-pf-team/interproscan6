"""Test the python script scripts/hmmer/test_parser_out.py

These test are intend to be run from the root of the repository using:
pytest -v
"""

import os
import sys
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../interproscan/scripts/hmmer')))
from interproscan.scripts.hmmer.parser_out import parse

HMMER_OUT_MULTIPLE_SEQUENCES = "tests/unit_tests/test_inputs/0.0._.pfam._.out"
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

HMMER_OUT_MULTIPLE_DOMAINS = "tests/unit_tests/test_inputs/0.0._.gene3d._.out"
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

HMMER_OUT_MULTIPLESEQ_MULTIPLEDOM = "tests/unit_tests/test_inputs/0.0._.panther._.out"
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


def test_check_output_types():
    result = parse(HMMER_OUT_MULTIPLE_SEQUENCES)

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
def test_parse_multiple_sequences(expected_result):
    result = parse(HMMER_OUT_MULTIPLE_SEQUENCES)
    assert result == expected_result


@pytest.mark.parametrize("expected_result", [EXPECTED_RESULT_MULTIPLE_DOMAINS])
def test_parse_multiple_domains(expected_result):
    result = parse(HMMER_OUT_MULTIPLE_DOMAINS)
    assert result == expected_result


@pytest.mark.parametrize("expected_result", [EXPECTED_RESULT_MULTIPLESEQ_MULTIPLEDOM])
def test_parse_multiple_sequences_multiple_domains(expected_result):
    result = parse(HMMER_OUT_MULTIPLESEQ_MULTIPLEDOM)
    assert result == expected_result
