import pytest
from microprotein_lm.insdseq_parser import INSDSeqParser, reverse_complement

SAMPLE_INSDSEQ_XML = """<?xml version="1.0" encoding="UTF-8"?>
<INSDSet>
  <INSDSeq>
    <INSDSeq_locus>HUMAN_CDS</INSDSeq_locus>
    <INSDSeq_length>120</INSDSeq_length>
    <INSDSeq_moltype>mRNA</INSDSeq_moltype>
    <INSDSeq_definition>Homo sapiens microprotein CDS</INSDSeq_definition>
    <INSDSeq_primary-accession>NM_001000.1</INSDSeq_primary-accession>
    <INSDSeq_organism>Homo sapiens</INSDSeq_organism>
    <INSDSeq_taxonomy>Eukaryota; Metazoa; Chordata; Mammalia; Primates; Hominidae; Homo</INSDSeq_taxonomy>
    <INSDSeq_sequence>ACGTACGTACGTATGGCTTACTGAGGGCCCAAA</INSDSeq_sequence>
    <INSDSeq_feature-table>
      <INSDFeature>
        <INSDFeature_key>source</INSDFeature_key>
        <INSDFeature_location>1..120</INSDFeature_location>
        <INSDFeature_quals>
          <INSDQualifier>
            <INSDQualifier_name>organelle</INSDQualifier_name>
            <INSDQualifier_value>mitochondrion</INSDQualifier_value>
          </INSDQualifier>
          <INSDQualifier>
            <INSDQualifier_name>db_xref</INSDQualifier_name>
            <INSDQualifier_value>taxon:9606</INSDQualifier_value>
          </INSDQualifier>
        </INSDFeature_quals>
      </INSDFeature>
      <INSDFeature>
        <INSDFeature_key>CDS</INSDFeature_key>
        <INSDFeature_location>13..21</INSDFeature_location>
        <INSDFeature_intervals>
          <INSDInterval>
            <INSDInterval_from>13</INSDInterval_from>
            <INSDInterval_to>21</INSDInterval_to>
            <INSDInterval_strand>plus</INSDInterval_strand>
          </INSDInterval>
        </INSDFeature_intervals>
        <INSDFeature_quals>
          <INSDQualifier>
            <INSDQualifier_name>gene</INSDQualifier_name>
            <INSDQualifier_value>MT-ATP8</INSDQualifier_value>
          </INSDQualifier>
          <INSDQualifier>
            <INSDQualifier_name>protein_id</INSDQualifier_name>
            <INSDQualifier_value>NP_001000.1</INSDQualifier_value>
          </INSDQualifier>
          <INSDQualifier>
            <INSDQualifier_name>product</INSDQualifier_name>
            <INSDQualifier_value>ATP synthase protein 8</INSDQualifier_value>
          </INSDQualifier>
          <INSDQualifier>
            <INSDQualifier_name>transl_table</INSDQualifier_name>
            <INSDQualifier_value>2</INSDQualifier_value>
          </INSDQualifier>
          <INSDQualifier>
            <INSDQualifier_name>translation</INSDQualifier_name>
            <INSDQualifier_value>MA</INSDQualifier_value>
          </INSDQualifier>
        </INSDFeature_quals>
      </INSDFeature>
    </INSDSeq_feature-table>
  </INSDSeq>
</INSDSet>
"""


def test_parse_insdseq_xml():
    parser = INSDSeqParser()
    records = parser.parse_string(SAMPLE_INSDSEQ_XML)

    assert len(records) == 1
    rec = records[0]
    assert rec['accession'] == "NM_001000.1"
    assert rec['organism'] == "Homo sapiens"
    assert rec['organelle'] == "mitochondrion"
    assert rec['sequence'] == "ACGTACGTACGTATGGCTTACTGAGGGCCCAAA"

    assert len(rec['cds_features']) == 1
    cds = rec['cds_features'][0]
    assert cds['gene'] == "MT-ATP8"
    assert cds['protein_id'] == "NP_001000.1"
    assert cds['product'] == "ATP synthase protein 8"
    assert cds['transl_table'] == 2
    assert cds['translation'] == "MA"
    assert cds['strand'] == "+"
    assert cds['exon_coordinates'] == [(13, 21)]


def test_reverse_complement():
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("AUGC") == "GCAU"


def test_utr_flanks_extraction():
    parser = INSDSeqParser()
    records = parser.parse_string(SAMPLE_INSDSEQ_XML)
    cds = records[0]['cds_features'][0]

    # Sequence: ACGTACGTACGT ATGGCTTAC TGAGGGCCCAAA
    # Pos 1-12: ACGTACGTACGT (5' UTR) -> RNA: ACGUACGUACGU
    # Pos 13-21: ATGGCTTAC (CDS)
    # Pos 22+: TGAGGGCCCAAA (3' UTR) -> RNA: UGAGGGCCCAAA
    assert cds['flanking_5prime_utr_30nt'] == "ACGUACGUACGU"
    assert cds['flanking_3prime_utr_30nt'] == "UGAGGGCCCAAA"
