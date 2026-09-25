"""Native NCBI INSDSeq XML DOM Parser.

Extracts features, coordinates, qualifiers, organelles, strands, protein IDs,
gene symbols, product names, and 5'/3' UTR flanks directly from INSDSeq XML payloads.
"""

import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union


def reverse_complement(seq: str) -> str:
    """Computes reverse complement of a nucleotide string (DNA/RNA)."""
    if 'U' in seq.upper():
        comp = str.maketrans("ACGUacguNn", "UGCAugcaNn")
    else:
        comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]


class INSDSeqParser:
    """Parser for INSDSeq XML format payloads."""

    def parse_string(self, xml_content: str) -> List[Dict]:
        """Parses INSDSeq XML string into structured python dicts."""
        root = ET.fromstring(xml_content.strip())
        return self._parse_root(root)

    def parse_file(self, xml_path: Union[str, Path]) -> List[Dict]:
        """Parses INSDSeq XML file into structured python dicts."""
        tree = ET.parse(Path(xml_path))
        return self._parse_root(tree.getroot())

    def _parse_root(self, root: ET.Element) -> List[Dict]:
        records: List[Dict] = []
        # Support root being <INSDSet> or single <INSDSeq>
        seq_elements = root.findall('INSDSeq') if root.tag == 'INSDSet' else ([root] if root.tag == 'INSDSeq' else root.findall('.//INSDSeq'))

        for seq_elem in seq_elements:
            rec = self._parse_insdseq(seq_elem)
            records.append(rec)
        return records

    def _parse_insdseq(self, seq_elem: ET.Element) -> Dict:
        def get_text(tag: str) -> str:
            elem = seq_elem.find(tag)
            return elem.text.strip() if elem is not None and elem.text else ""

        accession = get_text('INSDSeq_primary-accession') or get_text('INSDSeq_accession-version')
        version = get_text('INSDSeq_accession-version') or accession
        locus = get_text('INSDSeq_locus')
        organism = get_text('INSDSeq_organism')
        taxonomy = get_text('INSDSeq_taxonomy')
        definition = get_text('INSDSeq_definition')
        sequence = get_text('INSDSeq_sequence').upper()
        length_val = get_text('INSDSeq_length')
        length = int(length_val) if length_val.isdigit() else len(sequence)

        # Parse feature table
        features: List[Dict] = []
        cds_features: List[Dict] = []
        organelle = ""

        feat_table_elem = seq_elem.find('INSDSeq_feature-table')
        if feat_table_elem is not None:
            for feat_elem in feat_table_elem.findall('INSDFeature'):
                key = feat_elem.findtext('INSDFeature_key', '').strip()
                location = feat_elem.findtext('INSDFeature_location', '').strip()

                # Parse intervals
                intervals: List[Tuple[int, int]] = []
                strand = '+'
                intervals_elem = feat_elem.find('INSDFeature_intervals')
                if intervals_elem is not None:
                    for int_elem in intervals_elem.findall('INSDInterval'):
                        int_from = int_elem.findtext('INSDInterval_from', '0')
                        int_to = int_elem.findtext('INSDInterval_to', '0')
                        int_strand = int_elem.findtext('INSDInterval_strand', '')
                        if int_strand.lower() in ('minus', 'complement', '-'):
                            strand = '-'
                        elif int_elem.find('INSDInterval_strand') is not None and int_elem.find('INSDInterval_strand').get('value') == 'minus':
                            strand = '-'

                        if int_from.isdigit() and int_to.isdigit():
                            intervals.append((int(int_from), int(int_to)))

                if 'complement' in location.lower():
                    strand = '-'

                # Parse qualifiers
                quals: Dict[str, List[str]] = {}
                quals_elem = feat_elem.find('INSDFeature_quals')
                if quals_elem is not None:
                    for q_elem in quals_elem.findall('INSDQualifier'):
                        q_name = q_elem.findtext('INSDQualifier_name', '').strip()
                        q_val = q_elem.findtext('INSDQualifier_value', '').strip()
                        if q_name:
                            quals.setdefault(q_name, []).append(q_val)

                if key == 'source':
                    organelle = quals.get('organelle', [''])[0] or quals.get('mol_type', [''])[0]

                protein_id = quals.get('protein_id', [''])[0]
                gene = quals.get('gene', [''])[0]
                gene_synonyms = quals.get('gene_synonym', [])
                product = quals.get('product', [''])[0]
                transl_table_str = quals.get('transl_table', ['1'])[0]
                transl_table = int(transl_table_str) if transl_table_str.isdigit() else 1
                translation = quals.get('translation', [''])[0]

                # Extract 5' and 3' UTR flanks if full sequence context is available
                flank_5prime = ""
                flank_3prime = ""
                if sequence and intervals:
                    flank_5prime, flank_3prime = self._extract_utr_flanks(sequence, intervals, strand)

                feat_dict = {
                    'key': key,
                    'location': location,
                    'strand': strand,
                    'intervals': intervals,
                    'exon_coordinates': intervals,
                    'qualifiers': quals,
                    'protein_id': protein_id,
                    'gene': gene,
                    'gene_synonyms': gene_synonyms,
                    'product': product,
                    'transl_table': transl_table,
                    'translation': translation,
                    'flanking_5prime_utr_30nt': flank_5prime,
                    'flanking_3prime_utr_30nt': flank_3prime,
                }
                features.append(feat_dict)
                if key == 'CDS':
                    cds_features.append(feat_dict)

        return {
            'accession': accession,
            'version': version,
            'locus': locus,
            'organism': organism,
            'taxonomy': taxonomy,
            'definition': definition,
            'sequence': sequence,
            'length': length,
            'organelle': organelle,
            'features': features,
            'cds_features': cds_features,
        }

    def _extract_utr_flanks(
        self,
        sequence: str,
        intervals: List[Tuple[int, int]],
        strand: str,
        flank_length: int = 30
    ) -> Tuple[str, str]:
        """Extracts 5' UTR and 3' UTR nucleotide flanks from sequence given CDS intervals."""
        if not sequence or not intervals:
            return "", ""

        seq_len = len(sequence)
        if strand == '+':
            min_from = min(i[0] for i in intervals)
            max_to = max(i[1] for i in intervals)

            # 5' UTR is before min_from (1-based)
            start_5 = max(0, min_from - 1 - flank_length)
            end_5 = max(0, min_from - 1)
            utr_5 = sequence[start_5:end_5]

            # 3' UTR is after max_to (1-based)
            start_3 = min(seq_len, max_to)
            end_3 = min(seq_len, max_to + flank_length)
            utr_3 = sequence[start_3:end_3]
        else:
            # On minus strand, CDS is on reverse complement
            # The lowest coordinate in genomic sequence is the 3' end of CDS
            # The highest coordinate in genomic sequence is the 5' start of CDS
            min_from = min(i[0] for i in intervals)
            max_to = max(i[1] for i in intervals)

            # 5' UTR is downstream (higher genomic coordinate) on minus strand
            start_5 = max_to
            end_5 = min(seq_len, max_to + flank_length)
            utr_5_genomic = sequence[start_5:end_5]
            utr_5 = reverse_complement(utr_5_genomic)

            # 3' UTR is upstream (lower genomic coordinate) on minus strand
            start_3 = max(0, min_from - 1 - flank_length)
            end_3 = max(0, min_from - 1)
            utr_3_genomic = sequence[start_3:end_3]
            utr_3 = reverse_complement(utr_3_genomic)

        # Convert T to U for RNA representation
        return utr_5.replace('T', 'U'), utr_3.replace('T', 'U')
