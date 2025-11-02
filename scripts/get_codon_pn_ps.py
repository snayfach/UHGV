#!/usr/bin/env python

import gzip
import json
import re
from collections import Counter
from functools import lru_cache
from itertools import permutations, product

from Bio.Data import CodonTable
from pydantic.dataclasses import dataclass


class CodonPnPs:
    def __init__(
        self,
        codon_counts_dict: dict[str, int],
        genetic_code: int,
        alphabet: tuple[str, str, str, str] = ("A", "T", "G", "C"),
    ):
        self.codon_counts_dict = Counter(codon_counts_dict)
        self.most_frequent_codon = self.codon_counts_dict.most_common(1)[0][0]
        self.coverage = sum(codon_counts_dict.values())
        self.genetic_code = genetic_code
        self.alphabet = alphabet
        self.codon_table = self.get_codon_table(genetic_code)
        self.aa_counts_dict = self.get_aa_counts_dict()
        self.expected_tables = self.get_codon_s_n()
        self.codons_paths = self.get_all_codons_paths(alphabet)

    def get_aa_counts_dict(self) -> Counter[str]:
        aa_counts = {self.codon_table[k]: 0 for k in self.codon_counts_dict}
        for codon, count in self.codon_counts_dict.items():
            aa_counts[self.codon_table[codon]] += count
        return Counter(aa_counts)

    @staticmethod
    @lru_cache
    def get_codon_table(genetic_code) -> dict[str, str]:
        codon_table = CodonTable.generic_by_id[genetic_code].forward_table
        codon_table.update(
            {k: "*" for k in CodonTable.generic_by_id[genetic_code].stop_codons}
        )
        return {
            k: v
            for k, v in codon_table.items()
            if all(letter in ["A", "T", "G", "C"] for letter in k)
        }

    def is_synonymous(self, codon_1: str, codon_2: str) -> bool:
        return self.codon_table[codon_1] == self.codon_table[codon_2]

    @staticmethod
    @lru_cache
    def get_all_codons_paths(alphabet) -> dict[tuple[str, str], list[list[str]]]:
        all_codons = ["".join(c) for c in product(alphabet, repeat=3)]
        precomputed_paths = {}
        for codon_1 in all_codons:
            for codon_2 in all_codons:
                hamming_distance = sum(a != b for a, b in zip(codon_1, codon_2))
                if 1 <= hamming_distance <= 3:
                    differing_positions = [
                        i for i, (a, b) in enumerate(zip(codon_1, codon_2)) if a != b
                    ]
                    paths = []
                    for perm in permutations(differing_positions):
                        current_codon = list(codon_1)
                        path = [codon_1]
                        for position in perm:
                            current_codon[position] = codon_2[position]
                            path.append("".join(current_codon))
                        paths.append(path)
                    precomputed_paths[(codon_1, codon_2)] = paths
        return precomputed_paths

    def get_codon_s_n(self):
        n = {k: 0.0 for k in self.codon_table.keys()}
        for codon, aa in self.codon_table.items():
            for pos, base in enumerate(codon):
                if base not in self.alphabet:
                    raise ValueError(f"Invalid base {base} in codon {codon}")
                for letter in self.alphabet:
                    if letter != base:
                        mut = codon[:pos] + letter + codon[pos + 1 :]
                        if self.codon_table.get(mut, aa) != aa:
                            n[codon] += 1 / 3
        s = {k: 3 - v for (k, v) in n.items()}
        return (s, n)

    def get_codon_sd_nd(self, codon_1: str, codon_2: str) -> tuple[float, float]:
        if codon_1 == codon_2:
            return 0.0, 0.0
        hamming_distance = sum(a != b for a, b in zip(codon_1, codon_2))
        if hamming_distance == 1:
            return (1.0, 0.0) if self.is_synonymous(codon_1, codon_2) else (0.0, 1.0)
        paths = self.codons_paths[(codon_1, codon_2)]
        situations = len(paths)
        sd, nd = 0.0, 0.0
        for path in paths:
            for i in range(len(path) - 1):
                if self.is_synonymous(path[i], path[i + 1]):
                    sd += 1
                else:
                    nd += 1
        sd = sd / situations
        nd = nd / situations
        return sd, nd

    def _get_codon_pn_ps_specific_ref(self, ref_codon: str) -> tuple[float, float]:
        s_table, n_table = self.expected_tables
        s_table = s_table[ref_codon]
        n_table = n_table[ref_codon]
        ps, pn = 0.0, 0.0
        for obs in self.codon_counts_dict:
            if ref_codon == obs:
                continue
            _ps, _pn = self.get_codon_sd_nd(ref_codon, obs)
            ps += _ps * self.codon_counts_dict[obs] / self.coverage
            pn += _pn * self.codon_counts_dict[obs] / self.coverage
        ps = ps / s_table if s_table > 0 else 0.0
        pn = pn / n_table if n_table > 0 else 0.0
        return ps, pn

    def get_codon_pn_ps_single_ref(self) -> tuple[float, float]:
        return self._get_codon_pn_ps_specific_ref(self.most_frequent_codon)

    def get_codon_pn_ps_multiple_refs(self) -> tuple[float, float]:
        ps, pn = 0.0, 0.0
        for ref_codon, ref_codon_coverage in self.codon_counts_dict.items():
            _ps, _pn = self._get_codon_pn_ps_specific_ref(ref_codon)
            ps += _ps * ref_codon_coverage / self.coverage
            pn += _pn * ref_codon_coverage / self.coverage
        return ps, pn


@dataclass
class Row:
    entry_id: int
    unique_pos_identifier: int
    contig_name: str
    sample_id: str
    corresponding_gene_call: int
    codon_order_in_gene: int
    codon_number: int
    gene_length: int
    coverage: int
    AAA: int
    AAC: int
    AAG: int
    AAT: int
    ACA: int
    ACC: int
    ACG: int
    ACT: int
    AGA: int
    AGC: int
    AGG: int
    AGT: int
    ATA: int
    ATC: int
    ATG: int
    ATT: int
    CAA: int
    CAC: int
    CAG: int
    CAT: int
    CCA: int
    CCC: int
    CCG: int
    CCT: int
    CGA: int
    CGC: int
    CGG: int
    CGT: int
    CTA: int
    CTC: int
    CTG: int
    CTT: int
    GAA: int
    GAC: int
    GAG: int
    GAT: int
    GCA: int
    GCC: int
    GCG: int
    GCT: int
    GGA: int
    GGC: int
    GGG: int
    GGT: int
    GTA: int
    GTC: int
    GTG: int
    GTT: int
    TAA: int
    TAC: int
    TAG: int
    TAT: int
    TCA: int
    TCC: int
    TCG: int
    TCT: int
    TGA: int
    TGC: int
    TGG: int
    TGT: int
    TTA: int
    TTC: int
    TTG: int
    TTT: int
    reference: str
    consensus: str
    competing_codons: str
    departure_from_reference: float
    departure_from_consensus: float
    n2n1ratio: float
    entropy: float

    @property
    def competing_codons_tuple(self):
        return (self.competing_codons[0:3], self.competing_codons[3:6])

    @property
    def codon_dict(self):
        return Counter(
            {
                codon: getattr(self, codon)
                for codon in ["".join(c) for c in product("ACGT", repeat=3)]
            }
        )

    @property
    def codon_dict_dense(self):
        return Counter(
            {
                codon: getattr(self, codon)
                for codon in ["".join(c) for c in product("ACGT", repeat=3)]
                if getattr(self, codon) > 0
            }
        )

    @property
    def codon_dict_dense_filtered(self):
        return Counter(
            {
                codon: count
                for codon, count in self.codon_dict_dense.items()
                if count / getattr(self, self.consensus)
                > (1 / 3) ** (self.coverage ** (1 / 3) - 1.45) + 0.05
            }
        )


contig_sample_presence = set()
with open("contig_coverage.tsv") as fin:
    for line in fin:
        sample, contig, covered, length, mean_cov = line.strip("\n").split("\t")
        if int(covered) / int(length) >= 0.5 and float(mean_cov) >= 10:
            contig_sample_presence.add((contig, sample))

gene_coord_to_name_dict = {}
gene_name_to_genetic_code = {}
header_parser = re.compile(
    r"(.+)_(.+) # ([0-9]+) # ([0-9]+) # (-1|1) .+rbs_motif=(.+?)"
    r";.+;genetic_code=(.+?);gc_cont=(.+)"
)
with gzip.open("votus_mq_plus.faa.gz", "rt") as fin:
    for line in fin:
        if line.startswith(">"):
            header = line.strip("\n>")
            if match := header_parser.match(header):
                contig, gene, start, end, _, _, code, _ = match.groups()
            else:
                raise ValueError(f"Header does not match the expected format: {header}")
            gene = f"{contig}_{gene}"
            gene_coord_to_name_dict[(contig, int(start), int(end))] = gene
            gene_name_to_genetic_code[gene] = int(code)

gene_id_to_name_dict = {}
gene_name_to_coord_dict = {}
with open("votus_mq_plus_external_gene_call.tsv") as fin:
    next(fin)
    for line in fin:
        gene_id, contig, start, end, *_ = line.split()
        gene_name = gene_coord_to_name_dict[(contig, int(start) + 1, int(end))]
        gene_id_to_name_dict[int(gene_id)] = gene_name
        gene_name_to_coord_dict[gene_name] = (int(start) + 1, int(end))

with open("VARIABILITY_CDN.txt") as fin, open("CODON_PN_PS.tsv", "w") as fout:
    next(fin)
    fout.write(
        "\t".join(
            [
                "sample_id",
                "contig_name",
                "genetic_code",
                "gene_name",
                "gene_start",
                "gene_end",
                "site_coverage",
                "codon_order_in_gene",
                "n_competing_codons",
                "codon_counts",
                "aa_counts",
                "ps_multiple_refs",
                "pn_multiple_refs",
                "ps_single_ref",
                "pn_single_ref",
            ]
        )
        + "\n"
    )
    for line in fin:
        row = Row(*line.strip("\n").split("\t"))
        if not row.departure_from_consensus:
            continue
        if len(row.codon_dict_dense_filtered) <= 1:
            continue
        if row.coverage < 10:
            continue
        if (row.contig_name, row.sample_id) not in contig_sample_presence:
            continue
        gene_name = gene_id_to_name_dict[row.corresponding_gene_call]
        gene_start, gene_end = gene_name_to_coord_dict[gene_name]
        codon_pn_ps = CodonPnPs(
            row.codon_dict_dense_filtered, gene_name_to_genetic_code[gene_name]
        )
        if codon_pn_ps.coverage < 10:
            continue
        ps_multiple_refs, pn_multiple_refs = codon_pn_ps.get_codon_pn_ps_multiple_refs()
        ps_single_ref, pn_single_ref = codon_pn_ps.get_codon_pn_ps_single_ref()
        output_row = [
            row.sample_id,
            row.contig_name,
            codon_pn_ps.genetic_code,
            gene_name,
            gene_start,
            gene_end,
            codon_pn_ps.coverage,
            row.codon_order_in_gene,
            len(row.codon_dict_dense_filtered),
            json.dumps(dict(codon_pn_ps.codon_counts_dict.most_common())),
            json.dumps(dict(codon_pn_ps.aa_counts_dict.most_common())),
            ps_multiple_refs,
            pn_multiple_refs,
            ps_single_ref,
            pn_single_ref,
        ]
        fout.write("\t".join(map(str, output_row)) + "\n")
