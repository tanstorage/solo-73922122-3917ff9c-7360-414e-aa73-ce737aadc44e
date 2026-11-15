"""家系分相核查核心逻辑测试，仅依赖标准库：python -I -m unittest -v"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.bio import CheckError, parse_fasta  # noqa: E402
from app.trio import (  # noqa: E402
    MAYBE,
    MUST,
    NEVER,
    analyze_trio,
    parse_trio_vcf,
)

# 1-based: 1A 2C 3G 4T | 5A 6A 7A 8A | 9C 10G 11T 12T
SEQ = "ACGTAAAACGTT"
FASTA = f">chrT\n{SEQ}\n".encode()
HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"


def vcf(records, sample="S1", extra_header_cols=""):
    lines = [HEADER + sample + extra_header_cols]
    lines.extend(records)
    return ("\n".join(lines) + "\n").encode()


def rec(pos, ref, alt, gt, ps=None, chrom="chrT", fmt=None, sample=None):
    """构造一条记录；ps 给定时使用 GT:PS，fmt/sample 可完全自定义两列。"""
    if fmt is not None:
        return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\t{fmt}\t{sample}"
    if ps is not None:
        return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT:PS\t{gt}:{ps}"
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t{gt}"


def parse_person(records, sample="S1"):
    refobj = parse_fasta(FASTA)
    return parse_trio_vcf(vcf(records, sample), f"{sample}.vcf", refobj)


def analyze(father_recs, mother_recs, child_recs):
    refobj = parse_fasta(FASTA)
    father, _ = parse_person(father_recs, "F")
    mother, _ = parse_person(mother_recs, "M")
    child, _ = parse_person(child_recs, "C")
    return analyze_trio(refobj, father, mother, child)


class TrioParseTests(unittest.TestCase):
    def test_happy_path_returns_records_and_sample_name(self):
        records, sample = parse_person(
            [rec(2, "C", "G", "0|1", ps=7), rec(3, "G", "A", "1|1")], sample="MOTHER"
        )
        self.assertEqual(sample, "MOTHER")
        self.assertEqual(sorted(records), [2, 3])
        het = records[2]
        self.assertTrue(het.heterozygous)
        self.assertEqual(het.bases(), ("C", "G"))
        self.assertEqual(het.ps, "7")
        self.assertEqual(het.line_no, 3)
        hom = records[3]
        self.assertFalse(hom.heterozygous)
        self.assertEqual(hom.bases(), ("A", "A"))
        self.assertIsNone(hom.ps)

    def test_unphased_gt_rejected(self):
        with self.assertRaises(CheckError) as ctx:
            parse_person([rec(2, "C", "G", "0/1", ps=7)])
        self.assertIn("分相", str(ctx.exception))

    def test_missing_gt_rejected(self):
        for bad in (".|1", "0|.", ".|."):
            with self.subTest(gt=bad):
                with self.assertRaises(CheckError):
                    parse_person([rec(2, "C", "G", bad, ps=7)])

    def test_wrong_ploidy_rejected(self):
        for bad in ("0", "0|1|1"):
            with self.subTest(gt=bad):
                with self.assertRaises(CheckError):
                    parse_person([rec(2, "C", "G", bad, ps=7)])

    def test_allele_index_out_of_range_rejected(self):
        with self.assertRaises(CheckError):
            parse_person([rec(2, "C", "G", "0|2", ps=7)])

    def test_multi_alt_rejected(self):
        with self.assertRaises(CheckError) as ctx:
            parse_person([rec(2, "C", "G,T", "0|1", ps=7)])
        self.assertIn("单 ALT", str(ctx.exception))

    def test_indel_and_mnv_rejected(self):
        for ref, alt in (("A", "AA"), ("AC", "A"), ("AC", "GT")):
            with self.subTest(ref=ref, alt=alt):
                with self.assertRaises(CheckError) as ctx:
                    parse_person([rec(1, ref, alt, "0|1", ps=7)])
                self.assertIn("单碱基替换", str(ctx.exception))

    def test_ref_equals_alt_rejected(self):
        with self.assertRaises(CheckError):
            parse_person([rec(2, "C", "C", "0|1", ps=7)])

    def test_symbolic_and_non_acgt_alt_rejected(self):
        for alt in (".", "<DEL>", "g", "G*"):
            with self.subTest(alt=alt):
                with self.assertRaises(CheckError):
                    parse_person([rec(2, "C", alt, "0|1", ps=7)])

    def test_ref_mismatch_rejected(self):
        with self.assertRaises(CheckError) as ctx:
            parse_person([rec(2, "G", "T", "0|1", ps=7)])  # pos2 参考实为 C
        self.assertIn("REF", str(ctx.exception))

    def test_chrom_mismatch_rejected(self):
        with self.assertRaises(CheckError):
            parse_person([rec(2, "C", "G", "0|1", ps=7, chrom="chrOther")])

    def test_pos_out_of_range_rejected(self):
        with self.assertRaises(CheckError):
            parse_person([rec(99, "C", "G", "0|1", ps=7)])

    def test_duplicate_pos_rejected(self):
        data = [rec(2, "C", "G", "0|1", ps=7), rec(2, "C", "T", "0|1", ps=7)]
        with self.assertRaises(CheckError) as ctx:
            parse_person(data)
        self.assertIn("至多一条", str(ctx.exception))

    def test_het_without_ps_rejected(self):
        with self.assertRaises(CheckError) as ctx:
            parse_person([rec(2, "C", "G", "0|1")])  # FORMAT 只有 GT
        self.assertIn("PS", str(ctx.exception))

    def test_het_with_missing_ps_value_rejected(self):
        with self.assertRaises(CheckError):
            parse_person([rec(2, "C", "G", "0|1", fmt="GT:PS", sample="0|1:.")])

    def test_het_records_must_share_one_ps(self):
        ok = [rec(2, "C", "G", "0|1", ps=7), rec(3, "G", "A", "1|0", ps=7)]
        records, _ = parse_person(ok)
        self.assertEqual(len(records), 2)

        bad = [rec(2, "C", "G", "0|1", ps=7), rec(3, "G", "A", "1|0", ps=8)]
        with self.assertRaises(CheckError) as ctx:
            parse_person(bad)
        self.assertIn("PS 不一致", str(ctx.exception))

    def test_homozygous_records_do_not_need_ps(self):
        records, _ = parse_person(
            [rec(2, "C", "G", "1|1"), rec(3, "G", "A", "0|1", ps=7)]
        )
        self.assertEqual(len(records), 2)

    def test_homozygous_ps_not_part_of_consistency_check(self):
        # 纯合记录上的 PS 不参与“同一 PS”判定
        records, _ = parse_person(
            [rec(2, "C", "G", "1|1", ps=99), rec(3, "G", "A", "0|1", ps=7)]
        )
        self.assertEqual(len(records), 2)

    def test_homozygous_ref_record_kept(self):
        records, _ = parse_person([rec(2, "C", "G", "0|0")])
        self.assertEqual(records[2].bases(), ("C", "C"))

    def test_multi_sample_header_rejected(self):
        refobj = parse_fasta(FASTA)
        data = vcf([rec(2, "C", "G", "0|1", ps=7)], extra_header_cols="\tS2")
        with self.assertRaises(CheckError):
            parse_trio_vcf(data, "x.vcf", refobj)

    def test_missing_header_rejected(self):
        refobj = parse_fasta(FASTA)
        with self.assertRaises(CheckError):
            parse_trio_vcf((rec(2, "C", "G", "0|1", ps=7) + "\n").encode(), "x.vcf", refobj)

    def test_format_not_leading_gt_rejected(self):
        with self.assertRaises(CheckError):
            parse_person([rec(2, "C", "G", "0|1", fmt="PS:GT", sample="7:0|1")])

    def test_record_missing_sample_column_rejected(self):
        refobj = parse_fasta(FASTA)
        data = vcf(["chrT\t2\t.\tC\tG\t.\t.\t.\tGT"])
        with self.assertRaises(CheckError):
            parse_trio_vcf(data, "x.vcf", refobj)

    def test_ps_from_other_format_position(self):
        records, _ = parse_person(
            [rec(2, "C", "G", "0|1", fmt="GT:DP:PS", sample="0|1:30:7")]
        )
        self.assertEqual(records[2].ps, "7")


class TrioAnalyzeTests(unittest.TestCase):
    def test_single_site_no_interval(self):
        result = analyze(
            [rec(2, "C", "G", "0|1", ps=5)],
            [],
            [rec(2, "C", "G", "0|1", ps=9)],
        )
        self.assertEqual(result.min_switches, 0)
        self.assertEqual(result.intervals, [])
        self.assertEqual([s.pos for s in result.sites], [2])

    def test_empty_trio(self):
        result = analyze([], [], [])
        self.assertEqual(result.min_switches, 0)
        self.assertEqual(result.sites, [])
        self.assertEqual(result.intervals, [])

    def test_forced_maternal_switch(self):
        # 母亲 pos2 0|1（hap1=G）、pos3 1|0（hap0=T）；
        # 子代母源依次取 G、T → 母源必须在 (2,3) 切换；父源全自由 → 不发生
        result = analyze(
            [],
            [rec(2, "C", "G", "0|1", ps=7), rec(3, "G", "T", "1|0", ps=7)],
            [rec(2, "C", "G", "1|0", ps=9), rec(3, "G", "T", "1|0", ps=9)],
        )
        self.assertEqual(result.min_switches, 1)
        self.assertEqual(len(result.intervals), 1)
        iv = result.intervals[0]
        self.assertEqual((iv.pos_a, iv.pos_b), (2, 3))
        self.assertEqual(iv.maternal, MUST)
        self.assertEqual(iv.paternal, NEVER)

    def test_forced_paternal_switch(self):
        result = analyze(
            [rec(2, "C", "G", "0|1", ps=3), rec(3, "G", "A", "1|0", ps=3)],
            [],
            [rec(2, "C", "G", "0|1", ps=9), rec(3, "G", "A", "0|1", ps=9)],
        )
        self.assertEqual(result.min_switches, 1)
        iv = result.intervals[0]
        self.assertEqual(iv.maternal, NEVER)
        self.assertEqual(iv.paternal, MUST)

    def test_switch_may_fall_in_either_interval(self):
        # 母亲 pos2 强制 hap1、pos9 强制 hap0，中间 pos5 纯合（自由选择）：
        # 唯一的 1 次切换落在 (2,5) 或 (5,9) 都可以 → 两段均为“可能”
        result = analyze(
            [],
            [
                rec(2, "C", "G", "0|1", ps=7),
                rec(5, "A", "G", "1|1"),
                rec(9, "C", "T", "0|1", ps=7),
            ],
            [rec(2, "C", "G", "1|0", ps=9), rec(5, "A", "G", "1|0", ps=9)],
        )
        self.assertEqual(result.min_switches, 1)
        self.assertEqual(
            [(iv.maternal, iv.paternal) for iv in result.intervals],
            [(MAYBE, NEVER), (MAYBE, NEVER)],
        )

    def test_unrecorded_child_still_constrains_parent(self):
        # 子代完全未记录（按参考纯合）：母亲 pos2 须给 C（hap0）、pos3 须给 G（hap1）
        result = analyze(
            [],
            [rec(2, "C", "G", "0|1", ps=1), rec(3, "G", "A", "1|0", ps=1)],
            [],
        )
        self.assertEqual(result.min_switches, 1)
        self.assertEqual(result.intervals[0].maternal, MUST)

    def test_maternal_and_paternal_switches_add_up(self):
        # 同一区间母源、父源各强制 1 次切换（父母携带相同 ALT，子代 1|1 继承双方 ALT）
        result = analyze(
            [rec(2, "C", "G", "0|1", ps=8), rec(3, "G", "T", "1|0", ps=8)],
            [rec(2, "C", "G", "0|1", ps=7), rec(3, "G", "T", "1|0", ps=7)],
            [rec(2, "C", "G", "1|1", ps=9), rec(3, "G", "T", "1|1", ps=9)],
        )
        self.assertEqual(result.min_switches, 2)
        iv = result.intervals[0]
        self.assertEqual(iv.maternal, MUST)
        self.assertEqual(iv.paternal, MUST)

    def test_sites_sorted_by_coordinate(self):
        # 文件内记录顺序不影响按坐标排序
        result = analyze(
            [rec(3, "G", "A", "0|1", ps=8), rec(2, "C", "G", "0|1", ps=8)],
            [],
            [rec(3, "G", "A", "0|1", ps=9), rec(2, "C", "G", "0|1", ps=9)],
        )
        self.assertEqual([s.pos for s in result.sites], [2, 3])
        self.assertEqual((result.intervals[0].pos_a, result.intervals[0].pos_b), (2, 3))

    def test_multi_allelic_site_across_files_compatible(self):
        # 同一位置不同文件 ALT 不同：母亲 C>T、父亲 C>G，子代母源 T、父源 C
        result = analyze(
            [rec(2, "C", "G", "0|1", ps=2)],
            [rec(2, "C", "T", "0|1", ps=1)],
            [rec(2, "C", "T", "1|0", ps=3)],
        )
        self.assertEqual(result.min_switches, 0)

    def test_child_alt_not_in_parents_fails(self):
        with self.assertRaises(CheckError) as ctx:
            analyze(
                [rec(2, "C", "G", "0|1", ps=2)],
                [rec(2, "C", "T", "0|1", ps=1)],
                [rec(2, "C", "A", "1|0", ps=3)],  # 母源 A 无人可提供
            )
        self.assertIn("无法来自母亲", str(ctx.exception))

    def test_maternal_allele_impossible_fails(self):
        with self.assertRaises(CheckError) as ctx:
            analyze([], [], [rec(2, "C", "G", "1|0", ps=1)])  # 母亲未记录：C/C
        self.assertIn("无法来自母亲", str(ctx.exception))

    def test_paternal_allele_impossible_fails(self):
        with self.assertRaises(CheckError) as ctx:
            analyze([], [], [rec(2, "C", "G", "0|1", ps=1)])  # 父亲未记录：C/C
        self.assertIn("无法来自父亲", str(ctx.exception))

    def test_homozygous_parent_excludes_ref_child(self):
        # 母亲 1|1，子代未记录（母源=REF）→ 无法继承
        with self.assertRaises(CheckError):
            analyze([], [rec(2, "C", "G", "1|1")], [])

    def test_failure_message_contains_pos(self):
        with self.assertRaises(CheckError) as ctx:
            analyze([], [], [rec(9, "C", "T", "1|0", ps=1)])
        self.assertIn("POS=9", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
