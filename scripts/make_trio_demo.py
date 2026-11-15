"""从 samples/demo.fa 生成家系分相核查的演示 VCF（父、母、子各一份）。

运行：python3 scripts/make_trio_demo.py

演示场景（坐标与 REF 均由参考序列自动计算，保证合法）：
- 母源：pos5 与 pos20 之间夹着一个纯合位点，1 次切换落在哪段区间无法确定
  （三段相邻区间均为"可能发生"）；pos20 → pos21 两次相邻杂合强制 1 次切换
  （"必然发生"）。
- 父源：pos30 与 pos45 之间的 1 次切换落在哪段区间无法确定（"可能发生"）。
- 最小切换总次数 = 3（母源 2 次 + 父源 1 次）。
生成后用 app.trio.analyze_trio 自检上述期望，失败则报错退出。
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.bio import parse_fasta  # noqa: E402
from app.trio import analyze_trio, parse_trio_vcf  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
CHROM = "chrDemo"

EXPECTED_MIN_SWITCHES = 3
EXPECTED_INTERVALS = {
    # (pos_a, pos_b): (母源判定, 父源判定)
    (5, 10): ("maybe", "never"),
    (10, 12): ("maybe", "never"),
    (12, 20): ("maybe", "never"),
    (20, 21): ("must", "never"),
    (21, 25): ("never", "never"),
    (25, 30): ("never", "never"),
    (30, 35): ("never", "maybe"),
    (35, 45): ("never", "maybe"),
}


def header(sample: str) -> str:
    return (
        "##fileformat=VCFv4.2\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype">\n'
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n"
    )


def row(seq: str, pos: int, alt: str, gt: str, ps: int | None = None) -> str:
    ref = seq[pos - 1]
    if ps is None:
        return f"{CHROM}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}"
    return f"{CHROM}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT:PS\t{gt}:{ps}"


def main() -> None:
    ref = parse_fasta((ROOT / "samples" / "demo.fa").read_bytes())
    seq = ref.sequence

    mother_rows = [
        row(seq, 5, "C", "0|1", 100),   # 杂合：单倍型1 携带 C
        row(seq, 12, "G", "1|1"),       # 纯合：不限制单倍型选择
        row(seq, 20, "G", "1|0", 100),  # 杂合：单倍型0 携带 G
        row(seq, 21, "T", "0|1", 100),  # 杂合：单倍型1 携带 T
        row(seq, 35, "T", "0|1", 100),  # 杂合：单倍型1 携带 T
    ]
    father_rows = [
        row(seq, 10, "G", "0|1", 200),  # 杂合：单倍型1 携带 G
        row(seq, 25, "A", "1|1"),       # 纯合
        row(seq, 30, "C", "0|1", 200),  # 杂合：单倍型1 携带 C
        row(seq, 45, "C", "1|0", 200),  # 杂合：单倍型0 携带 C
    ]
    child_rows = [
        row(seq, 5, "C", "1|0", 300),   # 母源 C（母亲单倍型1）
        row(seq, 10, "G", "0|1", 300),  # 父源 G（父亲单倍型1）
        row(seq, 12, "G", "1|0", 300),  # 母源 G（母亲纯合，任选）
        row(seq, 20, "G", "1|0", 300),  # 母源 G（母亲单倍型0）
        row(seq, 21, "T", "1|0", 300),  # 母源 T（母亲单倍型1）
        row(seq, 25, "A", "0|1", 300),  # 父源 A（父亲纯合，任选）
        row(seq, 30, "C", "0|1", 300),  # 父源 C（父亲单倍型1）
        row(seq, 35, "T", "1|0", 300),  # 母源 T（母亲单倍型1）
        row(seq, 45, "C", "0|1", 300),  # 父源 C（父亲单倍型0）
    ]

    outputs = {
        "trio_father.vcf": ("FATHER", father_rows),
        "trio_mother.vcf": ("MOTHER", mother_rows),
        "trio_child.vcf": ("CHILD", child_rows),
    }
    parsed = {}
    for filename, (sample, rows) in outputs.items():
        text = header(sample) + "\n".join(rows) + "\n"
        (ROOT / "samples" / filename).write_text(text, encoding="utf-8")
        parsed[filename] = parse_trio_vcf(text.encode(), filename, ref)[0]
        print(f"wrote {len(rows)} records -> samples/{filename}")

    # 自检：演示数据必须复现期望的最小次数与区间判定
    result = analyze_trio(
        ref, parsed["trio_father.vcf"], parsed["trio_mother.vcf"], parsed["trio_child.vcf"]
    )
    assert result.min_switches == EXPECTED_MIN_SWITCHES, result.min_switches
    got = {(iv.pos_a, iv.pos_b): (iv.maternal, iv.paternal) for iv in result.intervals}
    assert got == EXPECTED_INTERVALS, got
    print(f"self-check ok: min switches = {result.min_switches}, "
          f"{len(result.intervals)} intervals as expected")


if __name__ == "__main__":
    main()
