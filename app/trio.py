"""家系分相核查的核心逻辑（仅依赖标准库）。

输入：同一条染色体参考 + 父、母、子三份单样本 VCF。仅接受：
- 单 ALT 的单碱基替换（REF / ALT 均为单个 A/C/G/T 碱基）；
- 二倍体分相 GT（形如 0|1、1|0、1|1、0|0）；
- 同一文件内所有杂合记录（0|1、1|0）带有相同的 PS（相位集）；
- 每人每个位点（POS）至多一条记录。

子代 GT 左拷贝视为母源、右拷贝视为父源；某人未记录的位点按参考纯合处理。
在所有满足孟德尔遗传（子代母源等位来自母亲两条单倍型之一、父源等位来自
父亲两条单倍型之一）的继承路径中，用前向-后向动态规划寻找母源 + 父源
单倍型切换总次数最少的路径；再对每段相邻位点区间，判定母源 / 父源切换
在所有最优路径中必然发生、可能发生还是不发生。
任一位点无法继承或输入出错即整次失败（抛 CheckError），不给局部结论。
"""

from __future__ import annotations

from dataclasses import dataclass

from .bio import CheckError, Reference

_BASES = frozenset("ACGT")

# 区间切换判定取值
MUST = "must"    # 必然发生：所有最优路径在该区间都切换
MAYBE = "maybe"  # 可能发生：部分最优路径在该区间切换
NEVER = "never"  # 不发生：没有最优路径在该区间切换


@dataclass(frozen=True)
class PersonRecord:
    """一份 VCF 中某个位点的分相调用。"""

    pos: int
    ref: str
    alt: str
    gt: str
    a1: int  # 左拷贝等位索引（0=REF，1=ALT）
    a2: int  # 右拷贝等位索引
    ps: str | None
    line_no: int

    @property
    def heterozygous(self) -> bool:
        return self.a1 != self.a2

    def bases(self) -> tuple[str, str]:
        """左、右拷贝对应的碱基。"""
        return (
            self.ref if self.a1 == 0 else self.alt,
            self.ref if self.a2 == 0 else self.alt,
        )


@dataclass(frozen=True)
class SiteRow:
    """结果中一个位点上三人的调用（None 表示未记录，按参考纯合处理）。"""

    pos: int
    ref: str
    father: PersonRecord | None
    mother: PersonRecord | None
    child: PersonRecord | None


@dataclass(frozen=True)
class IntervalVerdict:
    """一段相邻位点区间内，母源 / 父源切换在所有最优路径中的必然性。"""

    pos_a: int
    pos_b: int
    maternal: str  # MUST / MAYBE / NEVER
    paternal: str


@dataclass(frozen=True)
class TrioResult:
    sites: list[SiteRow]
    intervals: list[IntervalVerdict]
    min_switches: int


# --------------------------------------------------------------------------- #
# GT
# --------------------------------------------------------------------------- #
def _parse_phased_gt(gt: str, line_no: int) -> tuple[int, int]:
    """解析二倍体分相 GT（单 ALT，等位仅允许 0/1），返回 (左, 右) 等位索引。"""
    if not gt:
        raise CheckError(f"第 {line_no} 行：GT 为空。")
    if "/" in gt:
        raise CheckError(
            f"第 {line_no} 行：家系核查要求二倍体分相 GT（'|'），发现未分相写法：{gt}"
        )
    parts = gt.split("|")
    if len(parts) != 2:
        raise CheckError(f"第 {line_no} 行：仅接受二倍体分相 GT（两个等位），收到 {gt!r}。")

    alleles: list[int] = []
    for part in parts:
        if part in ("", "."):
            raise CheckError(f"第 {line_no} 行：GT 含缺失等位，本工具不接受缺失 GT：{gt}")
        try:
            idx = int(part)
        except ValueError:
            raise CheckError(f"第 {line_no} 行：GT 等位 {part!r} 不是整数索引。")
        if idx < 0:
            raise CheckError(f"第 {line_no} 行：GT 等位索引不能为负：{gt}")
        if idx > 1:
            raise CheckError(
                f"第 {line_no} 行：GT 等位索引 {idx} 超出 ALT 数量（单 ALT，仅允许 0/1）。"
            )
        alleles.append(idx)
    return alleles[0], alleles[1]


# --------------------------------------------------------------------------- #
# VCF
# --------------------------------------------------------------------------- #
def _parse_trio_record(fields: list[str], line_no: int, refobj: Reference) -> PersonRecord:
    """解析一条家系核查 VCF 数据行。"""
    if len(fields) < 8:
        raise CheckError(f"第 {line_no} 行：列数不足 8 列，不是合法 VCF 记录。")

    chrom, pos_s, _vid, ref, alt_s = fields[:5]
    if chrom != refobj.name:
        raise CheckError(
            f"第 {line_no} 行：CHROM={chrom!r} 与参考 FASTA 名称 {refobj.name!r} 不一致。"
        )

    try:
        pos = int(pos_s)
    except ValueError:
        raise CheckError(f"第 {line_no} 行：POS={pos_s!r} 不是整数。")
    if pos < 1 or pos > len(refobj.sequence):
        raise CheckError(
            f"第 {line_no} 行：POS={pos} 超出参考范围（参考长度 {len(refobj.sequence)}）。"
        )

    if not ref or any(ch not in _BASES for ch in ref):
        raise CheckError(f"第 {line_no} 行：REF 等位 {ref!r} 不是纯 A/C/G/T 序列。")
    if len(ref) != 1:
        raise CheckError(
            f"第 {line_no} 行：家系核查仅接受单碱基替换，REF 必须为单个碱基"
            f"（收到 {ref!r}，不支持缺失 / MNV）。"
        )
    if "," in alt_s:
        raise CheckError(
            f"第 {line_no} 行：家系核查仅接受单 ALT 记录，收到 ALT={alt_s!r}。"
        )
    if not alt_s or any(ch not in _BASES for ch in alt_s):
        raise CheckError(
            f"第 {line_no} 行：ALT {alt_s!r} 为空或含非 A/C/G/T 字符"
            "（本工具不接受 *、<DEL>、符号化等位及缺失占位）。"
        )
    if len(alt_s) != 1:
        raise CheckError(
            f"第 {line_no} 行：家系核查仅接受单碱基替换，ALT 必须为单个碱基"
            f"（收到 {alt_s!r}，不支持插入 / MNV）。"
        )
    if alt_s == ref:
        raise CheckError(f"第 {line_no} 行：REF 与 ALT 相同（{ref}），不是有效的单碱基替换。")

    ref_base = refobj.sequence[pos - 1]
    if ref_base != ref:
        raise CheckError(
            f"第 {line_no} 行：REF={ref} 与参考序列在 pos={pos} 处的 {ref_base} 不一致，"
            "整次核查失败。"
        )

    # FORMAT/GT/PS
    if len(fields) < 10:
        raise CheckError(f"第 {line_no} 行：缺少 FORMAT 与样本列，本工具要求单样本二倍体分相 GT。")
    fmt = fields[8]
    fmt_keys = fmt.split(":")
    if not fmt_keys or fmt_keys[0] != "GT":
        raise CheckError(f"第 {line_no} 行：FORMAT 必须以 GT 开头，实际为 {fmt!r}。")
    sample_values = fields[9].split(":")
    if len(sample_values) > len(fmt_keys):
        raise CheckError(f"第 {line_no} 行：样本列的子字段多于 FORMAT 定义。")
    gt = sample_values[0] if sample_values else ""
    a1, a2 = _parse_phased_gt(gt, line_no)

    ps: str | None = None
    if "PS" in fmt_keys:
        ps_idx = fmt_keys.index("PS")
        if ps_idx < len(sample_values):
            candidate = sample_values[ps_idx]
            if candidate and candidate != ".":
                ps = candidate
    if a1 != a2 and ps is None:
        raise CheckError(
            f"第 {line_no} 行：杂合记录（{gt}）缺少 PS；家系核查要求同一文件的"
            "杂合记录属于同一相位集（FORMAT 含 PS 且取值非缺失）。"
        )

    return PersonRecord(
        pos=pos, ref=ref, alt=alt_s, gt=gt, a1=a1, a2=a2, ps=ps, line_no=line_no
    )


def parse_trio_vcf(
    data: bytes, filename: str, refobj: Reference
) -> tuple[dict[int, PersonRecord], str]:
    """解析家系核查用的单样本 VCF，返回 ({POS: PersonRecord}, 样本名)。

    任何硬性错误（含同位点多条记录、杂合记录 PS 不一致）都抛 CheckError。
    """
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise CheckError(f"{filename} 不是 UTF-8 文本（字节 {exc.start} 处解码失败）。")

    lines = text.splitlines()
    header_idx = -1
    for i, line in enumerate(lines):
        if line.startswith("#CHROM"):
            header_idx = i
            break
    if header_idx < 0:
        raise CheckError(f"{filename}：缺少 '#CHROM' 表头行，不是合法 VCF。")

    header_fields = lines[header_idx].rstrip("\r").split("\t")
    if len(header_fields) < 10:
        raise CheckError(
            f"{filename}：'#CHROM' 表头只有 {len(header_fields)} 列，"
            "本工具要求恰好一个样本列（共 10 列）。"
        )
    if len(header_fields) > 10:
        raise CheckError(
            f"{filename}：'#CHROM' 表头有 {len(header_fields) - 9} 个样本列，"
            "本工具仅接受单样本 VCF。"
        )
    sample_name = header_fields[9]

    records: dict[int, PersonRecord] = {}
    het_ps: str | None = None
    het_ps_line: int | None = None
    for offset, raw in enumerate(lines[header_idx + 1:], start=header_idx + 2):
        line = raw.rstrip("\r")
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        rec = _parse_trio_record(fields, offset, refobj)

        previous = records.get(rec.pos)
        if previous is not None:
            raise CheckError(
                f"{filename}：位点 POS={rec.pos} 在第 {previous.line_no} 行与第 {offset} 行"
                "各有一条记录，每人每个位点至多一条，整次核查失败。"
            )
        if rec.heterozygous:
            if het_ps is None:
                het_ps, het_ps_line = rec.ps, offset
            elif rec.ps != het_ps:
                raise CheckError(
                    f"{filename}：杂合记录的 PS 不一致（第 {het_ps_line} 行为 {het_ps}，"
                    f"第 {offset} 行为 {rec.ps}），同一文件的杂合记录须属于同一 PS，"
                    "整次核查失败。"
                )
        records[rec.pos] = rec

    return records, sample_name


# --------------------------------------------------------------------------- #
# 最少切换路径
# --------------------------------------------------------------------------- #
def _verdict(flags: set[bool]) -> str:
    """根据最优边集合上是否出现切换，给出 必然 / 可能 / 不发生。"""
    if flags == {True}:
        return MUST
    if flags == {False}:
        return NEVER
    return MAYBE


def analyze_trio(
    refobj: Reference,
    father: dict[int, PersonRecord],
    mother: dict[int, PersonRecord],
    child: dict[int, PersonRecord],
) -> TrioResult:
    """在三人等位相容的继承路径中最小化亲本单倍型切换总次数。

    状态为每个位点上（母源单倍型, 父源单倍型）的选择（取值 0/1，即亲本
    分相 GT 的左 / 右拷贝）；相邻位点间每换一次单倍型计 1 次切换，
    母源与父源合计。任一位点等位无法继承时抛 CheckError（整次失败）。
    """
    positions = sorted(set(father) | set(mother) | set(child))

    def bases_at(records: dict[int, PersonRecord], pos: int) -> tuple[str, str]:
        rec = records.get(pos)
        if rec is None:  # 未记录：按参考纯合
            base = refobj.sequence[pos - 1]
            return (base, base)
        return rec.bases()

    # 每个位点允许的状态：子代母源等位必须落在所选母亲单倍型上，父源同理
    allowed: list[list[tuple[int, int]]] = []
    for pos in positions:
        mb = bases_at(mother, pos)
        fb = bases_at(father, pos)
        cb = bases_at(child, pos)
        mat_choices = [h for h in (0, 1) if mb[h] == cb[0]]
        pat_choices = [h for h in (0, 1) if fb[h] == cb[1]]
        if not mat_choices:
            raise CheckError(
                f"位点 POS={pos}（REF={refobj.sequence[pos - 1]}）：子代母源等位 {cb[0]} "
                f"无法来自母亲的两条单倍型（{mb[0]} / {mb[1]}），无法继承，整次核查失败。"
            )
        if not pat_choices:
            raise CheckError(
                f"位点 POS={pos}（REF={refobj.sequence[pos - 1]}）：子代父源等位 {cb[1]} "
                f"无法来自父亲的两条单倍型（{fb[0]} / {fb[1]}），无法继承，整次核查失败。"
            )
        allowed.append([(m, p) for m in mat_choices for p in pat_choices])

    sites = [
        SiteRow(
            pos=pos,
            ref=refobj.sequence[pos - 1],
            father=father.get(pos),
            mother=mother.get(pos),
            child=child.get(pos),
        )
        for pos in positions
    ]
    n = len(positions)
    if n < 2:  # 不足两个位点：没有相邻区间，切换次数为 0
        return TrioResult(sites=sites, intervals=[], min_switches=0)

    # 前向：到达每位点各状态的最少切换次数
    fwd: list[dict[tuple[int, int], int]] = [{s: 0 for s in allowed[0]}]
    for i in range(1, n):
        cur: dict[tuple[int, int], int] = {}
        for s2 in allowed[i]:
            cur[s2] = min(
                fwd[i - 1][s1] + (s1[0] != s2[0]) + (s1[1] != s2[1])
                for s1 in allowed[i - 1]
            )
        fwd.append(cur)
    optimal = min(fwd[-1].values())

    # 后向：从每位点各状态到末尾的最少切换次数
    bwd: list[dict[tuple[int, int], int]] = [dict() for _ in range(n)]
    bwd[-1] = {s: 0 for s in allowed[-1]}
    for i in range(n - 2, -1, -1):
        for s1 in allowed[i]:
            bwd[i][s1] = min(
                (s1[0] != s2[0]) + (s1[1] != s2[1]) + bwd[i + 1][s2]
                for s2 in allowed[i + 1]
            )

    # 落在某条最优路径上的边 (s1 -> s2) 决定该区间的切换判定
    intervals: list[IntervalVerdict] = []
    for i in range(n - 1):
        mat_flags: set[bool] = set()
        pat_flags: set[bool] = set()
        for s1 in allowed[i]:
            for s2 in allowed[i + 1]:
                cost = (s1[0] != s2[0]) + (s1[1] != s2[1])
                if fwd[i][s1] + cost + bwd[i + 1][s2] == optimal:
                    mat_flags.add(s1[0] != s2[0])
                    pat_flags.add(s1[1] != s2[1])
        intervals.append(
            IntervalVerdict(
                pos_a=positions[i],
                pos_b=positions[i + 1],
                maternal=_verdict(mat_flags),
                paternal=_verdict(pat_flags),
            )
        )

    return TrioResult(sites=sites, intervals=intervals, min_switches=optimal)
