"""FastAPI 入口：

- `/` 与 `/compare`：双 VCF 变异核对（参考 FASTA + 两份单样本 VCF）；
- `/trio` 与 `/trio/check`：家系分相核查（参考 FASTA + 父、母、子三份分相单样本 VCF），
  寻找亲本单倍型切换总次数最少的继承路径，并标注每段相邻位点区间的切换必然性。
"""

from pathlib import Path

from fastapi import FastAPI, File, Request, UploadFile
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.templating import Jinja2Templates

from .bio import (
    MAX_FASTA_BYTES,
    MAX_VCF_BYTES,
    CheckError,
    compare,
    parse_fasta,
    parse_vcf,
)
from .trio import MUST, MAYBE, analyze_trio, parse_trio_vcf

app = FastAPI(
    title="变异核对工具集",
    description="双 VCF 变异核对 + 家系分相核查：上传单条染色体参考 FASTA 与单样本 VCF，"
    "核对两套识别结果差异，或寻找亲本单倍型切换最少的继承路径。",
)

templates = Jinja2Templates(directory=str(Path(__file__).parent / "templates"))


def _read_limited(file: UploadFile, limit: int, label: str) -> bytes:
    data = file.file.read(limit + 1)
    if len(data) > limit:
        raise CheckError(
            f"{label} 超过大小上限 {limit} 字节"
            f"（FASTA 上限 {MAX_FASTA_BYTES // 1024} KiB，VCF 上限 {MAX_VCF_BYTES // 1024} KiB）。"
        )
    if not data:
        raise CheckError(f"{label} 为空文件。")
    return data


def _run_compare(fasta_bytes: bytes, vcf_a_bytes: bytes, vcf_b_bytes: bytes):
    refobj = parse_fasta(fasta_bytes)
    calls_a, sample_a = parse_vcf(vcf_a_bytes, "VCF A", refobj)
    calls_b, sample_b = parse_vcf(vcf_b_bytes, "VCF B", refobj)
    groups = compare(refobj, calls_a, calls_b)
    return refobj, sample_a, sample_b, groups


def _run_trio(fasta_bytes: bytes, father_bytes: bytes, mother_bytes: bytes, child_bytes: bytes):
    refobj = parse_fasta(fasta_bytes)
    father, sample_father = parse_trio_vcf(father_bytes, "父本 VCF", refobj)
    mother, sample_mother = parse_trio_vcf(mother_bytes, "母本 VCF", refobj)
    child, sample_child = parse_trio_vcf(child_bytes, "子代 VCF", refobj)
    result = analyze_trio(refobj, father, mother, child)
    return refobj, sample_father, sample_mother, sample_child, result


def _trio_form_context(**extra):
    context = {
        "fasta_limit_kib": MAX_FASTA_BYTES // 1024,
        "vcf_limit_kib": MAX_VCF_BYTES // 1024,
    }
    context.update(extra)
    return context


@app.get("/", response_class=HTMLResponse)
def index(request: Request):
    return templates.TemplateResponse(request, "index.html", _trio_form_context())


@app.post("/compare", response_class=HTMLResponse)
def compare_page(
    request: Request,
    fasta: UploadFile = File(...),
    vcf_a: UploadFile = File(...),
    vcf_b: UploadFile = File(...),
):
    try:
        fasta_bytes = _read_limited(fasta, MAX_FASTA_BYTES, "参考 FASTA")
        vcf_a_bytes = _read_limited(vcf_a, MAX_VCF_BYTES, "VCF A")
        vcf_b_bytes = _read_limited(vcf_b, MAX_VCF_BYTES, "VCF B")
        refobj, sample_a, sample_b, groups = _run_compare(
            fasta_bytes, vcf_a_bytes, vcf_b_bytes
        )
    except CheckError as exc:
        # 整次核对失败：只展示错误，不展示任何局部结论
        return templates.TemplateResponse(
            request,
            "index.html",
            _trio_form_context(error=str(exc)),
            status_code=422,
        )

    return templates.TemplateResponse(
        request,
        "result.html",
        {
            "ref_name": refobj.name,
            "ref_length": len(refobj.sequence),
            "sample_a": sample_a,
            "sample_b": sample_b,
            "groups": groups,
        },
    )


@app.post("/api/compare")
def compare_api(
    fasta: UploadFile = File(...),
    vcf_a: UploadFile = File(...),
    vcf_b: UploadFile = File(...),
):
    """JSON 接口，便于脚本化调用；失败返回 422 与错误信息。"""
    try:
        fasta_bytes = _read_limited(fasta, MAX_FASTA_BYTES, "参考 FASTA")
        vcf_a_bytes = _read_limited(vcf_a, MAX_VCF_BYTES, "VCF A")
        vcf_b_bytes = _read_limited(vcf_b, MAX_VCF_BYTES, "VCF B")
        refobj, sample_a, sample_b, groups = _run_compare(
            fasta_bytes, vcf_a_bytes, vcf_b_bytes
        )
    except CheckError as exc:
        return JSONResponse(status_code=422, content={"ok": False, "error": str(exc)})

    def dump_entry(entry):
        def dump_call(call):
            if call is None:
                return None
            return {
                "dosage": call.dosage,
                "line_numbers": call.line_numbers,
                "raw": call.raw_descriptions,
            }

        v = entry["variant"]
        return {
            "variant": {
                "pos": v.pos,
                "ref": v.ref,
                "alt": v.alt,
                "kind": v.kind,
                "label": v.label,
            },
            "left": dump_call(entry["left"]),
            "right": dump_call(entry["right"]),
        }

    return {
        "ok": True,
        "reference": {"name": refobj.name, "length": len(refobj.sequence)},
        "samples": {"left": sample_a, "right": sample_b},
        "result": {key: [dump_entry(e) for e in rows] for key, rows in groups.items()},
    }


# --------------------------------------------------------------------------- #
# 家系分相核查
# --------------------------------------------------------------------------- #
@app.get("/trio", response_class=HTMLResponse)
def trio_index(request: Request):
    return templates.TemplateResponse(request, "trio.html", _trio_form_context())


def _trio_result_context(refobj, sample_father, sample_mother, sample_child, result):
    intervals = result.intervals
    stats = {
        "m_must": sum(1 for iv in intervals if iv.maternal == MUST),
        "m_maybe": sum(1 for iv in intervals if iv.maternal == MAYBE),
        "p_must": sum(1 for iv in intervals if iv.paternal == MUST),
        "p_maybe": sum(1 for iv in intervals if iv.paternal == MAYBE),
    }
    return {
        "ref_name": refobj.name,
        "ref_length": len(refobj.sequence),
        "sample_father": sample_father,
        "sample_mother": sample_mother,
        "sample_child": sample_child,
        "result": result,
        "stats": stats,
    }


@app.post("/trio/check", response_class=HTMLResponse)
def trio_check_page(
    request: Request,
    fasta: UploadFile = File(...),
    vcf_father: UploadFile = File(...),
    vcf_mother: UploadFile = File(...),
    vcf_child: UploadFile = File(...),
):
    try:
        fasta_bytes = _read_limited(fasta, MAX_FASTA_BYTES, "参考 FASTA")
        father_bytes = _read_limited(vcf_father, MAX_VCF_BYTES, "父本 VCF")
        mother_bytes = _read_limited(vcf_mother, MAX_VCF_BYTES, "母本 VCF")
        child_bytes = _read_limited(vcf_child, MAX_VCF_BYTES, "子代 VCF")
        refobj, sample_father, sample_mother, sample_child, result = _run_trio(
            fasta_bytes, father_bytes, mother_bytes, child_bytes
        )
    except CheckError as exc:
        # 整次核查失败：只展示错误，不展示任何局部结论
        return templates.TemplateResponse(
            request,
            "trio.html",
            _trio_form_context(error=str(exc)),
            status_code=422,
        )

    return templates.TemplateResponse(
        request,
        "trio_result.html",
        _trio_result_context(refobj, sample_father, sample_mother, sample_child, result),
    )


@app.post("/api/trio/check")
def trio_check_api(
    fasta: UploadFile = File(...),
    vcf_father: UploadFile = File(...),
    vcf_mother: UploadFile = File(...),
    vcf_child: UploadFile = File(...),
):
    """家系分相核查的 JSON 接口；失败返回 422 与错误信息。"""
    try:
        fasta_bytes = _read_limited(fasta, MAX_FASTA_BYTES, "参考 FASTA")
        father_bytes = _read_limited(vcf_father, MAX_VCF_BYTES, "父本 VCF")
        mother_bytes = _read_limited(vcf_mother, MAX_VCF_BYTES, "母本 VCF")
        child_bytes = _read_limited(vcf_child, MAX_VCF_BYTES, "子代 VCF")
        refobj, sample_father, sample_mother, sample_child, result = _run_trio(
            fasta_bytes, father_bytes, mother_bytes, child_bytes
        )
    except CheckError as exc:
        return JSONResponse(status_code=422, content={"ok": False, "error": str(exc)})

    def dump_rec(rec):
        if rec is None:
            return None
        return {
            "pos": rec.pos,
            "ref": rec.ref,
            "alt": rec.alt,
            "gt": rec.gt,
            "ps": rec.ps,
            "line_number": rec.line_no,
        }

    return {
        "ok": True,
        "reference": {"name": refobj.name, "length": len(refobj.sequence)},
        "samples": {
            "father": sample_father,
            "mother": sample_mother,
            "child": sample_child,
        },
        "min_switches": result.min_switches,
        "sites": [
            {
                "pos": site.pos,
                "ref": site.ref,
                "father": dump_rec(site.father),
                "mother": dump_rec(site.mother),
                "child": dump_rec(site.child),
            }
            for site in result.sites
        ],
        "intervals": [
            {
                "pos_a": iv.pos_a,
                "pos_b": iv.pos_b,
                "maternal": iv.maternal,
                "paternal": iv.paternal,
            }
            for iv in result.intervals
        ],
    }
