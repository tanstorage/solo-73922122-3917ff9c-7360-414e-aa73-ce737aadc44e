# 变异核对工具集（FastAPI）

同一参考序列下的两个核对入口：

1. **双 VCF 变异核对**：同一样本、两套变异识别程序输出的单样本 VCF，按**归一化后的等位变异**
   列出两侧一致、仅单侧出现与剂量不符的结果，并保留每条结论的 VCF 原始行号，供人工复查。
2. **家系分相核查**：父、母、子三份分相单样本 VCF，在三人等位相容的继承路径中寻找
   **亲本单倍型切换总次数最少**的路径，给出最小次数，并对每段相邻位点区间标注
   母源 / 父源切换必然发生、可能发生或不发生，供复查相位切换区间。

## 双 VCF 变异核对：功能与判定规则

- **输入限制**
  - 单条染色体参考 FASTA，纯 A/C/G/T 序列，不超过 **1 MiB**；
  - 两份单样本 VCF（恰好 10 列、1 个样本），各不超过 **256 KiB**。
- **接受的变异**
  - A/C/G/T 构成的**单碱基替换**（SNV）；
  - 事件长度不超过 **20 bp** 的**纯插入 / 纯缺失**；
  - GT 必须是**非缺失、未分相的二倍体**（`0/1`、`1/1`，多 ALT 可为 `1/2`）。
  - 其他形态（MNV、符号化 ALT、`*`、分相 `|`、缺失 `.`、非二倍体等）一律判失败。
- **归一化**：对照参考序列对 indel 做左对齐（与 `bcftools norm` / `vt normalize` 等价），
  重复区域中位置写法不同的同一插入 / 缺失归为同一等位变异；SNV 位置不变。
- **剂量计算**：多 ALT 记录按 GT 分别计算每个 ALT 的 0/1/2 剂量，剂量为 0 的变异忽略；
  同一文件内同一归一变异的多条记录剂量累加，超过二倍体上限 2 即判定为**互相矛盾**。
- **整次失败**（返回 422，只显示错误、不展示任何局部结论）：
  - 任一记录 REF 与参考序列不符（含坐标越界、CHROM 不一致）；
  - 同一文件中同一归一变异存在互相矛盾的剂量；
  - 以及任何不满足上述输入规则的情形。
- **结果分组**（均按归一变异排序，并标注两侧原始行号与原始 `POS:REF>ALT` 写法）：
  两侧一致 / 仅 VCF A 或 B 出现 / 剂量不符。

## 家系分相核查：功能与判定规则

- **输入限制**：与双 VCF 核对相同的大小上限；参考为同一条染色体；
  父、母、子各一份单样本 VCF（恰好 10 列、1 个样本）。
- **接受的记录**
  - **单 ALT 的单碱基替换**（REF、ALT 均为单个 A/C/G/T 碱基，且互不相同）；
  - GT 为**二倍体分相写法**（`0|1`、`1|0`、`1|1`、`0|0`），不接受未分相 `/`、缺失 `.`；
  - 同一文件内所有**杂合记录**（`0|1`、`1|0`）必须带有相同的 **PS**（相位集），
    纯合记录不要求 PS；
  - 每人每个位点（POS）**至多一条记录**。
- **遗传模型**：子代 GT **左拷贝为母源、右拷贝为父源**；某人未记录的位点按参考纯合处理；
  同一位置不同文件的 ALT 可以不同（按各自记录的实际碱基判定）。
- **优化目标**：按坐标在三人等位相容（子代母源等位来自母亲两条单倍型之一、
  父源等位来自父亲两条单倍型之一）的继承路径中，最小化**母源 + 父源单倍型切换总次数**
  （前向-后向动态规划）。
- **区间判定**：对每段相邻位点区间，考察所有切换总数最少的最优路径——
  - **必然发生**：每条最优路径都在该区间切换；
  - **可能发生**：部分最优路径在该区间切换（切换位置无法唯一确定，值得复查）；
  - **不发生**：没有任何最优路径在该区间切换。
- **整次失败**（返回 422，只显示错误）：任一位点子代等位无法来自对应亲本，
  或任何输入违反上述规则。

## 目录结构

```
app/            FastAPI 应用与核心比对逻辑（bio.py / trio.py 仅用标准库）
  main.py       路由与上传校验（/、/compare、/trio、/trio/check 及 JSON 接口）
  bio.py        双 VCF：FASTA/VCF 解析、左对齐归一化、剂量与分组
  trio.py       家系：分相 VCF 解析（GT/PS 校验）、最少切换路径与区间判定
  templates/    上传页与结果页（Jinja2）
tests/          97 个单元 / 端到端测试
samples/        可直接试用的演示数据（demo.fa + 两份比对 VCF + 三份家系 VCF）
scripts/        演示数据生成器（坐标由参考自动计算，生成后自检）
Dockerfile, docker-compose.yml, requirements.txt
```

## 用 Docker Compose 启动（推荐）

```bash
docker compose up --build
```

浏览器入口：

- 双 VCF 核对：<http://localhost:8000/>
- 家系分相核查：<http://localhost:8000/trio>
- 交互式 API 文档（Swagger UI）：<http://localhost:8000/docs>

停止：`docker compose down`。容器以只读根文件系统、`no-new-privileges` 运行，
内置 `/` 健康检查。

## 本地运行（无 Docker 时）

```bash
python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

## 试用

**双 VCF 核对**：打开 <http://localhost:8000/>，分别选择
`samples/demo.fa`、`samples/caller_a.vcf`、`samples/caller_b.vcf`。
演示数据覆盖：SNV 一致、重复区内插入 / 缺失的两种等价写法归一、多 ALT 按 GT 拆分、
`0/0` 剂量 0 被忽略、剂量不符（1/1 vs 0/1）、仅单侧出现。

**家系分相核查**：打开 <http://localhost:8000/trio>，分别选择
`samples/demo.fa`、`samples/trio_father.vcf`、`samples/trio_mother.vcf`、`samples/trio_child.vcf`。
演示数据覆盖：最小切换总次数 3（母源 2 次 + 父源 1 次）、
纯合位点导致的切换位置不确定（可能发生）、相邻杂合强制的切换（必然发生）、
以及不发生切换的区间。

也可以走 JSON 接口：

```bash
# 双 VCF 核对
curl -F 'fasta=@samples/demo.fa' \
     -F 'vcf_a=@samples/caller_a.vcf' \
     -F 'vcf_b=@samples/caller_b.vcf' \
     http://localhost:8000/api/compare

# 家系分相核查
curl -F 'fasta=@samples/demo.fa' \
     -F 'vcf_father=@samples/trio_father.vcf' \
     -F 'vcf_mother=@samples/trio_mother.vcf' \
     -F 'vcf_child=@samples/trio_child.vcf' \
     http://localhost:8000/api/trio/check
```

## 测试

```bash
# 核心逻辑（仅标准库）
python3 -I -m unittest discover -s tests -t . -p 'test_bio.py' -v
python3 -I -m unittest discover -s tests -t . -p 'test_trio.py' -v

# 全部测试（Web 层需要先安装 requirements.txt）
python3 -I -m unittest discover -s tests -t . -v
```
