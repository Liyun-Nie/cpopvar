# cpopvar

[English](README.md)

`cpopvar` 用于叶绿体基因组（plastome）的群体变异分析。完成变异检测后，可以用它完成预处理、按样本出现频次过滤、每千碱基归一化，以及变异分布、分组比较、CDS / 基因间隔区（IGS）变异热点和染色体示意图分析。

目前通过 R 函数或 YAML 配置运行，没有图形界面。

cpopvar 目前是预发布版，尚未提供 CRAN 安装。

## 可以回答哪些问题

先完成共享的预处理与归一化，再按问题启用后续任务。不必一次全开；一次运行的结果可以用来收窄下一次的输入、分组或候选位点。启用任一 M01–M05 任务时，P01–P03 会一并运行。

- **处理基础。** P01 给变异加注释并整理区域信息；P02 按样本出现频次阈值过滤变异；P03 在统一的区域类型、变异类型和长度下做每千碱基归一化。
- **变异分布在哪里，是否不均匀？** 启用 `M01_distribution`。
- **哪些生物学分组存在差异？** 启用 `M02_comparative`。
- **哪些 CDS 基因或基因间隔区（IGS）是高变异候选？** 启用 `M03_hotspot` 和/或 `M04_igs_hotspot`。
- **这些候选位点落在叶绿体基因组图谱的何处？** 启用 `M05_hotspot_ideogram`（同一运行中需同时启用 `M03_hotspot` 和 `M04_igs_hotspot`）。

## 使用流程

1. 请准备 CSV/TSV 输入表（见下文「输入数据」）。
2. 按研究问题选择要启用的 M01–M05 任务。
3. 用 `cpopvar_config()` 在内存中创建配置，或准备内容等价的 YAML 文件。
4. 运行前请先用 `validate_cpopvar_config()` 检查配置，通过后再调用 `run_analysis()`。
5. 运行后查看返回的会话目录；如需调整输入或任务，改配置后再跑一次。

## 预发布版安装

cpopvar 目前是预发布版，尚未提供 CRAN 安装。请从本地源码压缩包安装。将下面命令中的 `VERSION` 换成归档文件名里的版本号。

**干净 R 库推荐。** 先安装 `remotes`，再安装压缩包并解析依赖：DESCRIPTION 中声明的硬依赖（`Depends`、`Imports`、`LinkingTo`）会从当前已配置的软件源解析或安装。

```r
install.packages("remotes")
remotes::install_local("cpopvar_VERSION.tar.gz", dependencies = NA, upgrade = "never")
library(cpopvar)
```

`dependencies = NA` 会安装上述硬依赖，不会安装 `Suggests` 中的包。`upgrade = "never"` 表示不升级已经装好的包。

**仅用基础 R 的备选方式。** 仅在所需依赖已经安装时使用：

```r
install.packages("cpopvar_VERSION.tar.gz", repos = NULL, type = "source")
library(cpopvar)
```

`repos = NULL` 不会从 CRAN 解析或安装缺失依赖。在干净的机器上，如果这些包尚未安装，该命令可能失败。

## 环境与依赖

- 需要 R 4.1.0 或更高版本（`Depends: R (>= 4.1.0)`）。
- 核心 R 包列在 DESCRIPTION 的 `Imports` 中。推荐的 `remotes::install_local(..., dependencies = NA)` 会安装或解析这些硬依赖。cpopvar 本身不会在运行时安装软件包。
- Snippy 是用于准备 `main_data` 的外部上游工具，cpopvar 既不安装也不调用它。
- 在 Windows 上，只有当 R 需要从源码编译某个依赖时才可能需要 Rtools；并非每次安装都必需。
- `Suggests` 中的包按功能选用，`dependencies = NA` 不会安装它们。缺少某个可选包时，请只为要用的功能单独安装。

## 快速开始

包内带有一套小型合成示例，可用来走通一次分析。下面这段代码会读取随包示例、构建配置，并在校验通过后运行 CDS 变异热点分析（`M03_hotspot`）：

```r
example_dir <- system.file("extdata", "synthetic", package = "cpopvar")
bundled <- yaml::read_yaml(file.path(example_dir, "config.yml"))
bundled$input_files <- lapply(
  bundled$input_files,
  function(path) file.path(example_dir, path)
)

cfg <- cpopvar::cpopvar_config(
  main_data = bundled$input_files$main_data,
  annotations = bundled$input_files$annotations,
  genome_regions = bundled$input_files$genome_regions,
  group_info = bundled$input_files$group_info,
  species_order = bundled$input_files$species_order,
  gene_function_mapping = bundled$input_files$gene_function_mapping,
  tasks = "M03_hotspot",
  task_options = list(M03_hotspot = bundled$analysis_tasks$M03_hotspot$parameters),
  visualization_settings = bundled$visualization_settings,
  output_settings = bundled$output_settings
)
stopifnot(cpopvar::validate_cpopvar_config(cfg)$valid)

result_dir <- cpopvar::run_analysis(
  cfg,
  output_dir = file.path(tempdir(), "cpopvar-example")
)
```

`result_dir` 是本次运行写入的会话目录，位于 `output_dir` 下的 `app_data/sessions/<session_id>/`。打开该目录即可查看输入副本、处理表、图和日志。

也可以把 YAML 配置文件路径传给 `run_analysis()`，参数名为 `config` 或 `config_file`。YAML 里的相对路径需要先改成绝对路径，或在该 YAML 所在目录下运行。

## 输入数据

当前 `cpopvar` 从已整理的 CSV/TSV 表开始，不直接读取组装 FASTA 或 GenBank 文件。

核心表的外部列名可在 YAML 的 `column_mappings`（或 `cpopvar_config()` 的对应参数）中映射到内部名称。分组、排序等辅助文件保留配置中的列名，由后续任务直接使用。

随包示例见 [`inst/extdata/synthetic/`](inst/extdata/synthetic/config.yml)。

### 需要准备的表格

- **`main_data`（必需）。** 每行是一条样本水平的变异记录；同一个样本通常有多行。推荐列名为：

  `sample_id`, `species`, `position`, `var_type`, `gene`, `region_type`,
  `ref_allele`, `alt_allele`

  其他表头可在配置中映射。示例：
  [`inst/extdata/synthetic/main_data.csv`](inst/extdata/synthetic/main_data.csv)。

- **`annotations`（完整 P01–P03 所需）。** 当前输入是 CSV/TSV，不是 `.gb`。需要字段 `Species`、`Type`、`Gene`、`Minimum`、`Maximum`、`Length`、`Intervals`。若 Geneious 导出使用 `Sequence Name` / `# Intervals`，请先改成上述表头。示例：
  [`inst/extdata/synthetic/annotations.csv`](inst/extdata/synthetic/annotations.csv)。

- **`genome_regions`（完整 P01–P03 所需）。** 每个物种一行，给出 LSC、倒位重复区（IR）、SSC 的起止坐标：

  `species`, `lsc_start`, `lsc_end`, `ir_start`, `ir_end`, `ssc_start`,
  `ssc_end`

  有基因组总长时请提供 `total_length`。缺少倒位重复区的基因组可在 `special_handling` 中标记为 `IR_lacking_genome`。示例：
  [`inst/extdata/synthetic/genome_regions.csv`](inst/extdata/synthetic/genome_regions.csv)。

- **可选辅助表。** 请对照随包示例准备，不要假定统一字段：
  - [`group_info.csv`](inst/extdata/synthetic/group_info.csv) — 物种水平分组标签，供 `M02_comparative` 使用。
  - [`species_order.csv`](inst/extdata/synthetic/species_order.csv) — 显示顺序与可选绘图标签。
  - [`gene_function_map.csv`](inst/extdata/synthetic/gene_function_map.csv) — 基因功能查找表，供变异热点富集视图使用。

YAML 示例：
[`inst/extdata/synthetic/config.yml`](inst/extdata/synthetic/config.yml)。

### 如何获得 main_data

以下步骤在 `cpopvar` 之外完成。当前包不读取 FASTA 或 GenBank，也不调用 Snippy。

对每个物种，请准备：该物种全部样本的组装叶绿体基因组 FASTA，以及一份带详细注释的参考 GenBank（通常来自其中一个样本）。组装序列与这份参考都需要一致地去掉一份倒位重复区（IR）拷贝，避免同一变异被计两次，并保持坐标一致。

然后对每个样本的组装，用该物种参考运行 Snippy，例如：

```bash
snippy --outdir SAMPLE --ctgs SAMPLE.single_IR.fasta --ref SPECIES_REFERENCE.gb --cpus 1
```

在每个样本的 `snps.tab` 前加上物种和样本标识；先在物种内拼接，再合并各物种表，并把列名映射为上面的 `main_data` 规范列。直接拼接的 Snippy 中间表还不能作为输入，请提供映射后的表。

## 输出结果

`run_analysis()` 把结果写到会话目录：

`<output_dir>/app_data/sessions/<session_id>/`

具体文件取决于本次启用的任务。常见布局如下。

- `raw/` — 本次使用的输入表副本（注释表在 `raw/annotations/`）。
- `results/processed_data/P01_preprocessed/` — 已注释的变异表（全基因组表与仅保留一个 IR 拷贝的表）以及 `region_info_complete.csv`。
- `results/processed_data/P02_filtered/` — `variants_filtered.csv`（按样本出现频次阈值过滤后的变异）。
- `results/processed_data/P03_normalized/` — `normalized_frequencies.csv` 和 `normalized_frequencies_summary.csv`。
- `results/plots/` — 已启用的 M01–M05 表格与图。可能包括 M01 分布摘要与图、M02 比较摘要与图、M03 CDS 变异热点候选、M04 基因间隔区（IGS）变异热点候选，以及 M05 染色体示意图与标记/热图表。
- `results/run_summary.rds` — 已执行任务摘要。若启用报告，还会写出 `results/Final_Analysis_Report.md`。本次实际使用的配置保存在 `results/config_used.yml`。
- `logs/` — 会话日志。

## 引用

安装后使用 `citation("cpopvar")`。文章元数据将在 DOI 可用时添加。

## 许可证

MIT © 2026 Liyun Nie。见 `LICENSE` 和 `LICENSE.md`。
