# mtfind 0.311

mtfind 根据长读长 RNA 测序中同一条 read 的 5′ 和 3′ 端点，发现线粒体 RNA 坐标候选。它以比对后生成的 SSC 文件为输入，候选发现阶段不读取注释。

输出的 GTF 是坐标候选展示轨道，不是经过验证的转录本注释。

English version: [README.md](README.md).

## 安装

需要 Python 3.10 或更高版本。

```bash
python3 -m pip install -r requirements.txt
```

## 一键示例

本仓库提供一个来自公开 SG-NEx A549 direct-RNA 样本的 MT-only SSC 派生输入。运行：

```bash
bash examples/sgnex_a549/run_demo.sh
```

脚本复用 mtfind 0.311 对该 A549 输入已记录的参数，并在 `examples/sgnex_a549/output/` 生成结果；若该目录已存在，脚本会拒绝覆盖。输出包括候选 GTF、审计表、运行参数和输入/代码副本。输入来源和引用要求见 [DATA_PROVENANCE_zh.md](examples/sgnex_a549/DATA_PROVENANCE_zh.md)。

## 内容

- `code/annotation_free_adaptive_depth.py`：mtfind 0.311 候选发现主程序。
- `code/gtf_coordinate_benchmark.py`：可选的运行后坐标比较模块，按实际执行版本原样保留。
- `tests/test_seed_backfill.py`：seed/backfill 聚类的小型回归测试。
- `examples/sgnex_a549/`：一个真实数据 SSC 示例及其准确运行命令。

## 快速检查

```bash
python3 code/annotation_free_adaptive_depth.py --version
python3 -m unittest discover -s tests -v
```

本包不含许可证文件；公开前请在 GitHub 仓库中选择许可证。
