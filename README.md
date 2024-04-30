# `easyTCGA`：让初学者也能感受"征服"`TCGA`的喜悦

## 为什么要写这个R包

生信数据挖掘必不可少要学习`TCGA`数据库，但是对于新手，经常卡在第一步：**下载和整理数据**。第一步完成了，又会卡在第二步，第三步：差异分析，生存分析......

有人会说[XENA](https://gdc.xenahubs.net)有整理好的数据，但这些数据下载后并不能直接用，还是要整理，初学者依然会卡在第一步！

对于R语言大神来说都不是问题，非常简单的R语言操作而已。但是对于初学者很难理解。

这几步操作又是必不可少的，我自己也经常需要重新下载整理数据。为了简化这几个流程，同时也是**让初学者也能感受到"征服"`TCGA`的喜悦**，我把自己常用的一些代码打包，写了这个R包。

## 使用注意

需要自己解决网络问题，比如访问`github，TCGA官网, google`等，如果你无法解决网络问题，那么生信数据挖掘可能不适合你......基本上你常见的生信数据库资源都是国外的，由于众所周知的原因，国外的数据很难下载，**网络问题我帮不了你**。

## 安装

首先安装依赖包：

```R
# 安装bioconductor上面的R包
# 首先要改镜像，下面是清华的镜像，有时会有问题，可更改其他镜像试试（自己百度下喽~）
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
if(!require("BiocManager")) install.packages("BiocManager")
if(!require("TCGAbiolinks")) BiocManager::install("TCGAbiolinks")
if(!require("SummarizedExperiment")) BiocManager::install("SummarizedExperiment")
if(!require("DESeq2")) BiocManager::install("DESeq2")
if(!require("edgeR")) BiocManager::install("edgeR")
if(!require("limma")) BiocManager::install("limma")

# 安装cran上面的R包
if(!require("survival")) install.packages("survival")
if(!require("broom")) install.packages("broom")
if(!require("devtools")) install.packages("devtools")
if(!require("reshape2")) install.packages("reshape2")
if(!require("data.table")) install.packages("data.table")
if(!require("ggplot2")) install.packages("ggplot2")
if(!require("ggpubr")) install.packages("ggpubr")
```

再安装`easyTCGA`包：

```R
devtools::install_github("ayueme/easyTCGA")
```

## 主要功能

> 解决TCGA（GTEx）数据下载和整理问题，顺便实现一些常见的分析和可视化

- `getmrnaexpr`
  - 只需要提供正确的`TCGA project`名字即可；
  - 自动下载并整理`mRNA`和`lncRNA`的`counts，tpm，fpkm`共6种表达矩阵（直接从官网的原始数据提取，未进行任何修改，所以是没有经过log转换的），以及对应的临床信息，临床信息样本顺序和表达矩阵样本顺序完全一致，无需再次整理；
  - 自动保存以上6种表达矩阵和临床信息到当前工作目录下的`output_mRNA_lncRNA_expr`文件夹下，并且同时保存`rdata`和`csv`两种文件格式；
  - 下载的数据为最新数据，和`GDC TCGA`[官网](https://portal.gdc.cancer.gov/)保持一致；
  - 支持通过手动下载的TCGA数据进行自动整理并完成以上过程（可参考b站教程：[easyTCGA：1行代码整理TCGA的6种表达矩阵和临床信息](https://www.bilibili.com/video/BV15c411J7bJ/?share_source=copy_web&vd_source=abc21f68a9e2a784892483fd768dbafa)）
  - lncRNA鉴别参考：[Biotypes (ensembl.org)](http://useast.ensembl.org/info/genome/genebuild/biotypes.html)
- `getmrnaexpr_xena`
  - 用于`XENA`网站下载的TCGA基因表达数据和临床信息的整理（仅限`gdchub`）；
  - 直接提供文件名即可，比如：`TCGA-ACC.htseq_counts.tsv.gz, TCGA-ACC.htseq_fpkm.tsv.gz`，`TCGA-ACC.GDC_phenotype.tsv.gz, TCGA-ACC.survival.tsv`；
  - 自动保存`mRNA`、`lncRNA`表达矩阵和临床信息到当前工作目录下的`output_mRNA_expr_xena`文件夹下；
  - id转换使用`gtf 22`，和`XENA`保持一致；
  - （单独使用XENA的表达谱数据和直接用GDC官网数据相比没有任何优势）
- `getmirnaexpr`
  - 只需要提供正确的`TCGA project`名字即可；
  - 自动下载并整理`miRNA`的`counts，rpm`2种表达矩阵；
  - 自动保存以上2种表达矩阵和对应的临床信息到当前工作目录下的`output_miRNA_expr`文件夹下，并且同时保存`rdata`和`csv`两种文件格式；
  - 下载的数据为最新数据，和`GDC TCGA`[官网](https://portal.gdc.cancer.gov/)保持一致
- `getsnvmaf`
  - 只需要提供正确的`TCGA project`名字即可；
  - 自动下载并整理`TCGA MAF`文件(masked somatic mutation)以及对应的临床信息，并自动保存到当前工作目录下的`output_snv`文件夹下；
  - 输出结果可以直接通过`maftools::read_maf()`函数读取，无需再次整理
- `getcnv`
  - 只需要提供正确的`TCGA project`名字即可；
  - 自动下载并整理`copy number variation`数据；数据保存到当前工作目录下的`output_cnv`文件夹下；
  - 下载的数据为最新数据，和`GDC TCGA`[官网](https://portal.gdc.cancer.gov/)保持一致
- `getmethybeta`
  - 只需要提供正确的`TCGA project`名字即可；
  - 自动下载并整理`450K`的`DNA methylation`的`beta值矩阵`，以及对应的临床信息，数量和顺序完全一致，无需再次整理；
  - 自动整理探针信息，比如探针对应的`gene symbol`等，基于`GRCh 38`；
  - 数据保存在当前工作目录下的`output_methy`文件夹下；
  - 下载的数据为最新数据，和`GDC TCGA`[官网](https://portal.gdc.cancer.gov/)保持一致
- `getclinical`
  - 下载XML格式的临床数据，包括各种常见的临床信息，如生存信息、病理分期、放化疗数据、化疗药物数据等
  - 与GDC TCGA[官网](https://portal.gdc.cancer.gov/)数据保持一致
  - 只需要提供正确的`TCGA project`名字即可

- `getpancancer_xena`
  - 实现对泛癌数据的整理，支持`TCGA`、`GTEx`，以及整合`TCGA+GTEx`
  - 原始文件是从XENA下载的；
  - 只需提供相应的表达矩阵文件和样本信息文件即可
  - 很费内存，可在公众号后台直接回复**pancancer**获取我整理好的

- `diff_analysis`
  - 与`getmrnaexpr`，`getmirnaexpr`，`getmrnaexpr_xena`函数无缝对接，直接使用其输出结果即可，无需任何整理（默认对tumor和normal组进行差异分析）；
  - 支持`count, tpm, fpkm`和`GEO`数据，如果是`count`则自动通过3个R包进行差异分析：`DESeq2, edgeR, limma`；如果是其他类型（`tpm, fpkm`和`基因表达芯片数据`）会自动判断是否需要`log2(x + 0.1)`转换，然后使用`limma`和`wilcoxon test`做差异分析；
  - 用`wilcoxon`秩和检验做差异分析的参考资料：[TCGA等大样本量差异分析该使用DEseq2还是edgeR呢？](https://mp.weixin.qq.com/s/PlQ9Sl6B12k9tSopSHQfyw),以及文中涉及的参考文献：https://doi.org/10.1186/s13059-022-02648-4
  - 支持输入自己的表达矩阵和**自定义分组**，分组信息需要因子型向量；
  - 输出结果默认为1个`list`，内含多种差异分析结果，支持保存`rdata`格式数据到本地
- `batch_survival`
  - 自动进行`logrank`检验和单因素`cox`分析，默认基于**最佳截点（P值最小）**；
  - 与`getmrnaexpr`，`getmirnaexpr`函数无缝对接，直接使用其输出结果即可，无需任何整理；
  - 支持`count，tpm，fpkm`3种格式的数据，如果是`counts`，则通过`DESeq2::vst()`进行转换，如果是`tpm/fpkm`，则进行`log2(x + 0.1)`转换；
  - 支持打印基因序号到屏幕，方便定位有问题的基因
- 可视化函数
  - 主要用来进行一些简单的探索；每个函数都会返回画图数据，方便你自己探索；
  - `plot_gene`：任意数量基因在任意癌种（TCGA33种其中之一都可以）的任意分组中的表达量箱线图；
  - `plot_gene_paired`：任意基因在某一癌种配对样本中的表达量配对箱线图；
  - `plot_km`：根据任意基因的表达量分组，并画出K-M生存曲线（支持最佳截点）


## 使用教程

文字版使用教程请关注公众号：**医学和生信笔记**。

- [easyTCGA：1行代码搞定TCGA的6种表达矩阵和临床信息](https://mp.weixin.qq.com/s/z1fgyXLZXwmoaI39f2ftYw)
- [easyTCGA：1行代码搞定TCGA突变maf文件下载和整理](https://mp.weixin.qq.com/s/GBkB8Hv45l06BVnyFNFzzw)
- [easyTCGA生存分析支持最佳截点，任意基因在不同组中的表达量箱线图](https://mp.weixin.qq.com/s/Qc9m6hX-qKVJt5GzrXY9bA)
- [TCGA、GTEx的泛癌数据也是1行代码整理](https://mp.weixin.qq.com/s/SzGB1wVH_DNBbXxvkBe5NA)
- [任意基因在泛癌中的表达量展示](https://mp.weixin.qq.com/s/MIDRG57oRSMTyX6Gm99-3w)
- [TCGA临床数据（化疗数据、用药反应等）和生存信息（4种临床结局）整理](https://mp.weixin.qq.com/s/FSdQ3UbsGSzop2ETVNkwkA)

视频版教程请关注哔哩哔哩：[阿越就是我](https://space.bilibili.com/42460432)，（视频教程滞后于包的更新速度）

- [easyTCGA：让初学者也能感受“征服”TCGA的喜悦](https://www.bilibili.com/video/BV1Jm4y1y7Qb/?share_source=copy_web&vd_source=abc21f68a9e2a784892483fd768dbafa)
- [easyTCGA：1行代码整理TCGA的6种表达矩阵和临床信息](https://www.bilibili.com/video/BV15c411J7bJ/?share_source=copy_web&vd_source=abc21f68a9e2a784892483fd768dbafa)
- [easyTCGA：1行代码整理TCGA的miRNA和突变数据](https://www.bilibili.com/video/BV1mP411m7Rc/?share_source=copy_web&vd_source=abc21f68a9e2a784892483fd768dbafa)
- [easyTCGA：1行代码完成TCGA的3种差异分析](https://www.bilibili.com/video/BV1Ao4y177TY/?share_source=copy_web&vd_source=abc21f68a9e2a784892483fd768dbafa)
- [easyTCGA：1行代码搞定两种TCGA批量生存分析](https://www.bilibili.com/video/BV18X4y1B7nB/?share_source=copy_web&vd_source=abc21f68a9e2a784892483fd768dbafa)

## 问题反馈

B站，公众号，Github，粉丝QQ群，都可以。

## TO DO

- [x] 支持`XENA`网站下载的`gene expression`和临床数据的整理
- [x] 支持`XENA`泛癌数据的整理，对电脑内存要求较高，正在优化代码中......
- [x] 增加对`miRNA`的差异分析支持
- [x] 增加对`miRNA`的批量生存分析支持
- [x] 增加对自定义表达矩阵/自定义分组差异分析的支持
- [ ] 增加对多分组差异分析的支持
- [x] 增加对`lncRNA`的差异分析和批量生存分析支持
- [x] 实现一些常见的分析和可视化
- [ ] 支持自定义生存信息
- [ ] ......

