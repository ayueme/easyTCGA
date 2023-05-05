# `easyTCGA`：让初学者也能感受"征服"`TCGA`的喜悦

## 为什么要写这个R包

生信数据挖掘必不可少要学习`TCGA`数据库，但是对于新手，经常卡在第一步：**下载和整理数据**。第一步完成了，又会卡在第二步，第三步：差异分析，生存分析......

对于R语言大神来说都不是问题，非常简单的R语言操作而已。但是对于初学者很难理解。

这几步操作又是必不可少的，我自己也经常需要重新下载整理数据。为了简化这几个流程，同时也是**让初学者也能感受到"征服"`TCGA`的喜悦**，我把自己常用的一些代码打包，写个R包玩玩。

## 安装

首先安装依赖包：

```R
# 安装bioconductor上面的R包
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
if(!require("cli")) install.packages("cli")
```

再安装`easyTCGA`包：

```R
devtools::install_github("ayueme/easyTCGA")
```

## 主要功能

> 1行代码实现1个常见分析！

- `getmrnaexpr`
  - 只需要提供正确的`TCGA project`名字即可；
  - 自动下载并整理`mRNA`和`lncRNA`的`counts，tpm，fpkm`共6种表达矩阵，以及对应的临床信息，临床信息样本顺序和表达矩阵样本顺序完全一致，无需再次整理；
  - 自动保存以上6种表达矩阵和临床信息到当前工作目录下的`output_mRNA_lncRNA_expr`文件夹下，并且同时保存`rdata`和`csv`两种文件格式；
  - 下载的数据为最新数据，和`GDC TCGA`[官网](https://portal.gdc.cancer.gov/)保持一致；
  - 支持通过手动下载的TCGA数据进行自动整理并完成以上过程
- `getmrnaexpr_xena`
  - 用于`XENA`网站下载的基因表达数据和临床信息的整理（`gdchub`）
  - 直接提供文件名即可，比如：`TCGA-ACC.htseq_counts.tsv.gz`，`TCGA-ACC.GDC_phenotype.tsv.gz`
  - 自动保存`mRNA`、`lncRNA`表达矩阵和临床信息到当前工作目录下的`output_mRNA_expr_xena`文件夹下
  - （单独使用和GDC官方数据没有任何优势）

- `getmirnaexpr`
  - 只需要提供正确的`TCGA project`名字即可；
  - 自动下载并整理`miRNA`的`counts，rpm`2种表达矩阵；
  - 自动保存以上2种表达矩阵和对应的临床信息到当前工作目录下的`output_miRNA_expr`文件夹下，并且同时保存`rdata`和`csv`两种文件格式；
  - 下载的数据为最新数据，和`GDC TCGA`[官网](https://portal.gdc.cancer.gov/)保持一致
- `getsnvmaf`
  - 只需要提供正确的`TCGA project`名字即可；
  - 自动下载并整理`TCGA MAF`文件(masked somatic mutation)以及对应的临床信息，并自动保存到当前工作目录下的`output_snv`文件夹下；
  - 输出结果可以直接通过`maftools::read_maf()`函数读取，无需再次整理
- `diff_analysis`
  - 与`getmrnaexpr`和`getmirnaexpr`函数无缝对接，直接使用其输出结果即可（只支持`counts`矩阵），无需任何整理；
  - 支持输入自己的表达矩阵和自定义分组；
  - 自动通过3个R包进行差异分析：`DESeq2, edgeR, limma`；
  - 输出结果默认为1个`list`，内含3种差异分析结果，支持保存`rdata`格式数据到本地
- `batch_survival`
  - 自动对大约20000个基因进行`logrank`检验和单因素`cox`分析，默认基于**最佳截点（P值最小）**；
  - 与`getmrnaexpr`函数无缝对接，直接使用其输出结果即可，无需任何整理；
  - 支持`counts，tpm，fpkm`3种格式的数据，如果是`counts`，则通过`DESeq2::vst()`进行转换，如果是`tpm/fpkm`，则进行`log2(x + 0.1)`转换；
  - 支持打印基因序号到屏幕，方便定位有问题的基因

## 使用教程

文字版使用教程请关注公众号：**医学和生信笔记**

视频版教程请关注哔哩哔哩：[阿越就是我](https://space.bilibili.com/42460432)

## 问题反馈

B站，公众号，Github，粉丝QQ群，都可以。

## 使用注意

需要自己解决网络问题，比如访问`github，TCGA官网, google`等，如果你无法解决网络问题，那么生信数据挖掘可能不适合你......

## TO DO

- [x] 支持`XENA`网站下载的`gene expression`和临床数据的整理
- [ ] 支持`XENA`泛癌数据的整理
- [x] 增加对`miRNA`的差异分析支持
- [x] 增加对`miRNA`的批量生存分析支持
- [x] 增加对自定义表达矩阵/自定义分组差异分析的支持
- [ ] 增加对多分组差异分析的支持
- [x] 增加对`lncRNA`的差异分析和批量生存分析支持
- [ ] ......

