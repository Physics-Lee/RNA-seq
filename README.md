此库为2023春季学期，USTC系统生物学课程的lab1中的RNA-seq。

此处，上游定义为从原始数据获得基因表达矩阵，下游定义为对基因表达矩阵进行各种机器学习的分析。

所谓的基因表达矩阵，用机器学习的话说，一个轴为feature，一个轴为sample。

显然，只有上游涉及生物学，下游完全可以交给机器学习的人干了。

上游的注意事项

* RNA-seq可以测rRNA，也可以测mRNA, tRNA，取决于你想要研究什么。
* 基因表达矩阵的最小单位被称为reads，是一段50-250 base pair的RNA序列。