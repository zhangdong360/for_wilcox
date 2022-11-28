# wilcox自动循环

本脚本为方便进行wilcox批量循环所写。

## 使用方法

使用`source()`函数导入此脚本，如放置在项目下的bio目录，则运行`source("bio/for_wilcox.R")`即可载入脚本直接使用。

也可以在Rstudio中打开此脚本，运行脚本内的`function()`函数进行载入。

# 数据要求

### data数据要求

输入数据的示例数据为`data1`，`data2`。其中数据格式为，行为**gene**，列为**sample**，请保证数据框内没有**非numeric数据**。

### group数据要求

输入分组示例数据为`group_AHR`，其所需要的两列分别为**sample_ID**和**group**列，其中**group**列是必需的，**sample_ID**列为建议。

**sample_ID**存放sample_ID，**group**存放对应的分组信息。

`for_wilcox()`函数将通过`rawname(group)`来读取sample_ID，通过`group[,by]`来读取group信息。

其中，`group[,by]`的参数**by**列名可使用`for_wilcox()`中的**by**参数进行选择。

by参数的默认值为*group*。

需注意：

+ sample_ID需与`data1`和`data2`的colnames对应。
+ group列建议转化为factor格式，并检查level的顺序，使得实验组在前，对照组在后。转化factor的函数为`factor(group_AHR$type,levels = c("AHR","NO_AHR"))`



### 其他

初期版本并未添加WARNING和ERROR提示，请严格按照数据要求使用。

## 测试代码

测试代码：

```R
test <- for_wilcox(data_1,data_2,
                   group = group_AHR,
                   by = 'type',
                   dir = "result/AHR/")
```

## 结果说明

`for_wilcox()`函数运行后会返回一个list结果，其数据结构为：

```R
  result <- list(all = result_all,
                 fill_pvalue = subset(result_all,result_all$pvalue < 0.05),
                 up = subset(result_all,result_all$pvalue < 0.05&result_all$logFC > 0),
                 down = subset(result_all,result_all$pvalue < 0.05&result_all$logFC < 0))
```

其中

+ `result$all`存放的是全部的差异分析结果
+ `result$fill_pvalue`存放的是根据p value < 0.05进行筛选后的结果
+ `result$up`存放的是p value < 0.05且logFC > 0的结果
+ `result$down`存放的是p value < 0.05且logFC < 0的结果

