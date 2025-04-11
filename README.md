# semg-decomposition-newest-method-
2025/2/16 完成算法的可视化与模块化

预处理部分：

基于fieldtrip-master包，10-500Hz的4阶巴特沃斯带通滤波，陷波滤波使用dftFilterGUI。滤波教程参考README_filter.md


分解部分：

Main_dynamic_decompose.m 是主程序

dynamic_decp_fun.m 是分解MU程序：ICKC [1]

Dynamic_for_real_stric_MUAP1_visual.m and Dynamic_for_real_stric_MUAP2_visual.m 是针对动态分解的模块：TF method [2]

Sampling rate: 2048 HZ

2025/4/2 算法更新模块 amplitude_transform
1 Integration of Motor Unit Filters for Enhanced Surface Electromyogram Decomposition During Varying Force Isometric Contraction 2024
创新点：用于将低MVC情况下的MU的权向量W应用与高MVC，得到更好的分解结果,对于方法本身没有改进。

2 对信号进行线性放大缩小进行处理 

Reference:

[1] Y.Zheng, et al., High-Density Surface EMG Decomposition by Combining Iterative Convolution Kernel Compensation With an Energy-Specific Peel-off Strategy, TNSRE 2023

[2] Y.Zheng, et al., Low-Energy Motor Unit Extraction During Dynamic Muscle Contractions Using Time-Domain Feature Fusion Strategy, TBME 2025 (major revision)
