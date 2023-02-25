# MASES
Implementation of MASES (The Maximum Separation Subspace in Sufficient Dimension Reduction with Categorical Response) with R languages. Improved with energy distance.


## Description
"method.R" file implemented some function for simulating the process of finding Maximum Separation Space, including basic SDR (Sufficient dimension reduction) methods (SIR, SAVE, DR), and also apply energy distance to find the MASES ("energy"). For detailed information, please read the report_for_MASES.pdf.

## Usage
The main function is in all.R file. Choose data generation method by "data_method"" variable and estimate methods by "estimate_method"" methods.

## Reference 
Zhang, X., Mai, Q., & Zou, H. (2020). The Maximum Separation Subspace in Sufficient Dimension Reduction with Categorical Response. Journal of Machine Learning Research, 21(29), 1-36. 