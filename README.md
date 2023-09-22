# small_sigmas_increase
### Repository Description:
Τhis repository contains the code for reproducing the results presented in the paper titled  [Small singular values can increase in lower precision](https://arxiv.org/abs/2303.03547) .The purpose of this code is to generate matrices $A \in \mathbb{R}^{m \times n}\$ with fixed singular values $\Sigma$ and subsequently compute the singular values in double, single, and half precision.

### Installation Requirements:
- Ensure you have Julia installed for performing major computations.
- Python is required for generating plots.
- If you're using Visual Studio Code, setting up a Python virtual environment is recommended. For instructions, please visit [this link](https://code.visualstudio.com/docs/python/python-tutorial).

Before running the code, please review the description and examples provided at the beginning of the **run.jl** script. To reproduce the examples from the paper, you can use the **run_fixed_cases.jl** script. The corresponding data and figures from the paper are located in the data/, figures/ folders respectively.
