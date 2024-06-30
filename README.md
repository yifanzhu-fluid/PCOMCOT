# PCOMCOT
PCOCMOT is a high-efficiency parallel computer program for simulating nonlinear dispersive tsunami waves. 
The motivation of this project is to enhance computational efficiency and accuracy for large-scale tsunami simulations. 
The main features of PCOMCOT include:
1) Accounting for wave dispersion by adding non-hydrostatic pressure to the shallow water model COMCOT;
2) Moving boundary technique for run-up and run-down;
3) Eddy-viscosity scheme for wave breaking;
4) One-way and two-way grid nesting for cross-scale tsunami modeling;
5) Parallel implementation with MPI library.

PCOMCOT is developed by Yifan Zhu and Prof. Chao An at Shanghai Jiao Tong University. 
We make the source code open without limitations on its redistribution and modification for research purposes. 
The current version of PCOMCOT (v2.0) is parallelized with MPI library, and a GPU-accelerated version is under development for faster performance.
We will keep updating the source code, manual, and examples when new features or modules are added.

To publish papers containing use of PCOMCOT or mentioning this model, please cite the articles about the algorithms and applications of PCOMCOT.
The articles on PCOMCOT include:
[1] Zhu, Y., C. An, H. Yu, W. Zhang, and X. Chen (2024). High-resolution tsunami hazard assessment for the Guangdong-Hong Kong-Macao Greater Bay Area based on a non-hydrostatic tsunami model. Science China Earth Sciences, 67. https://doi.org/10.1007/s11430-023-1300-9


If you find code bugs or plan to use PCOMCOT commercially, please contact us via the E-mails below.

zyftop@sjtu.edu.cn (Yifan Zhu)    
anchao@sjtu.edu.cn (Chao An)
