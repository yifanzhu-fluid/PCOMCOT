# PCOMCOT
PCOCMOT is a high-efficiency parallel computer program for simulating nonlinear dispersive tsunami waves.
It is developed by Yifan Zhu and Chao An at Shanghai Jiao Tong University, Shanghai, China.
The current version (2.1) is parallelized for CPU using the MPI library, and for GPU based on CUDA programming model.
The main features of PCOMCOT include:
1) Accounting for wave dispersion by adding non-hydrostatic pressure to the shallow water model COMCOT;
2) Moving boundary technique for run-up and run-down;
3) Eddy-viscosity scheme for wave breaking;
4) One-way and two-way grid nesting for cross-scale tsunami modeling;
5) Parallel implementation on both CPU and GPU.
   
We will keep updating the source code, manual, and examples when new features or modules are added.
To publish papers containing use of PCOMCOT or mentioning this model, please cite articles about the algorithms and applications of PCOMCOT.
The articles on PCOMCOT are: 
* Zhu, Y., An, C., Yu, H., Zhang, W., & Chen, X. (2024). High-resolution tsunami hazard assessment for the Guangdong-Hong Kong-Macao Greater Bay Area based on a non-hydrostatic tsunami model. Science China Earth Sciences, 67(7), 2326–2351. https://doi.org/10.1007/s11430-023-1300-9
* An, C., Liu, H., Ren, Z., & Yuan, Y. (2018). Prediction of tsunami waves by uniform slip models. Journal of Geophysical Research: Oceans, 123, 8366–8382. https://doi.org/10.1029/2018JC014363
* Wang, X., & Liu, P. L.-F. (2006). An analysis of 2004 Sumatra earthquake fault plane mechanisms and Indian ocean tsunami. Journal of
Hydraulic Research, 44(2), 147–154. https://doi.org/10.1080/00221686.2006.9521671

If you find code bugs or plan to use PCOMCOT commercially, please contact us via the E-mails below.

zyftop@sjtu.edu.cn (Yifan Zhu)    
anchao@sjtu.edu.cn (Chao An)
