#!/usr/bin/env python3
"""
DESCUR Python Implementation
============================

A Python rewrite of DESCUR for fitting 3D space curves using steepest descent
optimization with Fourier decomposition.

This version:
- Reads R, Z data from file only
- Outputs VMEC-compatible Fourier coefficients (rbc, zbs, rbs, zbc)

Author: Python port of original Fortran DESCUR code
"""

import numpy as np
from typing import Tuple, Optional
from dataclasses import dataclass
import time
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial


@dataclass
class DescurConfig:
    """Configuration parameters for DESCUR fitting."""
    
    # Maximum dimensions
    mu: int = 15           # Maximum poloidal modes
    nu: int = 400          # Maximum poloidal points
    nv: int = 72           # Maximum toroidal planes
    
    # Convergence parameters
    ftol: float = 1e-5     # Force tolerance
    niter: int = 500      # Maximum iterations
    nstep: int = 100       # Print interval
    
    # Spectral parameters
    mexp: int = 4          # Polar damping exponent
    pexp: float = 4.0      # Spectral width p-exponent
    qexp: float = 1.0      # Spectral width q-exponent
    
    # Other parameters
    HB_parameter: float = 1.0  # VMEC-compatible compression


class DescurFitter:
    """Main class for DESCUR curve fitting."""
    
    def __init__(self, config: Optional[DescurConfig] = None):
        """Initialize DESCUR fitter.
        
        Args:
            config: Configuration object. Uses defaults if None.
        """
        self.config = config or DescurConfig()
        
        # Will be set after reading input
        self.ntheta = 0
        self.nphi = 0
        self.nfp = 1
        self.mpol = 0
        self.mrho = 0
        
        # Arrays to be allocated
        self.twopi = 2.0 * np.pi
        self.r0n = None
        self.z0n = None
        self.raxis = None
        self.zaxis = None
        
        # Logger
        self.logger = None
        
    def setup_logger(self, log_file: str = 'descur.log', console: bool = True):
        """设置日志系统。
        
        Args:
            log_file: 日志文件名
            console: 是否同时输出到控制台
        """
        self.logger = logging.getLogger('DESCUR')
        self.logger.setLevel(logging.INFO)
        
        # 清除已有的 handlers
        self.logger.handlers.clear()
        
        # 文件 handler
        file_handler = logging.FileHandler(log_file, mode='w', encoding='utf-8')
        file_handler.setLevel(logging.INFO)
        file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(file_formatter)
        self.logger.addHandler(file_handler)
        
        # 控制台 handler
        if console:
            console_handler = logging.StreamHandler()
            console_handler.setLevel(logging.INFO)
            console_formatter = logging.Formatter('%(message)s')
            console_handler.setFormatter(console_formatter)
            self.logger.addHandler(console_handler)
        
        self.logger.info("DESCUR Python Implementation")
        self.logger.info("="*70)
    
    def read_input_file(self, filename: str) -> Tuple[np.ndarray, np.ndarray]:
        """Read R, Z data from input file.
        
        File format:
            Line 1: ntheta nphi nfp
            Following lines: R phi Z  (or R Z if only one phi plane)
            
        Args:
            filename: Path to input file
            
        Returns:
            rin, zin: Arrays of shape (ntheta, nphi)
        """
        with open(filename, 'r') as f:
            # Read header
            header = f.readline().split()
            self.ntheta = int(header[0])
            self.nphi = int(header[1])
            self.nfp = int(header[2])
            
            if self.ntheta > self.config.nu:
                raise ValueError(f"ntheta ({self.ntheta}) > nu ({self.config.nu})")
            if self.nphi > self.config.nv:
                raise ValueError(f"nphi ({self.nphi}) > nv ({self.config.nv})")
            
            # Read data points
            data = []
            for line in f:
                vals = [float(x) for x in line.split()]
                if len(vals) == 3:  # R, phi, Z format
                    data.append([vals[0], vals[2]])
                elif len(vals) == 2:  # R, Z format
                    data.append(vals)
                else:
                    raise ValueError(f"Invalid data line: {line}")
            
            data = np.array(data)
            
            # Reshape to (ntheta, nphi)
            rin = data[:, 0].reshape((self.ntheta, self.nphi), order='F')
            zin = data[:, 1].reshape((self.ntheta, self.nphi), order='F')
            
            # Handle single phi plane with symmetric extension
            if self.nphi == 1 and (np.all(zin >= 0) or np.all(zin <= 0)):
                rin_full = np.zeros((2*self.ntheta, 1))
                zin_full = np.zeros((2*self.ntheta, 1))
                rin_full[:self.ntheta, 0] = rin[:, 0]
                zin_full[:self.ntheta, 0] = zin[:, 0]
                # Mirror points
                for i in range(self.ntheta-1, -1, -1):
                    if zin[i, 0] != 0:
                        idx = 2*self.ntheta - i - 1
                        rin_full[idx, 0] = rin[i, 0]
                        zin_full[idx, 0] = -zin[i, 0]
                self.ntheta = 2 * self.ntheta
                rin = rin_full
                zin = zin_full
                
        return rin, zin
    
    def prepare_from_fieldlines(self, R_lines: np.ndarray, Z_lines: np.ndarray, 
                               Phi_lines: np.ndarray, lcfs_idx: int, 
                               nfp: int = 1, nphi_descur: int = 72) -> Tuple[np.ndarray, np.ndarray]:
        """直接从磁力线追踪结果准备 DESCUR 输入数据。
        
        此方法从磁力线追踪结果中提取 LCFS 数据，并将其转换为 DESCUR 拟合所需的格式，
        避免了先写入文件再读取的中间步骤。
        
        Args:
            R_lines: R坐标数组，形状为 (nlines, nturns, nphi)
            Z_lines: Z坐标数组，形状为 (nlines, nturns, nphi)
            Phi_lines: Phi坐标数组，形状为 (nlines, nturns, nphi)
            lcfs_idx: LCFS 磁力线的索引（从0开始）
            nfp: 场周期数，默认为1
            nphi_descur: 环向角点数，默认为72
            
        Returns:
            rin, zin: 准备好的 R, Z 数组，形状为 (ntheta, nphi)
            
        Example:
            >>> from descur_python import DescurFitter
            >>> fitter = DescurFitter()
            >>> # 假设已经从磁力线追踪得到 R_lines, Z_lines, Phi_lines
            >>> rin, zin = fitter.prepare_from_fieldlines(
            ...     R_lines, Z_lines, Phi_lines, lcfs_idx=50, nfp=1, nphi_descur=72
            ... )
            >>> results = fitter.fit(rin, zin)
        """
        # 提取 LCFS 数据并转置，使其形状为 (nphi_total, nturns)
        R_lcfs = R_lines[lcfs_idx].T
        Z_lcfs = Z_lines[lcfs_idx].T
        phi_lcfs = Phi_lines[lcfs_idx].T
        

        R_slice = R_lcfs[0, :]
        Z_slice = Z_lcfs[0, :]

        nturns = R_lcfs.shape[1]
        # 初始化排序索引和访问标记
        sorted_indices = np.zeros(nturns, dtype=int)
        visited = np.zeros(nturns, dtype=bool)

        # 从第一个点开始
        current_index = 0
        sorted_indices[0] = current_index
        visited[current_index] = True

        # 贪心最近邻搜索
        for k in range(1, nturns):
            current_R = R_slice[current_index]
            current_Z = Z_slice[current_index]
            
            min_dist_sq = np.inf
            next_index = -1
            
            # 寻找最近的未访问点
            for i in range(nturns):
                if not visited[i]:
                    dist_sq = (R_slice[i] - current_R)**2 + (Z_slice[i] - current_Z)**2
                    
                    if dist_sq < min_dist_sq:
                        min_dist_sq = dist_sq
                        next_index = i
            
            # 更新当前点并记录
            if next_index != -1:
                current_index = next_index
                sorted_indices[k] = current_index
                visited[current_index] = True
            else:
                break

        # 应用排序结果
        R_lcfs_reordered = R_lcfs[:, sorted_indices]
        Z_lcfs_reordered = Z_lcfs[:, sorted_indices]
        phi_lcfs_reordered = phi_lcfs[:, sorted_indices]

        # 确保第一个环向角点为0
        phi_lcfs_reordered[0, :] = 0
        
        # 获取总环向角点数
        nphi_total = R_lcfs.shape[0]
        
        # 选择要输出的环向角点索引（均匀采样）
        descur_idx = np.linspace(0, nphi_total, nphi_descur, endpoint=False).astype(int)
        
        # 提取对应的数据
        R_lcfs_descur = R_lcfs_reordered[descur_idx]
        Z_lcfs_descur = Z_lcfs_reordered[descur_idx]
        phi_lcfs_descur = phi_lcfs_reordered[descur_idx]
        
        # 获取极向转动次数（nturns）
        nturns = R_lcfs.shape[1]
        
        # 设置参数
        self.ntheta = nturns
        self.nphi = nphi_descur
        self.nfp = nfp
        
        # 检查维度
        if self.ntheta > self.config.nu:
            raise ValueError(f"ntheta ({self.ntheta}) > nu ({self.config.nu})")
        if self.nphi > self.config.nv:
            raise ValueError(f"nphi ({self.nphi}) > nv ({self.config.nv})")
        
        # 返回格式：rin, zin 的形状为 (ntheta, nphi)
        # R_lcfs_descur 当前形状为 (nphi, nturns)，需要转置为 (nturns, nphi)
        rin = R_lcfs_descur.T
        zin = Z_lcfs_descur.T
        
        if self.logger:
            self.logger.info(f"✓ 从磁力线数据准备 DESCUR 输入:")
            self.logger.info(f"  LCFS 索引: {lcfs_idx}")
            self.logger.info(f"  场周期数 (nfp): {nfp}")
            self.logger.info(f"  极向点数 (ntheta): {self.ntheta}")
            self.logger.info(f"  环向平面数 (nphi): {self.nphi}")
            self.logger.info(f"  数据形状: rin{rin.shape}, zin{zin.shape}")
        
        return rin, zin
    
    def fit(self, rin: np.ndarray, zin: np.ndarray, log_file: str = 'descur.log') -> dict:
        """Perform DESCUR curve fitting.
        
        Args:
            rin, zin: Input R, Z arrays of shape (ntheta, nphi)
            log_file: 日志文件名
            
        Returns:
            Dictionary containing Fourier coefficients and metadata
        """
        # 设置日志
        if self.logger is None:
            self.setup_logger(log_file)
        
        # Set dimensions
        self.mpol = self.config.mu
        self.mrho = self.mpol - 1
        
        if self.nphi % 2 != 0 and self.nphi != 1:
            raise ValueError("nphi must be even or 1 for non-symmetric systems")
        
        self.nphi2 = max(1, self.nphi // 2)
        
        # Initialize arrays
        self._initialize_arrays()
        
        # Compute initial guess
        self.logger.info("Computing initial guess...")
        angle = self._compute_initial_guess(rin, zin)
        
        # Main fitting loop
        self.logger.info("\nPerforming curve fitting...")
        result1 = self._scrunch(rin, zin, angle)
        
        # Perform toroidal Fourier transform
        self.logger.info("\nPerforming toroidal Fourier transform...")
        rbc, zbs, rbs, zbc, rmnaxis, zmnaxis = self._fftrans(result1)
        
        # Return results
        return {
            'rbc': rbc,
            'zbs': zbs,
            'rbs': rbs,
            'zbc': zbc,
            'raxis': rmnaxis,
            'zaxis': zmnaxis,
            'mpol': self.mpol,
            'ntor': max(1, self.nphi - 1),
            'nfp': self.nfp,
            'nphi2': self.nphi2
        }
    
    def _initialize_arrays(self):
        """Initialize fixed arrays and parameters."""
        # Mode number arrays
        self.ntor = max(1, self.nphi - 1)
        if self.nphi == 2:
            self.ntor = 2
            
        nn0 = 1 - (self.ntor + 1) // 2
        self.nn = np.array([(nn0 + n) * self.nfp for n in range(self.ntor)])
        self.mm = np.arange(self.mpol)
        
        # Create m, n mode lists
        m1_list, n1_list = [], []
        for m in range(self.mpol):
            for n in range(self.ntor):
                if m != 0 or self.nn[n] >= 0:
                    m1_list.append(m)
                    n1_list.append(self.nn[n])
        
        self.m1 = np.array(m1_list)
        self.n1 = np.array(n1_list)
        self.mpnt = len(m1_list)
        
        # Normalization
        self.dnorm = 2.0 / self.ntheta
        
        # Spectral width weights
        self.xmpq = np.zeros((self.mpol + 1, 2))
        for m in range(1, self.mpol + 1):
            self.xmpq[m, 0] = m ** (self.config.pexp + self.config.qexp)
            self.xmpq[m, 1] = m ** self.config.pexp
        
        # Mode derivative weights
        self.dm1 = np.arange(self.mpol + 1, dtype=float)
        
        # HB weights
        self.t1m = np.zeros(self.mrho + 1)
        self.t2m = np.zeros(self.mrho + 1)
        
        self.t1m[0] = 0.0
        self.t2m[0] = self.t1m[2] if self.mrho >= 2 else 0.0
        self.t1m[1] = 1.0
        
        for m in range(2, self.mrho + 1):
            self.t1m[m] = ((m - 1) / m) ** self.config.mexp
        
        for m in range(1, self.mrho):
            self.t2m[m] = ((m + 1) / m) ** self.config.mexp
        
        self.t2m[self.mrho] = 0.0
        
        # Storage arrays
        self.r0n = np.zeros(self.nphi)
        self.z0n = np.zeros(self.nphi)
        self.raxis = np.zeros(self.nphi)
        self.zaxis = np.zeros(self.nphi)
        
    def _compute_initial_guess(self, rin: np.ndarray, zin: np.ndarray) -> np.ndarray:
        """Compute initial guess for angles and centroid.
        
        Args:
            rin, zin: Input arrays
            
        Returns:
            angle: Array of shape (ntheta, nphi)
        """
        angle = np.zeros((self.ntheta, self.nphi))
        
        self.logger.info("Ordering surface points...")
        
        for i in range(self.nphi):
            # Compute centroid
            self.r0n[i] = np.mean(rin[:, i])
            self.z0n[i] = np.mean(zin[:, i])
            self.raxis[i] = self.r0n[i]
            self.zaxis[i] = self.z0n[i]
            
            # Order points and verify axis is inside
            rin[:, i], zin[:, i], inside = self._order_points(
                rin[:, i], zin[:, i], self.raxis[i], self.zaxis[i]
            )
            
            # If axis not inside, try different points
            if not inside:
                jskip = self.ntheta // 2
                while jskip > 0:
                    jskip -= 1
                    for j1 in range(self.ntheta):
                        j2 = (j1 + jskip) % self.ntheta
                        self.raxis[i] = 0.5 * (rin[j1, i] + rin[j2, i])
                        self.zaxis[i] = 0.5 * (zin[j1, i] + zin[j2, i])
                        
                        rin[:, i], zin[:, i], inside = self._order_points(
                            rin[:, i], zin[:, i], self.raxis[i], self.zaxis[i]
                        )
                        
                        if inside:
                            self.r0n[i] = self.raxis[i]
                            self.z0n[i] = self.zaxis[i]
                            break
                    
                    if inside:
                        break
                
                if not inside:
                    raise ValueError("Could not find internal axis point")
        
        # Compute optimized angles
        angle = self._getangle(rin, zin, angle)
        
        return angle
    
    def _order_points(self, rval: np.ndarray, zval: np.ndarray,
                     xaxis: float, yaxis: float) -> Tuple[np.ndarray, np.ndarray, bool]:
        """Order points monotonically and check axis containment.
        
        Args:
            rval, zval: Point coordinates
            xaxis, yaxis: Axis position
            
        Returns:
            rval, zval: Reordered points
            inside: Whether axis is inside curve
        """
        # Order points by nearest neighbor
        ordered_idx = [0]
        remaining = set(range(1, len(rval)))
        
        while remaining:
            current = ordered_idx[-1]
            dists = np.sqrt((rval[list(remaining)] - rval[current])**2 +
                           (zval[list(remaining)] - zval[current])**2)
            nearest_local = np.argmin(dists)
            nearest = list(remaining)[nearest_local]
            ordered_idx.append(nearest)
            remaining.remove(nearest)
        
        rval = rval[ordered_idx]
        zval = zval[ordered_idx]
        
        # Check if axis is inside using winding number
        residue = 0.0
        for i in range(len(rval)):
            i_next = (i + 1) % len(rval)
            x = 0.5 * (rval[i] + rval[i_next]) - xaxis
            y = 0.5 * (zval[i] + zval[i_next]) - yaxis
            dx = rval[i_next] - rval[i]
            dy = zval[i_next] - zval[i]
            residue += (x * dy - y * dx) / (x**2 + y**2 + 1e-10)
        
        inside = True
        
        # If counterclockwise (residue > 0), reverse to make clockwise
        if residue < -0.9 * self.twopi:
            rval = np.concatenate([[rval[0]], rval[:0:-1]])
            zval = np.concatenate([[zval[0]], zval[:0:-1]])
        elif abs(residue) < 0.9 * self.twopi:
            inside = False
        
        return rval, zval, inside
    
    def _getangle(self, rval: np.ndarray, zval: np.ndarray,
                  angle: np.ndarray) -> np.ndarray:
        """Compute optimized angle offset.
        
        Args:
            rval, zval: Point coordinates
            angle: Initial angle array
            
        Returns:
            Updated angle array
        """
        # Initialize with uniform spacing
        for i in range(self.nphi):
            angle[:, i] = np.linspace(0, self.twopi, self.ntheta, endpoint=False)
        
        # Iterate to find optimal angle offset
        for iterate in range(5):
            rcos = np.zeros(self.nphi)
            rsin = np.zeros(self.nphi)
            zcos = np.zeros(self.nphi)
            zsin = np.zeros(self.nphi)
            
            for i in range(self.nphi):
                xc = rval[:, i] - self.r0n[i]
                yc = zval[:, i] - self.z0n[i]
                rcos[i] = np.sum(np.cos(angle[:, i]) * xc)
                rsin[i] = np.sum(np.sin(angle[:, i]) * xc)
                zcos[i] = np.sum(np.cos(angle[:, i]) * yc)
                zsin[i] = np.sum(np.sin(angle[:, i]) * yc)
            
            # Compute elongation
            dnum = np.sum(zsin)
            denom = np.sum(rcos)
            elongate = dnum / denom if denom != 0 else 1e10
            
            # Compute angle corrections
            delangle = 0.0
            for i in range(self.nphi):
                phiangle = np.arctan2(
                    elongate * zcos[i] - rsin[i],
                    elongate * zsin[i] + rcos[i]
                )
                delangle = max(delangle, abs(phiangle))
                angle[:, i] += phiangle
            
            if delangle < 0.02:
                break
        
        self.logger.info(f"Average elongation = {elongate:.4e}")
        self.logger.info(f"Raxis = {self.raxis[0]:.6e}, Zaxis = {self.zaxis[0]:.6e}")
        self.logger.info(f"Number of Theta Points = {self.ntheta}")
        self.logger.info(f"Number of Phi Planes = {self.nphi}")
        self.logger.info(f"Max Poloidal Mode Number = {self.mrho}")
        self.logger.info(f"Max Toroidal Mode Number = {self.ntor}")
        
        return angle
    
    def _scrunch(self, rin: np.ndarray, zin: np.ndarray,
                angle: np.ndarray, n_workers: int = None) -> np.ndarray:
        """Main optimization loop for each toroidal plane (parallel version).
        
        Args:
            rin, zin: Input arrays
            angle: Angle array
            n_workers: Number of parallel workers (None = CPU count)
            
        Returns:
            result1: Array of optimized rho coefficients (2*mpol, nphi)
        """
        result1 = np.zeros((2 * self.mpol, self.nphi))
        
        time_on = time.time()
        
        # 确定工作进程数
        if n_workers is None:
            import os
            n_workers = os.cpu_count()
        
        self.logger.info(f"\n{'='*60}")
        self.logger.info(f"并行拟合 {self.nphi} 个环向平面，使用 {n_workers} 个工作进程")
        self.logger.info(f"{'='*60}\n")
        
        # 准备所有平面的输入数据
        plane_tasks = []
        for nplane in range(self.nphi):
            plane_tasks.append((
                nplane,
                rin[:, nplane].copy(),
                zin[:, nplane].copy(),
                angle[:, nplane].copy(),
                self.r0n[nplane],
                self.z0n[nplane]
            ))
        
        # 使用进程池并行处理
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            # 提交所有任务
            futures = {
                executor.submit(
                    self._fit_single_plane,
                    nplane, rin_plane, zin_plane, angle_plane, r0n_val, z0n_val
                ): nplane
                for nplane, rin_plane, zin_plane, angle_plane, r0n_val, z0n_val in plane_tasks
            }
            
            # 收集结果
            completed = 0
            for future in as_completed(futures):
                nplane, result, log_str = future.result()
                result1[:, nplane] = result
                
                # 输出该平面的日志
                self.logger.info(log_str)
                
                completed += 1
                self.logger.info(f"\n进度: {completed}/{self.nphi} 个平面已完成\n")
        
        time_off = time.time()
        
        self.logger.info(f"\n{'='*60}")
        self.logger.info(f"ANGLE CONSTRAINTS WITH POLAR DAMPING EXPONENT = {self.config.mexp}")
        self.logger.info(f"SPECTRUM COMPUTED WITH P = {self.config.pexp:.2f}, Q = {self.config.qexp:.2f}")
        self.logger.info(f"TOTAL TIME: {time_off - time_on:.2e} SEC")
        self.logger.info(f"AVERAGE TIME PER PLANE: {(time_off - time_on)/self.nphi:.2e} SEC")
        self.logger.info(f"{'='*60}\n")
        
        return result1
    
    def _fit_single_plane(self, nplane: int, rin_plane: np.ndarray, 
                         zin_plane: np.ndarray, angle_plane: np.ndarray,
                         r0n_val: float, z0n_val: float) -> Tuple[int, np.ndarray, str]:
        """拟合单个环向平面（用于并行计算）。
        
        Args:
            nplane: 平面索引
            rin_plane: 该平面的 R 坐标
            zin_plane: 该平面的 Z 坐标
            angle_plane: 该平面的角度
            r0n_val: R 轴位置
            z0n_val: Z 轴位置
            
        Returns:
            nplane: 平面索引
            result: 优化后的 rho 系数
            log_str: 日志字符串
        """
        log_lines = []
        log_lines.append(f"\n{'='*60}")
        log_lines.append(f"Fitting toroidal plane # {nplane + 1}")
        log_lines.append(f"{'='*60}")
        log_lines.append(f"\n{'ITERATIONS':>10} {'RMS ERROR':>16} {'GRADIENT':>16} "
                        f"{'<M>':>10} {'MAX m':>8} {'DELT':>10}")
        
        n2 = 2 * (self.mrho + 1) + self.ntheta
        
        # Initialize state vector
        xvec = np.zeros(n2)
        gvec = np.zeros(n2)
        xdot = np.zeros(n2)
        xstore = np.zeros(n2)
        
        # Set initial mode amplitudes and get r10
        xvec, r10 = self._amplitud_single(r0n_val, z0n_val, angle_plane, 
                                          rin_plane, zin_plane, log_lines)
        
        # Iteration parameters
        imodes = self.mrho
        delt = 1.0
        gmin = 1e10
        too_large = 1e6
        gtrig = 1e-4
        gradient_threshold = 1e-5  # 添加梯度阈值
        
        xstore[:] = xvec[:]
        
        for iter in range(1, self.config.niter + 1):
            gvec[:] = 0.0
            
            # Compute forces and gradient - pass r10 as parameter
            fsq, gnorm, specw = self._funct(xvec, gvec, rin_plane, zin_plane, imodes, r10)
            
            if iter == 1:
                gmin = gnorm
                imodes = min(4, self.mrho)
                g11 = gnorm
            elif iter == 2:
                g11 = gnorm
            else:
                gmin = min(gmin, gnorm)
                
                # Time evolution
                g11, xvec, xdot = self._evolve(xvec, gvec, xdot, delt, gnorm, g11)
                
                # Time step control
                if gnorm / gmin > too_large:
                    irst = 2
                else:
                    irst = 1
                
                if irst == 2 or gmin == gnorm:
                    xdot[:], xvec[:], delt = self._restart(xstore, xvec, xdot, delt, irst)
                    xstore[:] = xvec[:]
                
                # Increase mode number if converged
                if gnorm < gtrig and imodes < self.mrho:
                    imodes += 1
                    xdot[:] = 0.0
                    delt = min(0.98, delt / 0.975)
            
            # Print progress
            if iter % self.config.nstep == 0 or iter <= 2 or gnorm < self.config.ftol**2:
                gout = np.sqrt(gnorm)
                modeno = imodes if iter > 1 else self.mpol - 1
                specw_out = specw ** (1.0 / self.config.qexp) if self.config.qexp > 0 else specw
                log_lines.append(f"{iter:10d} {fsq:16.3e} {gout:16.3e} "
                               f"{specw_out:10.2f} {modeno:8d} {delt:10.2e}")
            
            # 检查是否满足梯度收敛条件
            gout = np.sqrt(gnorm)
            if gout < gradient_threshold:
                log_lines.append(f"\n收敛: 梯度 {gout:.3e} < {gradient_threshold:.3e} (第 {iter} 次迭代),RMS ERROR = {fsq:.3e}")
                log_lines.append(f"{iter:10d} {fsq:16.3e} {gout:16.3e} "
                               f"{specw_out:10.2f} {modeno:8d} {delt:10.2e}")
                break
        
        result = xvec[:2*self.mpol]
        return nplane, result, '\n'.join(log_lines)
    
    def _amplitud_single(self, rcenter: float, zcenter: float, angin: np.ndarray,
                        xin: np.ndarray, yin: np.ndarray, 
                        log_lines: list = None) -> Tuple[np.ndarray, float]:
        """Initialize mode amplitudes for single plane (parallel version).
        
        Args:
            rcenter, zcenter: Center coordinates
            angin: Angle array
            xin, yin: Input R, Z
            log_lines: Optional list to append log messages
            
        Returns:
            xvec: Initialized state vector
            r10: Normalization factor
        """
        n2 = 2 * (self.mrho + 1) + self.ntheta
        xvec = np.zeros(n2)
        
        # Set axis
        xvec[0] = rcenter  # r0c
        xvec[1 + self.mrho] = zcenter  # z0c
        
        # Set angles
        xvec[2*self.mrho + 2:] = angin
        
        # Compute Fourier components
        mrz = self.mpol - 1
        xmult = 2.0 / self.ntheta
        
        r1c = np.zeros(self.config.mu)
        r1s = np.zeros(self.config.mu)
        z1c = np.zeros(self.config.mu)
        z1s = np.zeros(self.config.mu)
        
        for m in range(1, mrz + 1):
            arg = m * angin
            xi = xmult * (xin - rcenter)
            yi = xmult * (yin - zcenter)
            r1c[m-1] = np.sum(np.cos(arg) * xi)
            r1s[m-1] = np.sum(np.sin(arg) * xi)
            z1c[m-1] = np.sum(np.cos(arg) * yi)
            z1s[m-1] = np.sum(np.sin(arg) * yi)
        
        r10 = np.sqrt(r1c[0]**2 + r1s[0]**2 + z1c[0]**2 + z1s[0]**2)
        
        if log_lines is not None:
            log_lines.append(f"\nRAXIS = {rcenter:.3e}, ZAXIS = {zcenter:.3e}, R10 = {r10:.3e}")
        
        # Initialize rho coefficients (same logic as _amplitud)
        for m in range(min(mrz, self.mrho)):
            if m <= 1:
                if m + 1 < len(self.t1m):
                    t1 = self.t1m[m + 1]
                    if t1 != 0 and m < len(r1c) and m < len(z1s):
                        xvec[1 + m] = 0.5 * (r1c[m] + z1s[m]) / t1
                        xvec[2 + self.mrho + m] = 0.5 * (r1s[m] - z1c[m]) / t1
            else:
                if m + 1 < len(self.t1m) and m - 1 >= 0:
                    t1 = self.t1m[m + 1]
                    t2 = self.t2m[m - 1]
                    tnorm = t1**2 + t2**2
                    if tnorm != 0 and m < len(r1c) and m >= 2:
                        tnorm = 0.5 / tnorm
                        xvec[1 + m] = tnorm * ((r1c[m] + z1s[m]) * t1 +
                                               (r1c[m-2] - z1s[m-2]) * t2)
                        xvec[2 + self.mrho + m] = tnorm * ((r1s[m] - z1c[m]) * t1 +
                                                           (r1s[m-2] + z1c[m-2]) * t2)
        
        return xvec, r10
    
    def _funct(self, xvec: np.ndarray, gvec: np.ndarray,
              xin: np.ndarray, yin: np.ndarray, mrho_in: int, 
              r10: float = None) -> Tuple[float, float, float]:
        """Compute forces and residuals.
        
        Args:
            xvec: State vector
            gvec: Gradient vector (output)
            xin, yin: Input data
            mrho_in: Number of rho modes
            r10: Normalization factor (optional, uses self.r10 if None)
            
        Returns:
            fsq: RMS error
            gnorm: Gradient norm
            specw: Spectral width
        """
        # Use provided r10 or fall back to instance attribute
        if r10 is None:
            r10 = self.r10
        
        # Extract components
        r0c = xvec[0]
        z0c = xvec[1 + self.mrho]
        rhoc = xvec[1:1 + self.mrho]
        rhos = xvec[2 + self.mrho:2 + 2*self.mrho]
        xpts = xvec[2 + 2*self.mrho:]
        
        # Compute cos and sin arrays
        cosa = np.zeros((self.ntheta, mrho_in + 1))
        sina = np.zeros((self.ntheta, mrho_in + 1))
        
        cosa[:, 0] = 1.0
        sina[:, 0] = 0.0
        cosa[:, 1] = np.cos(xpts)
        sina[:, 1] = np.sin(xpts)
        
        for m in range(2, mrho_in + 1):
            cosa[:, m] = cosa[:, m-1] * cosa[:, 1] - sina[:, m-1] * sina[:, 1]
            sina[:, m] = sina[:, m-1] * cosa[:, 1] + cosa[:, m-1] * sina[:, 1]
        
        # Compute R, Z from rho representation
        r1 = -xin.copy()
        z1 = -yin.copy()
        rt1 = np.zeros(self.ntheta)
        zt1 = np.zeros(self.ntheta)
        
        denom = 0.0
        specw = 0.0
        
        for m in range(mrho_in + 1):
            rmc_p, rms_p, zmc_p, zms_p = self._getrz(
                r0c, z0c, rhoc, rhos, m, mrho_in
            )
            
            # Spectral width
            t2 = rmc_p**2 + zmc_p**2 + rms_p**2 + zms_p**2
            denom += t2 * self.xmpq[m, 1]
            specw += self.xmpq[m, 0] * t2
            
            # Add to R, Z
            r1 += rmc_p * cosa[:, m] + rms_p * sina[:, m]
            z1 += zmc_p * cosa[:, m] + zms_p * sina[:, m]
            
            # Derivatives
            rt1 += self.dm1[m] * (rms_p * cosa[:, m] - rmc_p * sina[:, m])
            zt1 += self.dm1[m] * (zms_p * cosa[:, m] - zmc_p * sina[:, m])
        
        specw = specw / denom if denom > 0 else 0.0
        
        # Angle forces
        gtt = rt1**2 + zt1**2
        gpts = r1 * rt1 + z1 * zt1
        gpts = 0.5 * gpts / gtt
        
        # Limit angle motion
        t1 = np.max(np.abs(gpts))
        t2 = 1e-3
        if t1 > t2:
            gpts *= t2 / t1
        
        # RMS error - use r10 parameter
        fsq = 0.5 * self.dnorm * np.sum(r1**2 + z1**2)
        fsq = np.sqrt(fsq) / r10
        
        # Compute gradient components
        
        # Axis forces
        gvec[0] = self.dnorm * np.sum(cosa[:, 0] * r1)  # gr0c
        gvec[1 + self.mrho] = self.dnorm * np.sum(cosa[:, 0] * z1)  # gz0c
        
        # Rho forces
        for m in range(mrho_in):
            if m <= 0:
                t1 = self.dnorm / max(self.t1m[m + 1], 0.1)
                gvec[1 + m] = t1 * np.sum(cosa[:, m+1] * r1 + sina[:, m+1] * z1)
                gvec[2 + self.mrho + m] = t1 * np.sum(sina[:, m+1] * r1 - cosa[:, m+1] * z1)
            else:
                t1 = self.t1m[m + 1]
                t2 = self.t2m[m - 1]
                if t1 == 0 and t2 == 0:
                    continue
                tnorm = self.dnorm / (t1**2 + t2**2)
                t1 *= tnorm
                t2 *= tnorm
                gvec[1 + m] = np.sum(
                    (cosa[:, m+1] * r1 + sina[:, m+1] * z1) * t1 +
                    (cosa[:, m-1] * r1 - sina[:, m-1] * z1) * t2
                )
                gvec[2 + self.mrho + m] = np.sum(
                    (sina[:, m+1] * r1 - cosa[:, m+1] * z1) * t1 +
                    (sina[:, m-1] * r1 + cosa[:, m-1] * z1) * t2
                )
        
        gvec[2 + self.mrho] = 0.0  # Enforce toroidal angle constraint
        
        # Angle forces
        gvec[2 + 2*self.mrho:] = gpts
        
        # Compute gradient norm - use r10 parameter
        gnorm = (np.sum(gvec[1:1+mrho_in]**2) +
                np.sum(gvec[2+self.mrho:2+self.mrho+mrho_in]**2) +
                gvec[0]**2 + gvec[1+self.mrho]**2)
        gnorm /= r10**2
        gnorm += self.dnorm * np.sum(gpts**2)
        
        return fsq, gnorm, specw
    
    def _getrz(self, r0c: float, z0c: float, rhoc: np.ndarray,
              rhos: np.ndarray, m: int, mrho_in: int) -> Tuple[float, float, float, float]:
        """Get R, Z Fourier components from rho representation.
        
        Args:
            r0c, z0c: Axis positions
            rhoc, rhos: Rho coefficients
            m: Mode number
            mrho_in: Max rho mode
            
        Returns:
            rmc, rms, zmc, zms: Fourier components
        """
        if m == 0:
            rmc = r0c + self.t2m[0] * rhoc[0] if len(rhoc) > 0 else r0c
            zms = 0.0
            rms = 0.0
            zmc = z0c + self.t2m[0] * rhos[0] if len(rhos) > 0 else z0c
        elif m < mrho_in:
            # Check bounds for rhoc[m+1]
            rhoc_m1 = rhoc[m-1] if m > 0 else 0.0
            rhoc_p1 = rhoc[m+1] if m+1 < len(rhoc) else 0.0
            rhos_m1 = rhos[m-1] if m > 0 else 0.0
            rhos_p1 = rhos[m+1] if m+1 < len(rhos) else 0.0
            
            rmc = self.t1m[m] * rhoc_m1 + self.t2m[m] * rhoc_p1
            zms = self.t1m[m] * rhoc_m1 - self.t2m[m] * rhoc_p1
            rms = self.t1m[m] * rhos_m1 + self.t2m[m] * rhos_p1
            zmc = -(self.t1m[m] * rhos_m1 - self.t2m[m] * rhos_p1)
        else:
            rhoc_m1 = rhoc[m-1] if m > 0 and m-1 < len(rhoc) else 0.0
            rhos_m1 = rhos[m-1] if m > 0 and m-1 < len(rhos) else 0.0
            rmc = self.t1m[m] * rhoc_m1 * self.config.HB_parameter
            zms = rmc * self.config.HB_parameter
            rms = self.t1m[m] * rhos_m1 * self.config.HB_parameter
            zmc = -rms * self.config.HB_parameter
        
        return rmc, rms, zmc, zms
    
    def _evolve(self, xvec: np.ndarray, gvec: np.ndarray, xdot: np.ndarray,
               delt: float, gnorm: float, g11: float) -> Tuple[float, np.ndarray, np.ndarray]:
        """Time evolution step.
        
        Args:
            xvec: State vector
            gvec: Gradient vector
            xdot: Time derivative
            delt: Time step
            gnorm: Current gradient norm
            g11: Previous gradient norm
            
        Returns:
            g11: Updated gradient norm
            xvec: Updated state
            xdot: Updated derivative
        """
        bmax = 0.15
        
        ftest = gnorm / g11 if g11 > 0 else 1.0
        dtau = abs(1.0 - ftest)
        g11 = gnorm
        
        otav = dtau / delt
        dtau = delt * otav + 1e-3
        dtau = min(bmax, dtau)
        
        b1 = 1.0 - 0.5 * dtau
        fac = 1.0 / (1.0 + 0.5 * dtau)
        
        xdot = fac * (xdot * b1 - delt * gvec)
        xvec = xvec + xdot * delt
        
        return g11, xvec, xdot
    
    def _restart(self, xstore: np.ndarray, xvec: np.ndarray,
                xdot: np.ndarray, delt: float, irst: int) -> Tuple[np.ndarray, np.ndarray, float]:
        """Reset state or store checkpoint.
        
        Args:
            xstore: Stored state
            xvec: Current state
            xdot: Time derivative
            delt: Time step
            irst: Restart flag
            
        Returns:
            xdot, xvec, delt: Updated values
        """
        if irst == 1:
            xstore[:] = xvec[:]
        elif irst == 2:
            xdot[:] = 0.0
            xvec[:] = xstore[:]
            delt *= 0.975
        
        return xdot, xvec, delt
    
    def _fftrans(self, result1: np.ndarray) -> Tuple[np.ndarray, ...]:
        """Perform toroidal Fourier transform.
        
        Args:
            result1: Array of rho coefficients (2*mpol, nphi)
            
        Returns:
            rbc, zbs, rbs, zbc: Fourier coefficients
            rmnaxis, zmnaxis: Axis coefficients
        """
        rbc = np.zeros((self.mpol, 2*self.nphi2 + 1))
        zbs = np.zeros((self.mpol, 2*self.nphi2 + 1))
        rbs = np.zeros((self.mpol, 2*self.nphi2 + 1))
        zbc = np.zeros((self.mpol, 2*self.nphi2 + 1))
        rmnaxis = np.zeros(self.nphi2 + 1)
        zmnaxis = np.zeros(self.nphi2 + 1)
        
        delphi = 1.0 / self.nphi
        intgrate = np.full(self.nphi, delphi)
        argi = self.twopi * np.arange(self.nphi) / (self.nphi * self.nfp);
        
        for idx in range(len(self.m1)):
            mreal = self.m1[idx]
            nreal = self.n1[idx] // self.nfp
            dn = float(self.n1[idx])
            
            for i in range(self.nphi):
                r0c = result1[0, i]
                z0c = result1[self.mpol, i]
                rhoc = result1[1:self.mpol, i]
                rhos = result1[self.mpol+1:2*self.mpol, i]
                
                rmc_p, rms_p, zmc_p, zms_p = self._getrz(
                    r0c, z0c, rhoc, rhos, mreal, self.mrho
                )
                
                arg = dn * argi[i]
                tcosn = np.cos(arg)
                tsinn = np.sin(arg)
                
                n_idx = nreal + self.nphi2
                rbc[mreal, n_idx] += intgrate[i] * (tcosn * rmc_p + tsinn * rms_p)
                zbs[mreal, n_idx] += intgrate[i] * (tcosn * zms_p - tsinn * zmc_p)
                zbc[mreal, n_idx] += intgrate[i] * (tcosn * zmc_p + tsinn * zms_p)
                rbs[mreal, n_idx] += intgrate[i] * (tcosn * rms_p - tsinn * rmc_p)
            
            # Axis coefficients
            if mreal == 0 and nreal == 0:
                rmnaxis[0] = np.sum(intgrate * self.raxis[:self.nphi])
                zmnaxis[0] = 0.0
            elif mreal == 0 and nreal > 0:
                n_idx = nreal + self.nphi2
                rbc[0, n_idx] *= 2
                rbs[0, n_idx] *= 2
                zbc[0, n_idx] *= 2
                zbs[0, n_idx] *= 2
                
                rmnaxis[nreal] = np.sum(2 * intgrate * self.raxis[:self.nphi] *
                                       np.cos(dn * argi))
                zmnaxis[nreal] = -np.sum(2 * intgrate * self.zaxis[:self.nphi] *
                                        np.sin(dn * argi))
        
        return rbc, zbs, rbs, zbc, rmnaxis, zmnaxis
    
    def write_output(self, results: dict, filename: str = 'outcurve'):
        """Write VMEC-compatible output file.
        
        Args:
            results: Results dictionary from fit()
            filename: Output filename
        """
        rbc = results['rbc']
        zbs = results['zbs']
        rbs = results['rbs']
        zbc = results['zbc']
        rmnaxis = results['raxis']
        zmnaxis = results['zaxis']
        
        tol = 1e-6 * abs(rbc[1, self.nphi2])
        
        with open(filename, 'w') as f:
            f.write("DESCUR Python Output - VMEC Compatible Fourier Coefficients\n")
            f.write("="*70 + "\n\n")
            
            f.write(f"MPOL = {self.mpol}\n")
            f.write(f"NTOR = {results['ntor']}\n")
            f.write(f"NFP  = {self.nfp}\n\n")
            
            f.write(f"{'MB':>5} {'NB':>4} {'RBC':>12} {'RBS':>12} "
                   f"{'ZBC':>12} {'ZBS':>12} {'RAXIS':>12} {'ZAXIS':>12}\n")
            f.write("-"*80 + "\n")
            
            for m in range(self.mpol):
                for n_idx in range(2*self.nphi2 + 1):
                    n = n_idx - self.nphi2
                    
                    if (abs(rbc[m, n_idx]) < tol and abs(zbs[m, n_idx]) < tol and
                        abs(rbs[m, n_idx]) < tol and abs(zbc[m, n_idx]) < tol):
                        continue
                    
                    if m == 0 and n >= 0:
                        f.write(f"{m:5d} {n:4d} {rbc[m,n_idx]:12.4e} {rbs[m,n_idx]:12.4e} "
                               f"{zbc[m,n_idx]:12.4e} {zbs[m,n_idx]:12.4e} "
                               f"{rmnaxis[n]:12.4e} {zmnaxis[n]:12.4e}\n")
                    else:
                        f.write(f"{m:5d} {n:4d} {rbc[m,n_idx]:12.4e} {rbs[m,n_idx]:12.4e} "
                               f"{zbc[m,n_idx]:12.4e} {zbs[m,n_idx]:12.4e}\n")
            
            # VMEC format
            f.write("\n" + "="*70 + "\n")
            f.write("VMEC NAMELIST FORMAT:\n")
            f.write("="*70 + "\n\n")
            f.write("&indata\n")
            f.write("lasym = .T.\n")
            for m in range(self.mpol):
                for n_idx in range(2*self.nphi2 + 1):
                    n = n_idx - self.nphi2
                    
                    if (abs(rbc[m, n_idx]) < tol and abs(zbs[m, n_idx]) < tol and
                        abs(rbs[m, n_idx]) < tol and abs(zbc[m, n_idx]) < tol):
                        continue
                    
                    f.write(f"  RBC({n:3d},{m:2d}) = {rbc[m,n_idx]:14.6e}   "
                           f"RBS({n:3d},{m:2d}) = {rbs[m,n_idx]:14.6e}   "
                           f"ZBC({n:3d},{m:2d}) = {zbc[m,n_idx]:14.6e}   "
                           f"ZBS({n:3d},{m:2d}) = {zbs[m,n_idx]:14.6e}\n")
            f.write("/\n")
        self.logger.info(f"\nOutput written to {filename}")


def main():
    """Main function for command-line usage."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='DESCUR Python - Fit 3D curves with Fourier decomposition'
    )
    parser.add_argument('input_file', help='Input file with R, Z data')
    parser.add_argument('-o', '--output', default='descur_output.txt',
                       help='Output filename (default: descur_output.txt)')
    parser.add_argument('--log', default='descur.log',
                       help='Log filename (default: descur.log)')
    parser.add_argument('--ftol', type=float, default=1e-8,
                       help='Force tolerance (default: 1e-8)')
    parser.add_argument('--niter', type=int, default=500,
                       help='Maximum iterations (default: 1500)')
    parser.add_argument('--mexp', type=int, default=4,
                       help='Polar damping exponent (default: 4)')
    
    args = parser.parse_args()
    
    # Create config
    config = DescurConfig(
        ftol=args.ftol,
        niter=args.niter,
        mexp=args.mexp
    )
    
    # Create fitter
    fitter = DescurFitter(config)
    fitter.setup_logger(args.log)
    
    # Read input
    fitter.logger.info(f"Reading input from {args.input_file}...")
    rin, zin = fitter.read_input_file(args.input_file)
    
    # Perform fit
    fitter.logger.info("\nStarting DESCUR fitting...\n")
    results = fitter.fit(rin, zin, log_file=args.log)
    
    # Write output
    fitter.write_output(results, args.output)
    
    fitter.logger.info("\nDESCUR fitting complete!")


if __name__ == '__main__':
    main()
