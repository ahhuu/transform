import os
import sys
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict
import tkinter as tk
from tkinter import filedialog, ttk
import threading
import re
from collections import defaultdict

# 设置中文字体（全局）
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题


class GNSSAnalyzer:
    """GNSS观测数据分析与问题检测器"""

    def __init__(self):
        # 定义GNSS信号频率 (Hz)
        self.frequencies = {
            'G': {'L1C': 1575.42e6, 'L5Q': 1176.45e6},  # GPS
            'R': {'L1C': 1602e6, 'L5Q': 1246e6},  # GLONASS
            'E': {'L1B': 1575.42e6, 'L1C': 1575.42e6, 'L5Q': 1176.45e6, 'L7Q': 1207.14e6},  # Galileo
            'C': {'L2I': 1561.098e6, 'L1P': 1575.42e6, 'L1D': 1575.42e6,  'L5P': 1176.45e6},  # BeiDou
            'J': {'L1C': 1575.42e6, 'L5Q': 1176.45e6},  # QZSS
            'I': {'L5Q': 1176.45e6, 'S': 2492.028e6},  # IRNSS/NavIC
            'S': {'L1C': 1575.42e6, 'L5Q': 1176.45e6},  # SBAS
        }

        # GLONASS PRN到k值的映射表
        self.glonass_k_map = {
            'R01': +1, 'R02': -4, 'R03': +5, 'R04': +6,
            'R05': +1, 'R06': -4, 'R07': +5, 'R08': +6,
            'R09': -2, 'R10': -7, 'R11': 0, 'R12': -1,
            'R13': -2, 'R14': -7, 'R15': 0, 'R16': -1,
            'R17': +4, 'R18': -3, 'R19': +3, 'R20': +2,
            'R21': +4, 'R22': -3, 'R23': +3, 'R24': +2
        }

        # 计算对应波长 (m)
        self.wavelengths = {}
        self.speed_of_light = 299792458  # m/s
        for system, freqs in self.frequencies.items():
            self.wavelengths[system] = {
                freq: self.speed_of_light / f for freq, f in freqs.items()
            }

        # 存储以米为单位的观测值
        self.observations_meters = {}  # 结构: {sat_id: {freq: {'times': [], 'code': [], 'phase': []}}}

        # 存储分析结果
        self.results = {
            'code_carrier_inconsistency': {},
            'observation_inconsistency': {},
            'phase_stagnation': {},
            'observable_derivatives': {},
            'code_phase_differences': {},
            'phase_prediction_errors': {}
        }

        # 结果根目录
        self.results_root = "results"
        self.results_dir = self.results_root
        os.makedirs(self.results_root, exist_ok=True)
        # 当前输入文件路径及对应结果子目录
        self.input_file_path = None
        self.current_result_dir = None
        # 定义各类图表的子文件夹名称
        self.plot_categories = {
            'raw_observations': '原始观测值',
            'derivatives': '观测值一阶差分',
            'code_phase_diffs': '伪距相位差值之差',
            'code_phase_diff_raw': '伪距相位原始差值',
            'phase_pred_errors': '相位预测误差',
            'double_differences': '历元间双差'
        }

        # 进度管理相关属性
        self.progress_callback = None
        self.current_stage = 0  # 当前阶段
        self.total_stages = 9  # 总阶段数（读取、计算差分、计算相位停滞、计算差值、计算预测误差、计算历元双差、保存新文件、保存报告、保存图表）
        self.stage_progress = 0  # 当前阶段的进度
        self.stage_max = 100  # 当前阶段的最大进度

        # 剔除粗差后的观测值文件
        self.output_format = "rinex"  # 输出文件格式，保持与原文件一致
        self.cleaned_observations = {}  # 存储清洗后的观测数据

    def set_progress_callback(self, callback):
        """设置进度回调函数"""
        self.progress_callback = callback

    def update_progress(self, value):
        """更新当前阶段的进度"""
        self.stage_progress = value
        overall_progress = (self.current_stage / self.total_stages) + \
                           (value / self.stage_max / self.total_stages)
        if self.progress_callback:
            self.progress_callback(overall_progress)

    def start_stage(self, stage_index, stage_name, max_units=100):
        """开始一个新的处理阶段"""
        self.current_stage = stage_index
        self.stage_progress = 0
        self.stage_max = max_units
        print(f"开始阶段: {stage_name}")
        if self.progress_callback:
            self.progress_callback((self.current_stage) / self.total_stages)

    def read_rinex_obs(self, file_path: str) -> Dict:
        """读取RINEX观测文件并解析数据"""
        self.input_file_path = file_path
        filename = os.path.basename(file_path)
        self.current_result_dir = os.path.join(self.results_root, filename.split('.')[0])
        os.makedirs(self.current_result_dir, exist_ok=True)

        data = {
            'header': {},
            'epochs': []
        }

        self.start_stage(0, "读取RINEX文件", 100)
        self.observations_meters = {}
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = [line.rstrip('\n') for line in f]

        # 解析头部
        header_end = 0
        for i, line in enumerate(lines):
            if 'END OF HEADER' in line:
                header_end = i + 1
                break
            # 提取头部信息
            if 'RINEX VERSION' in line:
                version_str = line[:9].strip()
                data['header']['version'] = float(version_str) if version_str else None
            elif 'MARKER NAME' in line:
                data['header']['marker'] = line[:60].strip()
            elif 'OBS TYPES' in line:
                system = line[0]
                obs_types = line[6:60].split()
                data['header'][f'obs_types_{system}'] = obs_types

        # 更新头部解析进度
        header_progress = min(int(header_end / len(lines) * 20), 20)
        self.update_progress(header_progress)

        # 解析观测数据
        current_epoch = None
        current_satellites = {}
        i = header_end
        total_lines = len(lines)

        # 为每个卫星维护独立的波长记录
        satellite_wavelengths = {}

        while i < total_lines:
            # 每处理1%的行更新一次进度
            if i % (total_lines // 80) == 0:
                progress = 20 + int((i - header_end) / (total_lines - header_end) * 80)
                self.update_progress(progress)
            line = lines[i]
            if not line.strip():
                i += 1
                continue
            if line.startswith('>'):
                # 新历元
                if current_epoch is not None:
                    data['epochs'].append({
                        'time': current_epoch,
                        'satellites': current_satellites.copy()
                    })
                try:
                    parts = line[1:].split()
                    if len(parts) >= 6:
                        year = int(parts[0])
                        month = int(parts[1])
                        day = int(parts[2])
                        hour = int(parts[3])
                        minute = int(parts[4])
                        second = int(float(parts[5]))
                        current_epoch = pd.Timestamp(
                            year=year, month=month, day=day,
                            hour=hour, minute=minute, second=second
                        )
                        current_satellites = {}
                    else:
                        print(f"警告: 时间行格式错误 (行 {i + 1}): {line}")
                        i += 1
                        continue
                except (ValueError, IndexError) as e:
                    print(f"时间解析错误 (行 {i + 1}): {line}")
                    raise ValueError(f"无效的时间格式: {line}") from e
                i += 1
            else:
                # 解析卫星观测值
                if current_epoch is None:
                    i += 1
                    continue
                if len(line) < 3:
                    i += 1
                    continue
                sat_system = line[0]
                sat_prn = line[1:3].strip()
                if not sat_prn:
                    i += 1
                    continue
                sat_id = f"{sat_system}{sat_prn}"

                # 创建当前卫星的局部频率和波长字典
                current_freqs = self.frequencies.get(sat_system, {}).copy()
                current_wavelengths = self.wavelengths.get(sat_system, {}).copy()

                # 修改后的GLONASS频率计算逻辑
                if sat_system == 'R' and sat_prn.isdigit():
                    prn = f"R{sat_prn.zfill(2)}"
                    # 从GLONASS k值映射表获取实时k值
                    k = self.glonass_k_map.get(prn, 0)
                    # 验证k值有效性
                    if not (-7 <= k <= 6):  # 根据ICD文档调整范围
                        print(f"警告: GLONASS {prn} 的k值{k}超出有效范围")
                        k = 0  # 使用默认频道
                    l1c_freq = 1602e6 + k * 0.5625e6
                    current_freqs['L1C'] = l1c_freq
                    current_wavelengths['L1C'] = self.speed_of_light / l1c_freq

                # 保存当前卫星的波长信息
                satellite_wavelengths[sat_id] = current_wavelengths.copy()

                # 获取该系统的观测类型
                obs_types = data['header'].get(f'obs_types_{sat_system}', [])
                if not obs_types:
                    i += 1
                    continue

                # 解析观测值
                observations = {}
                sat_data = line[3:]  # 跳过卫星标识
                field_width = 16
                expected_fields = len(obs_types)
                actual_fields = (len(sat_data) + field_width - 1) // field_width
                for j in range(expected_fields):
                    if j < actual_fields:
                        start_idx = j * field_width
                        end_idx = start_idx + field_width
                        field = sat_data[start_idx:end_idx].strip()
                        obs_type = obs_types[j]
                        try:
                            observations[obs_type] = float(field) if field else None
                        except ValueError:
                            observations[obs_type] = None
                    else:
                        obs_type = obs_types[j]
                        observations[obs_type] = None

                current_satellites[sat_id] = observations

                # 存储原始观测值（伪距、相位、多普勒）
                if sat_id not in self.observations_meters:
                    self.observations_meters[sat_id] = {}

                # 使用已保存的卫星特定波长值
                sat_wavelengths = satellite_wavelengths.get(sat_id, {})

                for freq in current_freqs:
                    code_obs_type = f'C{freq[1:]}'
                    phase_obs_type = f'L{freq[1:]}'
                    doppler_obs_type = f'D{freq[1:]}'

                    # 获取观测值
                    code_val = observations.get(code_obs_type)
                    phase_val = observations.get(phase_obs_type)
                    doppler_val = observations.get(doppler_obs_type)
                    wavelength = sat_wavelengths.get(freq)

                    # 初始化数据结构（确保wavelength字段存在）
                    if freq not in self.observations_meters[sat_id]:
                        self.observations_meters[sat_id][freq] = {
                            'times': [],
                            'code': [],
                            'phase': [],
                            'phase_cycle': [],
                            'doppler': [],
                            'wavelength': []  # 强制初始化波长列表
                        }
                    # 确保wavelength字段已创建
                    if 'wavelength' not in self.observations_meters[sat_id][freq]:
                        self.observations_meters[sat_id][freq]['wavelength'] = []

                    # 存储时间和伪距
                    self.observations_meters[sat_id][freq]['times'].append(current_epoch)
                    self.observations_meters[sat_id][freq]['code'].append(code_val)

                    if wavelength is None:
                        global_wavelength = self.wavelengths.get(sat_system, {}).get(freq)
                        if global_wavelength is not None:
                            wavelength = global_wavelength
                            print(f"警告: 卫星 {sat_id} 频率 {freq} 局部波长缺失，使用全局值 {wavelength}")
                        else:
                            print(f"错误: 卫星 {sat_id} 频率 {freq} 未找到波长定义")

                    # 存储相位（米）
                    if phase_val is not None and wavelength is not None:
                        self.observations_meters[sat_id][freq]['phase'].append(phase_val * wavelength)
                    else:
                        self.observations_meters[sat_id][freq]['phase'].append(None)

                    # 存储相位（周）和波长（关键修改：无论相位是否存在，都存储波长）
                    self.observations_meters[sat_id][freq]['phase_cycle'].append(phase_val)
                    self.observations_meters[sat_id][freq]['wavelength'].append(wavelength)  # 强制存储波长

                    # 存储多普勒（米/秒）
                    if doppler_val is not None and wavelength is not None:
                        self.observations_meters[sat_id][freq]['doppler'].append(-doppler_val * wavelength)
                    else:
                        self.observations_meters[sat_id][freq]['doppler'].append(None)

                i += 1  # 处理下一行

        # 添加最后一个历元
        if current_epoch is not None and current_satellites:
            data['epochs'].append({
                'time': current_epoch,
                'satellites': current_satellites.copy()
            })

        # # 输出调试信息
        # self.print_observation_debug(data)
        self.update_progress(100)
        return data

    def print_observation_debug(self, data: Dict) -> None:
        """输出各卫星/频率的原始观测值调试信息并保存为txt文件"""
        print("\n=== 原始观测数据调试输出 ===")
        debug_file_path = os.path.join(self.current_result_dir, "observation_debug.txt")

        with open(debug_file_path, 'w') as debug_file:
            debug_file.write("=== 原始观测数据调试输出 ===\n\n")

            for epoch_idx, epoch in enumerate(data['epochs'], 1):
                # 历元标题
                epoch_header = f"\n>> 历元 {epoch_idx} - {epoch['time']} <<\n"
                debug_file.write(epoch_header)
                print(epoch_header, end='')

                for sat_id, obs in epoch['satellites'].items():
                    system = sat_id[0]
                    # 获取当前卫星的频率（使用存储的局部频率，而非全局）
                    available_freqs = self.observations_meters.get(sat_id, {}).keys()
                    # 卫星标题
                    sat_header = f"\n * 卫星 {sat_id}\n"
                    debug_file.write(sat_header)
                    print(sat_header, end='')

                    for freq in available_freqs:
                        # 获取观测类型和值
                        code_obs_type = f'C{freq[1:]}'
                        phase_obs_type = f'L{freq[1:]}'
                        doppler_obs_type = f'D{freq[1:]}'
                        snr_obs_type = f'S{freq[1:]}'
                        code_val = obs.get(code_obs_type)
                        phase_val = obs.get(phase_obs_type)
                        doppler_val = obs.get(doppler_obs_type)
                        snr_val = obs.get(snr_obs_type)

                        # 跳过全空观测
                        if all(v is None for v in [code_val, phase_val, doppler_val, snr_val]):
                            continue

                        # 从存储数据获取计算值和波长（关键修改：使用该卫星存储的波长）
                        stored_data = self.observations_meters.get(sat_id, {}).get(freq, {})
                        phase_in_meters = stored_data.get('phase', [None])[-1]
                        doppler_in_mps = stored_data.get('doppler', [None])[-1]
                        # 优先从存储的波长列表中取当前历元的值（若存在）
                        wavelength = stored_data.get('wavelength', [None])[-1] if 'wavelength' in stored_data else None

                        # 构建频率块输出
                        output = [
                            f"  频率 {freq}:",
                            f"    - 伪距 ({code_obs_type}): {code_val:.3f} m" if code_val is not None else "- 伪距: None",
                            f"    - 相位 ({phase_obs_type}): {phase_val:.3f} 周" if phase_val is not None else "- 相位: None",
                            f"    - 波长: {wavelength:.6f} m" if wavelength is not None else "- 波长: None",
                            f"    - 相位(米): {phase_in_meters:.3f} m" if phase_in_meters is not None else "- 相位(米): None",
                            f"    - 多普勒 ({doppler_obs_type}): {doppler_val:.3f} Hz" if doppler_val is not None else "- 多普勒: None",
                            f"    - 多普勒速度: {doppler_in_mps:.3f} m/s" if doppler_in_mps is not None else "- 多普勒速度: None",
                            f"    - 信噪比 ({snr_obs_type}): {snr_val:.1f} dBHz" if snr_val is not None else "- 信噪比: None",
                            ""
                        ]

                        # 输出到控制台和文件
                        block_output = "\n".join(output)
                        print(block_output)
                        debug_file.write(block_output + "\n")

                print(f"\n调试信息已保存至: {debug_file_path}")

    def calculate_observable_derivatives(self, data: Dict) -> Dict:
        """计算每个卫星每个频率伪距、相位与多普勒观测的一阶差分"""
        self.start_stage(1, "计算观测值一阶差分", 100)

        derivatives = {}
        total_sats = len(data['epochs'][0]['satellites'].keys())
        processed_sats = 0

        # 遍历所有卫星
        for sat_id in data['epochs'][0]['satellites'].keys():
            system = sat_id[0]
            freq_derivatives = {}
            available_freqs = self.frequencies.get(system, {})

            # 初始化频率数据结构
            for freq in available_freqs:
                freq_derivatives[freq] = {
                    'times': [],
                    'pr_derivative': [],  # (m/s)
                    'ph_derivative': [],  # (m/s)
                    'doppler': []  # (m/s)
                }

            # 收集该卫星在各历元的观测值
            for epoch in data['epochs']:
                if sat_id in epoch['satellites']:
                    obs = epoch['satellites'][sat_id]
                    time = epoch['time']
                    for freq in available_freqs:
                        code_obs_type = f'C{freq[1:]}'
                        phase_obs_type = f'L{freq[1:]}'
                        doppler_obs_type = f'D{freq[1:]}'
                        code_value = obs.get(code_obs_type)
                        phase_value = obs.get(phase_obs_type)
                        doppler_value = obs.get(doppler_obs_type)
                        if code_value is not None or phase_value is not None or doppler_value is not None:
                            freq_derivatives[freq]['times'].append(time)
                            freq_derivatives[freq]['pr_derivative'].append(code_value)
                            freq_derivatives[freq]['ph_derivative'].append(phase_value)
                            freq_derivatives[freq]['doppler'].append(doppler_value)

            # 计算一阶差分
            for freq in available_freqs:
                times = freq_derivatives[freq]['times']
                pr_values = freq_derivatives[freq]['pr_derivative']
                ph_values = freq_derivatives[freq]['ph_derivative']
                doppler_values = freq_derivatives[freq]['doppler']
                wavelength = self.wavelengths[system].get(freq)
                if wavelength is None:
                    continue

                # 伪距一阶差分 (m/s)
                pr_derivatives = []
                for i in range(1, len(times)):
                    time_diff = (times[i] - times[i - 1]).total_seconds()
                    if time_diff > 0 and pr_values[i] is not None and pr_values[i - 1] is not None:
                        pr_derivatives.append((pr_values[i] - pr_values[i - 1]) / time_diff)
                    else:
                        pr_derivatives.append(None)

                # 相位一阶差分 (m/s)
                ph_derivatives = []
                for i in range(1, len(times)):
                    time_diff = (times[i] - times[i - 1]).total_seconds()
                    if time_diff > 0 and ph_values[i] is not None and ph_values[i - 1] is not None:
                        phase_rate_cycles = (ph_values[i] - ph_values[i - 1]) / time_diff
                        phase_rate_meters = phase_rate_cycles * wavelength  # 转换为 m/s
                        ph_derivatives.append(phase_rate_meters)
                    else:
                        ph_derivatives.append(None)

                # 多普勒值转换为 m/s
                doppler_meters = []
                for doppler in doppler_values:
                    if doppler is not None:
                        doppler_meters.append(-doppler * wavelength)
                    else:
                        doppler_meters.append(None)

                # 保存结果（移除第一个时间点，因为差分计算少一个数据点）
                freq_derivatives[freq]['pr_derivative'] = pr_derivatives
                freq_derivatives[freq]['ph_derivative'] = ph_derivatives
                freq_derivatives[freq]['doppler'] = doppler_meters
                freq_derivatives[freq]['times'] = times[1:]

            # 保存该卫星的所有频率差分结果
            derivatives[sat_id] = freq_derivatives

            processed_sats += 1
            self.update_progress(int(processed_sats / total_sats * 100))

        # 存储结果
        self.results['observable_derivatives'] = derivatives
        return derivatives

    def detect_phase_stagnation(self, data: Dict, threshold_cycles: float = 0.1, min_consecutive: int = 5) -> Dict:
        """
        检测载波相位停滞
        参数:
            data: RINEX观测数据
            threshold_cycles: 相位变化阈值（周）
            min_consecutive: 最小连续停滞历元数
        返回:
            停滞检测结果字典
        """
        self.start_stage(2, "检测载波相位停滞", 100)

        stagnation_results = {}
        total_sats = len(self.observations_meters)
        processed_sats = 0

        for sat_id, freq_data in self.observations_meters.items():
            freq_stagnation = {}

            for freq, obs_data in freq_data.items():
                times = obs_data['times']
                phase_cycles = obs_data['phase_cycle']  # 相位（周）

                # 初始化停滞检测结果
                stagnant_epochs = []  # 停滞历元索引
                current_streak = 0  # 当前连续停滞计数
                max_streak = 0  # 最大连续停滞数

                for i in range(1, len(phase_cycles)):
                    # 检查当前历元和前一历元的相位变化
                    if (phase_cycles[i] is not None and
                            phase_cycles[i - 1] is not None and
                            abs(phase_cycles[i] - phase_cycles[i - 1]) < threshold_cycles):
                        current_streak += 1
                    else:
                        # 相位变化超过阈值，重置计数
                        current_streak = 0

                    # 更新最大连续停滞数
                    if current_streak > max_streak:
                        max_streak = current_streak

                    # 记录停滞历元（当连续停滞数达到最小阈值时）
                    if current_streak >= min_consecutive:
                        stagnant_epochs.append(i)

                # 存储该频率的停滞结果
                freq_stagnation[freq] = {
                    'is_stagnant': max_streak >= min_consecutive,
                    'max_stagnant_epochs': max_streak,
                    'stagnant_epochs': stagnant_epochs,
                    'threshold': threshold_cycles,
                    'min_consecutive': min_consecutive
                }

            # 存储该卫星的所有频率停滞结果
            stagnation_results[sat_id] = freq_stagnation

            processed_sats += 1
            self.update_progress(int(processed_sats / total_sats * 100))

        # 存储结果
        self.results['phase_stagnation'] = stagnation_results
        return stagnation_results

    def calculate_code_phase_differences(self, data: Dict) -> Dict:
        """计算每个卫星每个频率的伪距与相位观测值的差值及差值变化率"""
        self.start_stage(3, "计算伪距相位差值", 100)

        # 先执行相位停滞检测
        if not self.results.get('phase_stagnation'):
            self.detect_phase_stagnation(data)

        differences = {}
        total_sats = len(self.observations_meters)
        processed_sats = 0

        # 遍历所有卫星
        for sat_id, freq_data in self.observations_meters.items():
            freq_differences = {}

            # 获取该卫星的相位停滞结果
            sat_stagnation = self.results['phase_stagnation'].get(sat_id, {})

            # 为每个频率计算差值
            for freq, obs_data in freq_data.items():
                times = obs_data['times']
                code_values = obs_data['code']  # 米
                phase_values = obs_data['phase']  # 米

                # 获取该频率的停滞历元索引
                stagnant_epochs = sat_stagnation.get(freq, {}).get('stagnant_epochs', [])

                # 初始化结果存储
                freq_differences[freq] = {
                    'times': [],
                    'code_phase_diff': [],  # 伪距相位原始差值
                    'diff_changes': [],  # 差值的历元间变化
                    'original_epochs': len(times),
                    'filtered_epochs': 0,
                    'stagnant_epochs_removed': len(stagnant_epochs),
                    'missing_epochs': 0  # 新增：缺失观测值的历元数
                }

                prev_diff = None  # 上一历元的差值，用于计算变化
                missing_obs = 0  # 记录缺失观测值的历元数

                # 计算差值，跳过停滞历元
                for i in range(len(times)):
                    # 跳过停滞历元或观测值缺失的历元
                    if (i in stagnant_epochs or
                            code_values[i] is None or
                            phase_values[i] is None):
                        if code_values[i] is None or phase_values[i] is None:
                            missing_obs += 1
                        continue

                    diff = code_values[i] - phase_values[i]  # 直接计算（单位：米）
                    freq_differences[freq]['times'].append(times[i])
                    freq_differences[freq]['code_phase_diff'].append(diff)

                    # 计算差值的变化（仅当前后两个历元都有有效值时）
                    if prev_diff is not None:
                        diff_change = abs(diff - prev_diff)
                        freq_differences[freq]['diff_changes'].append(diff_change)
                    else:
                        freq_differences[freq]['diff_changes'].append(None)  # 第一个历元没有变化

                    prev_diff = diff

                # 更新统计信息
                freq_differences[freq]['filtered_epochs'] = len(freq_differences[freq]['code_phase_diff'])
                freq_differences[freq]['missing_epochs'] = missing_obs

            # 保存该卫星的所有频率差值结果
            differences[sat_id] = freq_differences

            processed_sats += 1
            self.update_progress(int(processed_sats / total_sats * 100))

        # 存储结果
        self.results['code_phase_differences'] = differences
        return differences

    def calculate_phase_prediction_errors(self, data: Dict) -> Dict:
        """计算相位预测误差"""
        self.start_stage(4, "计算相位预测误差", 100)

        errors = {}
        total_sats = len(self.observations_meters)
        processed_sats = 0

        for sat_id, freq_data in self.observations_meters.items():
            freq_errors = {}

            for freq, obs_data in freq_data.items():
                times = obs_data['times']
                phase_values = obs_data['phase_cycle']  # 周
                doppler_values = obs_data['doppler']  # Hz

                # 获取对应频率(Hz)
                frequency = self.frequencies[sat_id[0]].get(freq)
                wavelength = self.wavelengths[sat_id[0]].get(freq)

                # 初始化结果存储
                freq_errors[freq] = {
                    'times': [],
                    'actual_phase': [],  # 周
                    'predicted_phase': [],  # 周
                    'prediction_error': [],  # 米
                    'doppler': []  # Hz
                }

                # 计算预测误差
                for i in range(1, len(times)):
                    # 检查数据有效性
                    if (phase_values[i - 1] is not None and
                            doppler_values[i - 1] is not None and
                            i < len(doppler_values) and
                            doppler_values[i] is not None and
                            phase_values[i] is not None and
                            frequency is not None and
                            wavelength is not None):
                        time_diff = (times[i] - times[i - 1]).total_seconds()
                        doppler_mean = (doppler_values[i - 1] + doppler_values[i]) / 2

                        # 正确计算相位变化率(周/秒)
                        phase_rate = -doppler_mean / frequency  # 周/秒
                        predicted_phase = phase_values[i - 1] + phase_rate * time_diff  # 周

                        # 计算误差(米)
                        error = (phase_values[i] - predicted_phase) * wavelength  # 米

                        # 保存结果
                        freq_errors[freq]['times'].append(times[i])
                        freq_errors[freq]['actual_phase'].append(phase_values[i])
                        freq_errors[freq]['predicted_phase'].append(predicted_phase)
                        freq_errors[freq]['prediction_error'].append(error)
                        freq_errors[freq]['doppler'].append(doppler_mean)

                errors[sat_id] = freq_errors
                processed_sats += 1
                self.update_progress(int(processed_sats / total_sats * 100))

        self.results['phase_prediction_errors'] = errors
        return errors

    def calculate_epoch_double_differences(self):
        """计算各卫星各频率的历元间双差（伪距、相位、多普勒）"""
        self.start_stage(5, "计算历元间双差", 100)
        double_diffs = {}  # 存储双差结果 {sat_id: {freq: {dd_code: [], dd_phase: [], dd_doppler: []}}}

        for sat_id, freq_data in self.observations_meters.items():
            double_diffs[sat_id] = {}
            for freq, data in freq_data.items():
                code = np.array(data['code'], dtype=float)  # 伪距（米）
                phase = np.array(data['phase'], dtype=float)  # 载波相位（米）
                doppler = np.array(data['doppler'], dtype=float)  # 多普勒（周/秒）

                # 计算双差（i>2时，dd = x[i+2] - 2x[i+1] + x[i]）
                n = len(code)
                if n < 3:
                    continue  # 至少需要3个历元才能计算双差

                dd_code = np.zeros(n - 2)
                dd_phase = np.zeros(n - 2)
                dd_doppler = np.zeros(n - 2)

                for i in range(n - 2):
                    dd_code[i] = code[i + 2] - 2 * code[i + 1] + code[i]
                    dd_phase[i] = phase[i + 2] - 2 * phase[i + 1] + phase[i]
                    dd_doppler[i] = doppler[i + 2] - 2 * doppler[i + 1] + doppler[i]

                # 存储结果（剔除前两个历元，双差结果长度为n-2）
                double_diffs[sat_id][freq] = {
                    'times': data['times'][2:],  # 双差对应的时间为第3个历元起
                    'dd_code': dd_code.tolist(),
                    'dd_phase': dd_phase.tolist(),
                    'dd_doppler': dd_doppler.tolist()
                }

        self.results['double_differences'] = double_diffs  # 保存双差结果
        self.update_progress(100)
        return double_diffs

    def calculate_triple_median_error(self, double_diffs):
        """计算双差结果的三倍中误差并检测超限值（包含伪距、相位、多普勒）"""
        triple_errors = {}  # 存储三倍中误差及超限值 {sat_id: {freq: {threshold: float, outliers: list}}}

        for sat_id, freq_data in double_diffs.items():
            triple_errors[sat_id] = {}
            for freq, dd_data in freq_data.items():
                # 提取双差观测值（过滤无效值）
                code = np.array([v for v in dd_data['dd_code'] if v is not None and not np.isnan(v)])
                phase = np.array([v for v in dd_data['dd_phase'] if v is not None and not np.isnan(v)])
                doppler = np.array([v for v in dd_data['dd_doppler'] if v is not None and not np.isnan(v)])

                # 计算中误差（使用标准差）
                def std_error(arr):
                    if len(arr) < 2:  # 至少需要两个观测值计算标准差
                        return 0
                    return np.std(arr, ddof=1)  # 样本标准差（自由度n-1）

                # 伪距双差
                if len(code) > 1:
                    sigma_code = std_error(code)
                    triple_sigma_code = 3 * sigma_code if sigma_code != 0 else 0.1  # 三倍中误差
                    outliers_code = np.where(np.abs(dd_data['dd_code']) > triple_sigma_code)[0].tolist()  # 超限索引
                else:
                    triple_sigma_code, outliers_code = 0, []

                # 载波相位双差
                if len(phase) > 1:
                    sigma_phase = std_error(phase)
                    triple_sigma_phase = 3 * sigma_phase if sigma_phase != 0 else 0.01  # 相位精度更高，阈值更小
                    outliers_phase = np.where(np.abs(dd_data['dd_phase']) > triple_sigma_phase)[0].tolist()
                else:
                    triple_sigma_phase, outliers_phase = 0, []

                # 多普勒双差
                if len(doppler) > 1:
                    sigma_doppler = std_error(doppler)
                    triple_sigma_doppler = 3 * sigma_doppler if sigma_doppler != 0 else 0.05  # 多普勒阈值
                    outliers_doppler = np.where(np.abs(dd_data['dd_doppler']) > triple_sigma_doppler)[0].tolist()
                else:
                    triple_sigma_doppler, outliers_doppler = 0, []

                # 合并结果
                triple_errors[sat_id][freq] = {
                    'code': {'threshold': triple_sigma_code, 'outliers': outliers_code},
                    'phase': {'threshold': triple_sigma_phase, 'outliers': outliers_phase},
                    'doppler': {'threshold': triple_sigma_doppler, 'outliers': outliers_doppler}
                }

        self.results['triple_median_errors'] = triple_errors
        return triple_errors

    # 剔除粗差保存新文件
    def remove_outliers_and_save(self, double_diffs, triple_errors):
        """
        基于伪距、相位、多普勒双差检测结果，修改RINEX文件中的异常观测值
        """
        # 各观测类型的最大阈值限制（根据观测精度设置）
        max_threshold_limit = {
            'code': 10.0,  # 伪距（米）
            'phase': 3.0,  # 相位（米）
            'doppler': 5.0  # 多普勒（米/秒）
        }

        # 初始化日志内容
        log_content = [
            "=" * 70 + "\n",
            "RINEX 粗差处理详细日志\n",
            "=" * 70 + "\n\n",
            f"观测值最大阈值限制: 伪距={max_threshold_limit['code']}m, 相位={max_threshold_limit['phase']}m, 多普勒={max_threshold_limit['doppler']}m/s\n\n"
        ]

        # 1. 整理每个卫星每个频率的阈值
        # 存储结构：{sat_id: {freq: {obs_type: threshold}}}
        satellite_freq_thresholds = defaultdict(lambda: defaultdict(dict))

        for sat_id, freq_data in triple_errors.items():
            for freq, errors in freq_data.items():
                for obs_type in ['code', 'phase', 'doppler']:
                    # 从triple_errors中获取已计算的三倍中误差
                    triple_sigma = errors.get(obs_type, {}).get('threshold', 0)

                    if triple_sigma <= 0:
                        # 无有效三倍中误差时使用最大阈值限制
                        threshold = max_threshold_limit[obs_type]
                        log_content.append(
                            f"卫星 {sat_id} 频率 {freq} {obs_type} 无有效三倍中误差，使用最大阈值: {threshold:.4f}m\n"
                        )
                    else:
                        # 确保阈值不超过最大限制
                        threshold = min(triple_sigma, max_threshold_limit[obs_type])
                        # 确保阈值为正数
                        threshold = max(threshold, 0.01)  # 避免阈值过小导致误判
                        log_content.append(
                            f"卫星 {sat_id} 频率 {freq} {obs_type} 使用三倍中误差作为阈值: "
                            f"计算值={triple_sigma:.4f}m, 应用阈值={threshold:.4f}m\n"
                        )
                    satellite_freq_thresholds[sat_id][freq][obs_type] = threshold

        # 2. 读取RINEX文件内容
        with open(self.input_file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()

        # 3. 解析历元时间戳
        epoch_timestamps = {}
        current_epoch = 0
        for line in lines:
            if line.startswith('>'):
                current_epoch += 1
                parts = line[1:].split()
                if len(parts) >= 6:
                    timestamp = ' '.join(parts[:6])
                    epoch_timestamps[current_epoch] = timestamp
        log_content.append(f"成功解析 {len(epoch_timestamps)} 个历元时间戳\n\n")

        # 4. 解析观测类型（同时处理伪距、相位、多普勒）
        system_obs_info = {}  # {系统: {'obs_types': [], 'freq_to_indices': {频率: {观测类型: 字段索引}}}}
        for line in lines:
            if 'SYS / # / OBS TYPES' in line:
                system = line[0]
                obs_types = line.split()[2:]
                freq_to_indices = defaultdict(dict)  # {频率: {'code': 索引, 'phase': 索引, 'doppler': 索引}}

                for idx, obs in enumerate(obs_types):
                    if obs.startswith('C'):  # 伪距
                        freq = f"L{obs[1:]}"
                        freq_to_indices[freq]['code'] = idx
                    elif obs.startswith('L'):  # 相位
                        freq = f"L{obs[1:]}"
                        freq_to_indices[freq]['phase'] = idx
                    elif obs.startswith('D'):  # 多普勒
                        freq = f"L{obs[1:]}"
                        freq_to_indices[freq]['doppler'] = idx

                system_obs_info[system] = {
                    'obs_types': obs_types,
                    'freq_to_indices': freq_to_indices
                }
                log_content.append(f"解析系统 {system} 观测类型: {len(obs_types)} 种\n")

        # 5. 识别异常历元（使用每个卫星每个频率的阈值）
        outlier_epochs = defaultdict(lambda: defaultdict(list))  # {sat_id: {历元: [(观测类型, 频率)]}}
        outlier_details = defaultdict(list)  # {sat_id: [异常详情]}

        for sat_id, freq_data in double_diffs.items():
            for freq, dd_data in freq_data.items():
                # 检查三种观测值的双差超限情况
                for obs_type in ['code', 'phase', 'doppler']:
                    dd_key = f"dd_{obs_type}"
                    if dd_key not in dd_data:
                        continue  # 跳过无数据的观测类型

                    # 获取该卫星该频率的阈值（已考虑最大限制）
                    threshold = satellite_freq_thresholds[sat_id][freq].get(
                        obs_type, max_threshold_limit[obs_type]
                    )

                    # 检测超限值
                    valid_dd = [(i, d) for i, d in enumerate(dd_data[dd_key])
                                if d is not None and not np.isnan(d)]
                    for orig_idx, dd_value in valid_dd:
                        if abs(dd_value) > threshold:
                            epoch_idx = orig_idx + 2  # 双差结果对应原历元（滞后2个）
                            timestamp = epoch_timestamps.get(epoch_idx, f"未知时间戳(历元{epoch_idx})")
                            outlier_info = {
                                'obs_type': obs_type,
                                'freq': freq,
                                'dd_value': dd_value,
                                'threshold_used': threshold,
                                'epoch_idx': epoch_idx,
                                'timestamp': timestamp
                            }
                            outlier_details[sat_id].append(outlier_info)
                            outlier_epochs[sat_id][epoch_idx].append((obs_type, freq))

        # 6. 修改异常历元的观测值（同时处理三种观测类型）
        modified_count = defaultdict(int)  # {观测类型: 修改数量}
        modified_satellites = set()
        satellite_modify_details = []

        for sat_id, epoch_obs_map in outlier_epochs.items():
            sat_system = sat_id[0]
            sat_prn = sat_id[1:].zfill(2)
            system_info = system_obs_info.get(sat_system, {})
            freq_indices = system_info.get('freq_to_indices', {})  # 频率到字段索引的映射

            if not freq_indices:
                continue  # 跳过无观测类型信息的卫星

            satellite_modifications = []

            # 处理每个异常历元
            for epoch_idx, obs_freq_list in epoch_obs_map.items():
                timestamp = epoch_timestamps.get(epoch_idx, f"未知时间戳(历元{epoch_idx})")

                # 定位历元行
                epoch_start = -1
                for i, line in enumerate(lines):
                    if line.startswith('>') and timestamp in line:
                        epoch_start = i
                        break
                if epoch_start < 0:
                    continue

                # 定位该卫星在历元中的数据行
                sat_line_idx = -1
                j = epoch_start + 1
                while j < len(lines) and not lines[j].startswith('>'):
                    line = lines[j]
                    if len(line) >= 3 and line[0] == sat_system and line[1:3].strip().zfill(2) == sat_prn:
                        sat_line_idx = j
                        break
                    j += 1
                if sat_line_idx < 0:
                    continue

                # 修改异常观测值字段
                original_line = lines[sat_line_idx]
                modified_line = list(original_line)
                field_modified = False
                modified_fields = []

                for obs_type, freq in obs_freq_list:
                    # 定位该观测值在数据行中的位置
                    if freq not in freq_indices or obs_type not in freq_indices[freq]:
                        continue  # 跳过无索引的观测值
                    field_idx = freq_indices[freq][obs_type]
                    start_pos = 3 + field_idx * 16  # 3: 卫星标识长度
                    end_pos = start_pos + 16

                    if end_pos > len(modified_line):
                        continue  # 防止越界

                    # 清除该字段（设为空格）
                    original_field = original_line[start_pos:end_pos].strip()
                    if original_field:
                        modified_line[start_pos:end_pos] = ' ' * 16  # RINEX标准字段宽度为16字符
                        modified_count[obs_type] += 1
                        field_modified = True
                        modified_fields.append(f"{freq}({obs_type})")

                if field_modified:
                    lines[sat_line_idx] = ''.join(modified_line)
                    modified_satellites.add(sat_id)
                    satellite_modifications.append(
                        f"  历元 {epoch_idx} ({timestamp}): 已修改 {', '.join(modified_fields)}")

            if satellite_modifications:
                satellite_modify_details.append(f"卫星 {sat_id} 的修改详情:")
                satellite_modify_details.extend(satellite_modifications)
                satellite_modify_details.append("")

        # 7. 保存修改后的文件
        results_dir = self.current_result_dir
        file_name = os.path.basename(self.input_file_path)
        modified_file_name = f"mod-{file_name}"
        output_path = os.path.join(results_dir, modified_file_name)

        with open(output_path, 'w', encoding='utf-8') as f:
            f.writelines(lines)

        # 8. 生成详细日志
        total_modified = sum(modified_count.values())
        log_content.append("\n一、修改统计摘要\n")
        log_content.append("-" * 70 + "\n")
        log_content.append(f"总计修改卫星数: {len(modified_satellites)}\n")
        log_content.append(f"总计修改观测值: {total_modified}\n")
        log_content.append(
            f"修改分类: 伪距={modified_count['code']}, 相位={modified_count['phase']}, 多普勒={modified_count['doppler']}\n\n")

        # 添加异常历元详情到日志
        log_content.append("\n二、异常历元检测详情\n")
        log_content.append("-" * 70 + "\n")

        # 按系统组织卫星信息
        system_satellites = defaultdict(list)
        for sat_id in outlier_details.keys():
            system_satellites[sat_id[0]].append(sat_id)

        # 按系统顺序输出
        for system, satellites in system_satellites.items():
            log_content.append(f"卫星系统 {system}:\n")

            for sat_id in sorted(satellites):
                details = outlier_details[sat_id]
                log_content.append(f"  卫星 {sat_id} ({len(details)}个异常观测值):\n")

                # 按历元分组异常信息
                epoch_groups = defaultdict(list)
                for detail in details:
                    epoch_groups[detail['epoch_idx']].append(detail)

                for epoch_idx in sorted(epoch_groups.keys()):
                    details = epoch_groups[epoch_idx]
                    timestamp = details[0]['timestamp']
                    log_content.append(f"    历元 {epoch_idx} ({timestamp}):\n")

                    for detail in details:
                        log_content.append(f"      - {detail['obs_type']}@{detail['freq']}: "
                                           f"双差值={detail['dd_value']:.6f}m, "
                                           f"阈值={detail['threshold_used']:.6f}m, "
                                           f"状态={'已修改' if sat_id in modified_satellites else '未修改'}\n")
                log_content.append("\n")

        # 写入日志文件
        debug_file_name = "cleaning_debug.log"
        debug_path = os.path.join(results_dir, debug_file_name)
        with open(debug_path, 'w', encoding='utf-8') as f:
            f.writelines(log_content)

        # 控制台输出
        print(f"--成功解析 {len(epoch_timestamps)} 个历元时间戳")
        print(f"--检测到 {len(outlier_epochs)} 颗卫星存在异常历元")
        print(f"--总计修改 {total_modified} 个观测值，涉及 {len(modified_satellites)} 颗卫星")
        print(f"   - 伪距: {modified_count['code']}")
        print(f"   - 相位: {modified_count['phase']}")
        print(f"   - 多普勒: {modified_count['doppler']}")
        print(f"--剔除粗差后的文件保存至: {output_path}")
        print(f"--剔除粗差详细信息保存至: {debug_path}")

        return output_path

    def plot_raw_observations(self, sat_id: str, save=True) -> None:
        """绘制指定卫星的原始伪距和载波相位随历元变化图"""
        if sat_id not in self.observations_meters:
            print(f"错误: 未找到卫星 {sat_id} 的观测数据")
            return

        satellite_data = self.observations_meters[sat_id]
        valid_freqs = [freq for freq, data in satellite_data.items()
                       if len(data['code']) > 0 or len(data['phase']) > 0]

        if not valid_freqs:
            print(f"错误: 卫星 {sat_id} 没有有效的观测数据")
            return

        plt.figure(figsize=(12, 6))

        # 定义样式字典，为每个频率分配不同的颜色、线型和标记
        style_dict = {
            'L1C': {'color': 'blue', 'linestyle': '-', 'marker': 's', 'label_code': 'L1C 伪距',
                    'label_phase': 'L1C 相位', 'phase_marker': 'o'},
            'L1D': {'color': 'cyan', 'linestyle': '-', 'marker': '^', 'label_code': 'L1D 伪距',
                    'label_phase': 'L1D 相位', 'phase_marker': 'x'},
            'L1P': {'color': 'red', 'linestyle': '-', 'marker': 'D', 'label_code': 'L1P 伪距',
                    'label_phase': 'L1P 相位', 'phase_marker': '*'},
            'L2I': {'color': 'blue', 'linestyle': '-', 'marker': 'o', 'label_code': 'L2I 伪距',
                    'label_phase': 'L2I 相位', 'phase_marker': '^'},
            'L5Q': {'color': 'cyan', 'linestyle': '-', 'marker': '^', 'label_code': 'L5Q 伪距',
                    'label_phase': 'L5Q 相位', 'phase_marker': 's'},
            'L7Q': {'color': 'magenta', 'linestyle': '-', 'marker': '^', 'label_code': 'L7Q 伪距',
                    'label_phase': 'L7Q 相位', 'phase_marker': 'D'},
            'L5P': {'color': 'magenta', 'linestyle': '-', 'marker': 'D', 'label_code': 'L5P 伪距',
                    'label_phase': 'L5P 相位', 'phase_marker': '^'}
        }

        # 记录是否有绘制的线条（用于判断是否添加图例）
        has_plotted_lines = False

        # 绘制所有频率的伪距和相位
        for idx, freq in enumerate(valid_freqs):
            data = satellite_data[freq]
            times = data['times']
            code_values = data['code']
            phase_values = data['phase']

            if not times or (not code_values and not phase_values):
                continue

            epochs = list(range(1, len(times) + 1))
            style = style_dict.get(freq, {'color': 'gray', 'linestyle': '-', 'marker': 'o', 'phase_marker': 'o'})

            # 计算调整常数（统计分析法）
            system = sat_id[0]
            wavelength = self.wavelengths[system].get(freq)
            if wavelength is None:
                continue

            # 过滤掉None值
            valid_indices = []
            for i in range(len(code_values)):
                if code_values[i] is not None and phase_values[i] is not None:
                    valid_indices.append(i)

            if len(valid_indices) < 10:
                print(f"--警告: 卫星 {sat_id} 频率 {freq} 的有效数据点太少，无法计算调整常数")
                continue

            # 创建仅包含有效值的数组
            valid_code_values = np.array([code_values[i] for i in valid_indices])
            valid_phase_values = np.array([phase_values[i] for i in valid_indices])

            # 计算调整常数
            differences = valid_code_values - valid_phase_values * wavelength
            adjustment_constant = np.mean(differences)

            # 调整载波相位值
            adjusted_phase_values = [None] * len(phase_values)
            for i in valid_indices:
                adjusted_phase_values[i] = phase_values[i] * wavelength + adjustment_constant

            # ---------------------- 新增标记间隔控制逻辑 ----------------------
            # 生成每隔200历元的索引
            mark_indices = np.arange(0, len(epochs), 200)
            # ---------------------- 伪距绘图（控制标记显示） ----------------------
            plt.plot(
                epochs, code_values,
                linestyle=style['linestyle'],
                color=style['color'],
                label=style.get('label_code', f'{freq} 伪距 (m)'),
                linewidth=1,  # 保留线宽
                marker=style['marker'],  # 启用标记
                markevery=mark_indices,  # 仅在指定索引显示标记
                markersize=10,  # 适当增大标记尺寸
                markeredgewidth=1  # 标记边框宽度
            )
            has_plotted_lines = True  # 标记有线条绘制
            # ---------------------- 相位绘图（控制标记显示） ----------------------
            plt.plot(
                epochs, adjusted_phase_values,
                linestyle=style['linestyle'],
                color=style['color'],
                label=style.get('label_phase', f'{freq} 相位 (m)'),
                linewidth=1,
                marker=style['phase_marker'],  # 相位单独标记
                markevery=mark_indices,  # 仅在指定索引显示标记
                markersize=10,
                fillstyle='none',  # 空心标记
                markeredgewidth=1
            )
            has_plotted_lines = True  # 标记有线条绘制

        plt.xlabel('历元', fontsize=12)
        plt.ylabel('原始观测值 (米)', fontsize=12)
        plt.title(f'卫星 {sat_id} 原始伪距与载波相位', fontsize=14)
        plt.grid(True, linestyle='--', alpha=0.7)
        if has_plotted_lines:
            plt.legend(loc='upper left', fontsize=10)
        else:
            print(f"--警告: 卫星 {sat_id} 没有可绘制的有效数据，未生成图例")
        plt.tick_params(axis='both', which='major', labelsize=10)
        plt.tight_layout()

        if save:
            category_dir = os.path.join(self.current_result_dir, "raw_observations")
            os.makedirs(category_dir, exist_ok=True)
            filename = f"{sat_id}_raw_observations.png"
            full_path = os.path.join(category_dir, filename)
            plt.savefig(full_path, dpi=600, bbox_inches='tight')
            print(f"--原始观测值图表已保存至: {full_path}")
        else:
            plt.show()
        plt.close()

    def plot_observable_derivatives(self, derivatives: Dict, sat_id: str, freq: str, save=True) -> None:
        """绘制指定卫星和频率的观测值一阶差分变化图"""
        if sat_id in derivatives and freq in derivatives[sat_id]:
            data = derivatives[sat_id][freq]

            system = sat_id[0]
            wavelength = self.wavelengths[system].get(freq)
            if wavelength is None:
                print(f"错误: 无法获取 {sat_id} 的 {freq} 频率波长")
                return

            plt.figure(figsize=(14, 12))

            # 第一部分：将伪距、相位一阶差分和多普勒观测值绘制在一个图中
            plt.subplot(2, 1, 1)

            # 找出所有三种观测值都有效的历元索引
            valid_indices = []
            for i in range(len(data['times'])):
                if (data['pr_derivative'][i] is not None and
                        data['ph_derivative'][i] is not None and
                        data['doppler'][i] is not None):
                    valid_indices.append(i)

            # 使用相同的有效索引提取数据
            valid_pr_derivatives = [data['pr_derivative'][i] for i in valid_indices]
            valid_ph_derivatives = [data['ph_derivative'][i] for i in valid_indices]
            valid_doppler = [data['doppler'][i] for i in valid_indices]

            epochs = list(range(1, len(valid_indices) + 1))

            plt.plot(epochs, valid_pr_derivatives, 'b-', label='伪距一阶差分 (m/s)')
            plt.plot(epochs, valid_ph_derivatives, 'g-', label='相位一阶差分 (m/s)')
            plt.plot(epochs, valid_doppler, 'r-', label='多普勒观测值 (m/s)')

            plt.xlabel('历元')
            plt.ylabel('速度 (m/s)')
            plt.title(f'{sat_id} - {freq} 观测值一阶差分对比')
            plt.grid(True)
            plt.legend()

            # 第二部分：绘制多普勒观测值减去伪距和相位一阶差分后的值
            plt.subplot(2, 1, 2)

            # 计算 -λ·D - dP/dt
            doppler_minus_pr = [
                valid_doppler[i] - valid_pr_derivatives[i]
                for i in range(len(valid_indices))
            ]

            # 计算 -λ·D - dΦ/dt
            doppler_minus_ph = [
                valid_doppler[i] - valid_ph_derivatives[i]
                for i in range(len(valid_indices))
            ]

            plt.plot(epochs, doppler_minus_pr, 'm-', label='-λ·D - dP/dt')
            plt.plot(epochs, doppler_minus_ph, 'c-', label='-λ·D - dΦ/dt')

            plt.axhline(y=0, color='k', linestyle='--', alpha=0.3)
            plt.xlabel('历元')
            plt.ylabel('差值 (m/s)')
            plt.title(f'{sat_id} - {freq} 多普勒与观测值一阶差分差值')
            plt.grid(True)
            plt.legend()

            plt.tight_layout()

            if save:
                category_dir = os.path.join(self.current_result_dir, "derivatives")
                os.makedirs(category_dir, exist_ok=True)
                filename = f"{sat_id}_{freq}_derivatives_comparison.png"
                full_path = os.path.join(category_dir, filename)
                plt.savefig(full_path, dpi=300, bbox_inches='tight')
                print(f"--一阶差分图表已保存至: {full_path}")
            else:
                plt.show()
            plt.close()

    def plot_code_phase_raw_differences(self, differences: Dict, sat_id: str, freq: str, save=True) -> None:
        """绘制指定卫星和频率的伪距与载波相位原始差值图"""
        if sat_id in differences and freq in differences[sat_id]:
            data = differences[sat_id][freq]

            # 获取有效差值数据和对应的时间点
            valid_diffs = []
            valid_times = []
            for i, diff in enumerate(data['code_phase_diff']):
                if diff is not None:
                    valid_diffs.append(diff)
                    valid_times.append(data['times'][i])

            if not valid_diffs:
                print(f"--警告: {sat_id} - {freq} 没有有效的伪距相位差值数据")
                return

            plt.figure(figsize=(12, 6))

            epochs = list(range(1, len(valid_diffs) + 1))

            plt.plot(epochs, valid_diffs, 'b-', label='伪距-相位差')

            plt.xlabel('历元')
            plt.ylabel('差值 (m)')
            plt.title(f'{sat_id} - {freq} 伪距相位差')
            plt.grid(True)
            plt.legend()

            plt.tight_layout()

            if save:
                category_dir = os.path.join(self.current_result_dir, "code_phase_diff_raw")
                os.makedirs(category_dir, exist_ok=True)
                filename = f"{sat_id}_{freq}_code_phase_raw_diff.png"
                full_path = os.path.join(category_dir, filename)
                plt.savefig(full_path, dpi=300, bbox_inches='tight')
                print(f"--码相差值图表已保存至: {full_path}")
            else:
                plt.show()
            plt.close()

    def plot_code_phase_differences(self, differences: Dict, sat_id: str, freq: str, save=True) -> None:
        """绘制指定卫星和频率的伪距与载波相位差值变化图（横坐标为历元序号）"""
        if sat_id in differences and freq in differences[sat_id]:
            data = differences[sat_id][freq]

            # 绘制差值变化率
            valid_changes = [c for c in data['diff_changes'] if c is not None]
            epochs = list(range(1, len(valid_changes) + 1))  # 生成历元序号列表

            plt.figure(figsize=(12, 6))
            plt.plot(epochs, valid_changes, 'b-', label='伪距-相位差变化率')

            # # 添加阈值线（10米）
            # plt.axhline(y=10.0, color='r', linestyle='--', alpha=0.5, label='阈值线 (10米)')

            plt.xlabel('历元')
            plt.ylabel('差值变化 (m/历元)')
            plt.title(f'{sat_id} - {freq} 伪距相位差变化率')
            plt.grid(True)
            plt.legend()

            plt.tight_layout()

            if save:
                category_dir = os.path.join(self.current_result_dir, "code_phase_diffs")
                os.makedirs(category_dir, exist_ok=True)
                filename = f"{sat_id}_{freq}_code_phase_diff_changes.png"
                full_path = os.path.join(category_dir, filename)
                plt.savefig(full_path, dpi=300, bbox_inches='tight')
                print(f"--码相差值变化图表已保存至: {full_path}")
            else:
                plt.show()
            plt.close()

    def plot_phase_prediction_errors(self, errors: Dict, sat_id: str, freq: str, save=True) -> None:
        """绘制指定卫星和频率的载波相位预测误差图（横坐标为历元）"""
        if sat_id in errors and freq in errors[sat_id]:
            data = errors[sat_id][freq]

            plt.figure(figsize=(12, 8))

            # 绘制预测误差
            valid_errors = [e for e in data['prediction_error'] if e is not None]
            epochs = list(range(1, len(valid_errors) + 1))  # 生成历元序号

            plt.plot(epochs, valid_errors, 'b-', label='预测误差')
            plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
            plt.xlabel('历元')
            plt.ylabel('预测误差 (m)')
            plt.title(f'{sat_id} - {freq} 载波相位预测误差')
            plt.grid(True)
            plt.legend()

            plt.tight_layout()

            if save:
                category_dir = os.path.join(self.current_result_dir, "phase_pred_errors")
                os.makedirs(category_dir, exist_ok=True)
                filename = f"{sat_id}_{freq}_phase_pred_error.png"
                full_path = os.path.join(category_dir, filename)
                plt.savefig(full_path, dpi=300, bbox_inches='tight')
                print(f"--相位预测图表已保存至: {full_path}")
            else:
                plt.show()
            plt.close()

    # def plot_double_differences(self, double_diffs, triple_errors, sat_id, freq, save=True):
    #     """绘制所有观测值双差结果并显示阈值"""
    #     if sat_id not in double_diffs or freq not in double_diffs[sat_id]:
    #         print(f"错误: 卫星 {sat_id} 或频率 {freq} 无双差数据")
    #         return
    #
    #     data = double_diffs[sat_id][freq]
    #     errors = triple_errors[sat_id][freq]
    #     epochs = list(range(1, len(data['dd_code']) + 1))  # 历元序号
    #
    #     # 创建子图
    #     fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10))
    #
    #     # 1. 伪距双差
    #     ax1.plot(epochs, data['dd_code'], 'b-', label='伪距双差 (m)')
    #     ax1.axhline(y=errors['code']['threshold'], color='r', linestyle='--', alpha=0.7,
    #                 label=f'阈值: {errors["code"]["threshold"]:.2f} m')
    #     ax1.axhline(y=-errors['code']['threshold'], color='r', linestyle='--', alpha=0.7)
    #     ax1.set_xlabel('历元')
    #     ax1.set_ylabel('伪距双差 (m)')
    #     ax1.set_title(f'{sat_id} - {freq} 历元间双差')
    #     ax1.grid(True)
    #     ax1.legend()
    #
    #     # 2. 载波相位双差
    #     ax2.plot(epochs, data['dd_phase'], 'g-', label='载波相位双差 (m)')
    #     ax2.axhline(y=errors['phase']['threshold'], color='r', linestyle='--', alpha=0.7,
    #                 label=f'阈值: {errors["phase"]["threshold"]:.2f} m')
    #     ax2.axhline(y=-errors['phase']['threshold'], color='r', linestyle='--', alpha=0.7)
    #     ax2.set_xlabel('历元')
    #     ax2.set_ylabel('载波相位双差 (m)')
    #     ax2.grid(True)
    #     ax2.legend()
    #
    #     # 3. 多普勒双差
    #     ax3.plot(epochs, data['dd_doppler'], 'm-', label='多普勒双差 (cycle/s)')
    #     ax3.axhline(y=errors['doppler']['threshold'], color='r', linestyle='--', alpha=0.7,
    #                 label=f'阈值: {errors["doppler"]["threshold"]:.2f} cycle/s')
    #     ax3.axhline(y=-errors['doppler']['threshold'], color='r', linestyle='--', alpha=0.7)
    #     ax3.set_xlabel('历元')
    #     ax3.set_ylabel('多普勒双差 (cycle/s)')
    #     ax3.grid(True)
    #     ax3.legend()
    #
    #     plt.tight_layout()
    #
    #     if save:
    #         category_dir = os.path.join(self.current_result_dir, "double_differences")
    #         os.makedirs(category_dir, exist_ok=True)
    #         filename = f"{sat_id}_{freq}_double_differences.png"
    #         full_path = os.path.join(category_dir, filename)
    #         plt.savefig(full_path, dpi=300, bbox_inches='tight')
    #         print(f"双差图表已保存至: {full_path}")
    #     else:
    #         plt.show()
    #     plt.close()
    def plot_double_differences(self, double_diffs, triple_errors, sat_id, freq, save=True):
        """绘制伪距、相位和多普勒双差结果并显示各自的阈值"""
        if sat_id not in double_diffs or freq not in double_diffs[sat_id]:
            print(f"错误: 卫星 {sat_id} 或频率 {freq} 无双差数据")
            return

        data = double_diffs[sat_id][freq]
        errors = triple_errors[sat_id][freq]
        epochs = list(range(1, len(data['dd_code']) + 1))  # 历元序号

        # 创建包含三个子图的图表，修改 sharex 参数为 'all'
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 15), sharex='all')
        fig.suptitle(f'{sat_id} - {freq} 历元间双差分析', fontsize=14)

        # 1. 伪距双差
        ax1.plot(epochs, data['dd_code'], 'b-', label='伪距双差 (m)')
        ax1.axhline(y=errors['code']['threshold'], color='r', linestyle='--', alpha=0.7,
                    label=f'阈值: {errors["code"]["threshold"]:.2f} m')
        ax1.axhline(y=-errors['code']['threshold'], color='r', linestyle='--', alpha=0.7)
        ax1.set_ylabel('伪距双差 (m)')
        ax1.grid(True, linestyle='--', alpha=0.7)
        ax1.legend(loc='upper right')
        ax1.set_title('伪距双差', fontsize=12)

        # 2. 载波相位双差
        ax2.plot(epochs, data['dd_phase'], 'g-', label='相位双差 (m)')
        ax2.axhline(y=errors['phase']['threshold'], color='r', linestyle='--', alpha=0.7,
                    label=f'阈值: {errors["phase"]["threshold"]:.2f} m')
        ax2.axhline(y=-errors['phase']['threshold'], color='r', linestyle='--', alpha=0.7)
        ax2.set_ylabel('相位双差 (m)')
        ax2.grid(True, linestyle='--', alpha=0.7)
        ax2.legend(loc='upper right')
        ax2.set_title('载波相位双差', fontsize=12)

        # 3. 多普勒双差
        ax3.plot(epochs, data['dd_doppler'], 'm-', label='多普勒双差 (m/s)')
        ax3.axhline(y=errors['doppler']['threshold'], color='r', linestyle='--', alpha=0.7,
                    label=f'阈值: {errors["doppler"]["threshold"]:.2f} m/s')
        ax3.axhline(y=-errors['doppler']['threshold'], color='r', linestyle='--', alpha=0.7)
        ax3.set_xlabel('历元')
        ax3.set_ylabel('多普勒双差 (m/s)')
        ax3.grid(True, linestyle='--', alpha=0.7)
        ax3.legend(loc='upper right')
        ax3.set_title('多普勒双差', fontsize=12)

        # 调整子图间距
        plt.tight_layout(rect=[0, 0, 1, 0.96])

        if save:
            category_dir = os.path.join(self.current_result_dir, "double_differences")
            os.makedirs(category_dir, exist_ok=True)
            filename = f"{sat_id}_{freq}_double_differences.png"
            full_path = os.path.join(category_dir, filename)
            plt.savefig(full_path, dpi=300, bbox_inches='tight')
            print(f"--双差图表已保存至: {full_path}")
        else:
            plt.show()
        plt.close()

    def save_report(self) -> None:
        """保存分析报告到当前文件结果目录"""
        self.start_stage(7, "保存分析报告", 100)

        report = self.generate_report()
        filename = "analysis_report.txt"
        full_path = os.path.join(self.current_result_dir, filename)

        with open(full_path, 'w', encoding='utf-8') as f:
            f.write(report)

        self.update_progress(100)
        print(f"--报告已保存至: {full_path}")

    def save_all_plots(self) -> None:
        """保存所有分析结果的图表"""
        self.start_stage(8, "保存分析图表", 100)
        print("--开始保存各图表...")
        # 创建结果目录
        os.makedirs(self.current_result_dir, exist_ok=True)
        # 创建类别子文件夹
        for category in self.plot_categories:
            os.makedirs(os.path.join(self.current_result_dir, category), exist_ok=True)

        # 保存原始观测值图表
        if self.observations_meters:
            total_sats = len(self.observations_meters)
            processed_sats = 0

            for sat_id in self.observations_meters:
                has_data = False
                for freq, data in self.observations_meters[sat_id].items():
                    if (len([c for c in data['code'] if c is not None]) > 0 or
                            len([p for p in data['phase'] if p is not None]) > 0):
                        has_data = True
                        break
                if has_data:
                    self.plot_raw_observations(sat_id, save=True)
                    processed_sats += 1
                    self.update_progress(int(processed_sats / total_sats * 16))
                else:
                    print(f"--跳过保存: {sat_id} 原始观测数据不足")

        # 保存观测值一阶差分图表
        if 'observable_derivatives' in self.results:
            total_sats = len(self.results['observable_derivatives'])
            processed_sats = 0

            for sat_id in self.results['observable_derivatives']:
                has_data = False
                for freq, data in self.results['observable_derivatives'][sat_id].items():
                    valid_pr = [d for d in data['pr_derivative'] if d is not None]
                    valid_ph = [d for d in data['ph_derivative'] if d is not None]
                    valid_doppler = [d for d in data['doppler'] if d is not None]
                    if len(valid_pr) > 0 and len(valid_ph) > 0 and len(valid_doppler) > 0:
                        has_data = True
                        break
                if has_data:
                    for freq in self.results['observable_derivatives'][sat_id]:
                        self.plot_observable_derivatives(
                            self.results['observable_derivatives'],
                            sat_id, freq, save=True
                        )
                    processed_sats += 1
                    self.update_progress(16 + int(processed_sats / total_sats * 16))
                else:
                    print(f"--跳过保存: {sat_id} 观测值一阶差分数据不足")

        # 保存伪距相位原始差值图表
        if 'code_phase_differences' in self.results:
            total_sats = len(self.results['code_phase_differences'])
            processed_sats = 0

            for sat_id in self.results['code_phase_differences']:
                has_data = False
                for freq, data in self.results['code_phase_differences'][sat_id].items():
                    valid_diffs = [d for d in data['code_phase_diff'] if d is not None]
                    if len(valid_diffs) > 0:
                        has_data = True
                        break
                if has_data:
                    for freq in self.results['code_phase_differences'][sat_id]:
                        # 调用新的绘图函数
                        self.plot_code_phase_raw_differences(
                            self.results['code_phase_differences'],
                            sat_id, freq, save=True
                        )
                    processed_sats += 1
                    self.update_progress(32 + int(processed_sats / total_sats * 16))
                else:
                    print(f"--跳过保存: {sat_id} 伪距相位原始差值数据不足")

        # 保存伪距相位差值变化图表
        if 'code_phase_differences' in self.results:
            total_sats = len(self.results['code_phase_differences'])
            processed_sats = 0

            for sat_id in self.results['code_phase_differences']:
                has_data = False
                for freq, data in self.results['code_phase_differences'][sat_id].items():
                    valid_diffs = [d for d in data['code_phase_diff'] if d is not None]
                    if len(valid_diffs) > 0:
                        has_data = True
                        break
                if has_data:
                    for freq in self.results['code_phase_differences'][sat_id]:
                        self.plot_code_phase_differences(
                            self.results['code_phase_differences'],
                            sat_id, freq, save=True
                        )
                    processed_sats += 1
                    self.update_progress(47 + int(processed_sats / total_sats * 16))
                else:
                    print(f"--跳过保存: {sat_id} 伪距相位差值数据不足")

        # 保存相位预测误差图表
        if 'phase_prediction_errors' in self.results:
            total_sats = len(self.results['phase_prediction_errors'])
            processed_sats = 0

            for sat_id in self.results['phase_prediction_errors']:
                has_data = False
                for freq, data in self.results['phase_prediction_errors'][sat_id].items():
                    valid_errors = [e for e in data['prediction_error'] if e is not None]
                    if len(valid_errors) > 0:
                        has_data = True
                        break
                if has_data:
                    for freq in self.results['phase_prediction_errors'][sat_id]:
                        self.plot_phase_prediction_errors(
                            self.results['phase_prediction_errors'],
                            sat_id, freq, save=True
                        )
                    processed_sats += 1
                    self.update_progress(63 + int(processed_sats / total_sats * 16))
                else:
                    print(f"--跳过保存: {sat_id} 相位预测误差数据不足")

        # 保存双差图表
        if 'double_differences' in self.results and 'triple_median_errors' in self.results:
            double_diffs = self.results['double_differences']
            triple_errors = self.results['triple_median_errors']
            total_sats = len(double_diffs)
            processed_sats = 0

            for sat_id in double_diffs:
                has_data = any(len(freq_data['dd_code']) > 0 for freq_data in double_diffs[sat_id].values())
                if has_data:
                    for freq in double_diffs[sat_id]:
                        self.plot_double_differences(double_diffs, triple_errors, sat_id, freq, save=True)
                    processed_sats += 1
                    self.update_progress(79 + int(processed_sats / total_sats * 21))
                else:
                    print(f"--跳过保存: {sat_id} 双差数据不足")

        # 完成进度
        self.update_progress(100)
        print("所有图表保存完成")

    def generate_report(self) -> str:
        """生成检测结果报告"""
        report = "=== GNSS观测数据分析报告 ===\n\n"

        # 观测值一阶差分变化与多普勒值不一致性检测
        report += "1. 观测值一阶差分变化与多普勒值不一致性检测:\n"
        inconsistent_derivatives_doppler = {}
        for sat_id, freq_data in self.results['observable_derivatives'].items():
            for freq, data in freq_data.items():
                valid_indices = []
                for i in range(len(data['times'])):
                    if (data['pr_derivative'][i] is not None and
                            data['ph_derivative'][i] is not None and
                            data['doppler'][i] is not None):
                        valid_indices.append(i)
                valid_pr_derivatives = [data['pr_derivative'][i] for i in valid_indices]
                valid_ph_derivatives = [data['ph_derivative'][i] * self.wavelengths[sat_id[0]].get(freq) for i in
                                        valid_indices]
                valid_doppler = [-data['doppler'][i] * self.wavelengths[sat_id[0]].get(freq) for i in valid_indices]

                doppler_minus_pr = [valid_doppler[i] - valid_pr_derivatives[i] for i in range(len(valid_indices))]
                doppler_minus_ph = [valid_doppler[i] - valid_ph_derivatives[i] for i in range(len(valid_indices))]

                max_diff_pr = max(doppler_minus_pr) if doppler_minus_pr else 0
                max_diff_ph = max(doppler_minus_ph) if doppler_minus_ph else 0
                threshold_pr = 10  # 可根据实际情况调整阈值，单位m/s
                threshold_ph = 10  # 可根据实际情况调整阈值，单位m/s
                if max_diff_pr > threshold_pr or max_diff_ph > threshold_ph:
                    inconsistent_derivatives_doppler.setdefault(sat_id, {})[freq] = {
                        'max_diff_pr': max_diff_pr,
                        'max_diff_ph': max_diff_ph,
                        'threshold_pr': threshold_pr,
                        'threshold_ph': threshold_ph
                    }

        if inconsistent_derivatives_doppler:
            report += f"  发现 {len(inconsistent_derivatives_doppler)} 颗卫星在部分频率上存在观测值一阶差分变化与多普勒值不一致性:\n"
            for sat_id, freq_info in inconsistent_derivatives_doppler.items():
                for freq, info in freq_info.items():
                    report += f"    - 卫星 {sat_id}，频率 {freq}:\n"
                    report += f"      - 多普勒值减去伪距一阶差分之最大差值: {info['max_diff_pr']:.2f} m/s (阈值: {info['threshold_pr']} m/s)\n"
                    report += f"      - 多普勒值减去相位一阶差分之最大差值: {info['max_diff_ph']:.2f} m/s (阈值: {info['threshold_ph']} m/s)\n"
                    report += "      - 潜在影响: 可能表明卫星信号存在异常，或在信号传播过程中受到干扰，导致观测值与理论计算值偏差较大，影响定位精度。\n"
        else:
            report += "  未发现观测值一阶差分变化与多普勒值不一致性问题\n"
        report += "\n"

        # 伪距相位差值不一致检测
        report += "2. 伪距相位差值不一致检测:\n"
        # 设置阈值(单位:米)
        threshold = 10.0  # 差值变化阈值（米）
        inconsistent_code_phase = []
        for sat, freq_data in self.results['code_phase_differences'].items():
            for freq, data in freq_data.items():
                # 使用diff_changes而非code_phase_diff
                changes = [c for c in data['diff_changes'] if c is not None]
                if len(changes) < 3:  # 数据点太少，跳过
                    continue

                # 计算变化率的统计特征
                mean_change = np.mean(changes)
                max_change = np.max(changes)
                std_change = np.std(changes)

                # 数据过滤信息
                original_epochs = data.get('original_epochs', 0)
                filtered_epochs = data.get('filtered_epochs', 0)
                stagnant_removed = data.get('stagnant_epochs_removed', 0)
                missing_obs = data.get('missing_epochs', 0)

                # 判断是否异常
                if max_change > threshold and std_change > 1.0:
                    inconsistent_code_phase.append({
                        'sat': sat,
                        'freq': freq,
                        'max_change': max_change,
                        'mean_change': mean_change,
                        'std_change': std_change,
                        'threshold': threshold,
                        'original_epochs': original_epochs,
                        'filtered_epochs': filtered_epochs,
                        'stagnant_removed': stagnant_removed,
                        'missing_obs': missing_obs
                    })

        if inconsistent_code_phase:
            report += f"  发现 {len(inconsistent_code_phase)} 颗卫星存在伪距相位差值波动异常:\n"
            for item in inconsistent_code_phase:
                report += f"    - 卫星 {item['sat']}，频率 {item['freq']}:\n"
                report += f"      - 差值变化统计: 最大值={item['max_change']:.3f} m, 均值={item['mean_change']:.3f} m, 标准差={item['std_change']:.3f} m\n"
                report += f"      - 阈值: {item['threshold']} m\n"
                report += f"      - 数据过滤: 原始历元={item['original_epochs']}, 缺失观测值={item['missing_obs']}, 剔除停滞历元={item['stagnant_removed']}, 有效历元={item['filtered_epochs']}\n"

                # 分析差值异常可能原因
                if item['missing_obs'] > item['filtered_epochs'] * 0.3:  # 缺失率超过30%
                    report += "      - 潜在原因: 大量观测值缺失导致计算不稳定，可能是信号遮挡或接收机问题\n"
                elif item['max_change'] > 1000:  # 异常大的变化值
                    report += "      - 潜在原因: 整周模糊度变化、信号失锁或接收机钟跳变\n"
                else:
                    report += "      - 潜在原因: 多径效应、接收机噪声或卫星钟差异常\n"
                report += "      - 潜在影响: 可能导致定位结果偏差，影响导航和定位服务的可靠性\n"
        else:
            report += "  未发现伪距相位差值波动异常（所有波动均在10米以内）\n"
        report += "\n"

        # 相位停滞统计
        report += "\n3. 相位停滞统计:\n"

        if 'phase_stagnation' in self.results:
            stagnant_sats = set()
            total_stagnant_freqs = 0
            for sat_id, freq_data in self.results['phase_stagnation'].items():
                for freq, stats in freq_data.items():
                    if stats.get('is_stagnant', False):
                        stagnant_sats.add(sat_id)
                        total_stagnant_freqs += 1

            report += f"  检测到 {len(stagnant_sats)} 颗卫星存在相位停滞现象，涉及 {total_stagnant_freqs} 个频率通道\n"
            if stagnant_sats:
                report += "  停滞卫星列表:\n"
                for sat in stagnant_sats:
                    sat_freqs = [freq for freq, stats in self.results['phase_stagnation'][sat].items()
                                 if stats.get('is_stagnant', False)]
                    report += f"    - {sat}: 频率 {', '.join(sat_freqs)}\n"
            else:
                report += "  未检测到相位停滞现象\n"
        else:
            report += "  未执行相位停滞检测\n"

        # 粗差检测：历元间双差超限检测
        report += "\n4. 历元间双差超限检测:\n"
        double_diff_outliers = self.results.get('triple_median_errors', {})  # 从results中获取双差误差结果
        total_outliers = 0

        for sat_id, freq_data in double_diff_outliers.items():
            for freq, errors in freq_data.items():
                code_outliers = errors['code']['outliers']  # 伪距双差超限历元索引
                phase_outliers = errors['phase']['outliers']  # 载波相位双差超限历元索引
                doppler_outliers = errors['doppler']['outliers']  # 多普勒双差超限历元索引

                if code_outliers or phase_outliers or doppler_outliers:
                    total_outliers += 1
                    report += f"  - 卫星 {sat_id}, 频率 {freq}:\n"
                    if code_outliers:
                        report += f"    ▶ 伪距双差超限历元: {code_outliers}（阈值: {errors['code']['threshold']:.2f} m）\n"
                    if phase_outliers:
                        report += f"    ▶ 载波相位双差超限历元: {phase_outliers}（阈值: {errors['phase']['threshold']:.2f} m）\n"
                    if doppler_outliers:
                        report += f"    ▶ 多普勒双差超限历元: {doppler_outliers}（阈值: {errors['doppler']['threshold']:.2f} m/s）\n"
                    report += "    - 潜在原因: 多径效应、卫星信号失锁或接收机钟差异常\n"
                    report += "    - 潜在影响: 可能导致载波相位差分（RTK）解算失败或定位精度骤降\n"

        if total_outliers == 0:
            report += "  未检测到双差超限值（所有双差值均在3×MAD阈值内）\n"
        else:
            report += f"  共检测到 {total_outliers} 个频率通道存在双差超限值，建议结合图表分析具体历元\n"

        return report

    def reset_analysis(self):
        """重置分析器状态，准备处理新文件"""
        self.observations_meters = {}
        self.results = {
            'code_carrier_inconsistency': {},
            'observation_inconsistency': {},
            'phase_stagnation': {},
            'observable_derivatives': {},
            'code_phase_differences': {},
            'phase_prediction_errors': {}
        }
        # 重置当前文件相关路径
        self.input_file_path = None
        self.current_result_dir = None


# 主函数
def main():
    root = tk.Tk()
    root.withdraw()

    # 构建文件类型筛选器（.yyO/.yyo和.RNX/.rnx格式）
    file_types = [
        ("RINEX Files", "*.??O *.??o *.RNX *.rnx"),
        ("All Files", "*.*")
    ]

    # 调用多选文件对话框
    selected_files = filedialog.askopenfilenames(
        title="选择RINEX观测文件（支持多选）",
        filetypes=file_types
    )

    if not selected_files:
        print("未选择任何文件")
        return

    # 转换为列表并验证文件格式
    file_paths = list(selected_files)
    valid_files = []
    pattern = r'^.*?\.\d{2}[oO]$|^.*?\.(RNX|rnx)$'  # 匹配.yyO/.yyo或.RNX/.rnx
    for file_path in file_paths:
        if re.search(pattern, file_path, re.IGNORECASE):
            valid_files.append(file_path)
        else:
            print(f"警告: 跳过无效文件 {os.path.basename(file_path)}（需为RINEX格式）")

    if not valid_files:
        print("没有有效的RINEX文件被选择")
        return

    # 进度窗口
    progress_window = tk.Toplevel()
    progress_window.title("GNSS数据分析进度")
    progress_window.geometry("400x250")
    progress_window.resizable(False, False)

    # 进度条
    progress_frame = ttk.Frame(progress_window, padding="20")
    progress_frame.pack(fill=tk.BOTH, expand=True)

    ttk.Label(progress_frame, text="处理进度:").pack(anchor=tk.W)

    # 总体进度条
    total_progress_var = tk.DoubleVar()
    total_progress_label = ttk.Label(progress_frame, text="总体进度:")
    total_progress_label.pack(anchor=tk.W)
    total_progress_bar = ttk.Progressbar(progress_frame, orient=tk.HORIZONTAL,
                                         length=360, variable=total_progress_var, mode='determinate')
    total_progress_bar.pack(pady=5)

    # 当前文件进度条
    file_progress_var = tk.DoubleVar()
    file_progress_label = ttk.Label(progress_frame, text="当前文件进度:")
    file_progress_label.pack(anchor=tk.W)
    file_progress_bar = ttk.Progressbar(progress_frame, orient=tk.HORIZONTAL,
                                        length=360, variable=file_progress_var, mode='determinate')
    file_progress_bar.pack(pady=5)

    status_var = tk.StringVar()
    status_var.set(f"准备处理 {len(valid_files)} 个文件...")
    ttk.Label(progress_frame, textvariable=status_var).pack(anchor=tk.W)

    current_file_var = tk.StringVar()
    current_file_var.set("当前文件: 无")
    ttk.Label(progress_frame, textvariable=current_file_var).pack(anchor=tk.W)

    # 进度更新函数
    def update_file_progress(progress):
        file_progress_var.set(progress * 100)
        progress_window.update_idletasks()

    def update_total_progress(current_file, total_files, progress):
        total_progress = ((current_file - 1) + progress) / total_files
        total_progress_var.set(total_progress * 100)
        status_var.set(f"处理中: 文件 {current_file}/{total_files} ({int(total_progress * 100)}%)")
        progress_window.update_idletasks()

    # 处理函数（在单独线程中运行）
    def process_data():
        try:
            analyzer = GNSSAnalyzer()
            analyzer.set_progress_callback(update_file_progress)

            total_files = len(valid_files)

            for i, file_path in enumerate(valid_files, 1):
                current_file_var.set(f"当前文件: {os.path.basename(file_path)}")
                update_total_progress(i, total_files, 0)

                # 重置分析器状态
                analyzer.reset_analysis()

                # 读取文件
                status_var.set(f"正在读取RINEX文件 ({i}/{total_files})...")
                data = analyzer.read_rinex_obs(file_path)
                print(f"--成功读取文件: {file_path}")

                # 计算新增分析指标
                status_var.set(f"正在计算观测值一阶差分 ({i}/{total_files})...")
                derivatives = analyzer.calculate_observable_derivatives(data)

                status_var.set(f"正在计算伪距相位差值 ({i}/{total_files})...")
                code_phase_diffs = analyzer.calculate_code_phase_differences(data)
                code_phase_raw_diffs = analyzer.calculate_code_phase_differences(data)

                status_var.set(f"正在计算相位预测误差 ({i}/{total_files})...")
                phase_pred_errors = analyzer.calculate_phase_prediction_errors(data)

                status_var.set(f"正在计算历元间双差 ({i}/{total_files})...")
                double_diffs = analyzer.calculate_epoch_double_differences()  # 计算双差

                status_var.set(f"正在计算三倍中误差 ({i}/{total_files})...")
                triple_errors = analyzer.calculate_triple_median_error(double_diffs)  # 计算三倍中误差

                status_var.set(f"正在剔除异常值并保存文件 ({i}/{total_files})...")
                analyzer.remove_outliers_and_save(double_diffs, triple_errors)  # 新增：剔除异常值并保存

                # 保存报告和所有图表
                status_var.set(f"正在保存分析报告 ({i}/{total_files})...")
                analyzer.save_report()

                status_var.set(f"正在保存分析图表 ({i}/{total_files})...")
                analyzer.save_all_plots()

                update_total_progress(i, total_files, 1)

                # # 如果是最后一个文件，显示示例图表
                # if i == total_files:
                #     # 显示某卫星某频率图表示例
                #     sample_sat = "C40"
                #     sample_freq = "L1P"
                #     if sample_sat in analyzer.observations_meters:
                #         analyzer.plot_raw_observations(sample_sat, save=False)
                #     if sample_sat in derivatives and sample_freq in derivatives[sample_sat]:
                #         analyzer.plot_observable_derivatives(derivatives, sample_sat, sample_freq, save=False)
                #     if sample_sat in code_phase_raw_diffs and sample_freq in code_phase_raw_diffs[sample_sat]:
                #         analyzer.plot_code_phase_raw_differences(code_phase_raw_diffs, sample_sat, sample_freq,
                #                                                  save=False)
                #     if sample_sat in code_phase_diffs and sample_freq in code_phase_diffs[sample_sat]:
                #         analyzer.plot_code_phase_differences(code_phase_diffs, sample_sat, sample_freq, save=False)
                #     if sample_sat in phase_pred_errors and sample_freq in phase_pred_errors[sample_sat]:
                #         analyzer.plot_phase_prediction_errors(phase_pred_errors, sample_sat, sample_freq, save=False)

                # 处理完成后关闭进度窗口
                status_var.set(f"处理完成! 共处理 {total_files} 个文件")
                current_file_var.set("")
                progress_window.after(2000, progress_window.destroy)
                root.after(3000, root.quit)

        except Exception as e:
            import traceback
            error_info = {
                "错误类型": str(type(e).__name__),
                "错误信息": str(e),
                "错误位置": traceback.format_exc(),
                "当前处理阶段": status_var.get(),
                "当前文件": current_file_var.get()
            }

            # 在控制台打印详细错误信息
            print("\n==== 详细错误信息 ====")
            for key, value in error_info.items():
                print(f"{key}: {value}")
            print("====================\n")

            # 在GUI中显示简化的错误信息
            status_var.set(f"错误: {error_info['错误类型']} - {error_info['错误信息']}")

            # 保存完整的错误日志到文件
            error_log_path = os.path.join(os.path.dirname(file_paths[0]), "gnss_analyzer_error.log")
            with open(error_log_path, 'w', encoding='utf-8') as f:
                f.write("==== GNSS分析器错误日志 ====\n\n")
                for key, value in error_info.items():
                    f.write(f"{key}:\n{value}\n\n")
                f.write("==== 环境信息 ====\n")
                f.write(f"Python版本: {sys.version}\n")
                f.write(f"运行时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

            print(f"完整错误日志已保存至: {error_log_path}")

    processing_thread = threading.Thread(target=process_data)
    processing_thread.daemon = True
    processing_thread.start()
    root.mainloop()


if __name__ == "__main__":
    main()
