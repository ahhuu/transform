import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, least_squares
import os
import sys


# 定义模型函数
def elevation_model(elevation_rad, a, b):
    """高度角模型：σ² = a + b / sin²(el)"""
    sin_el = np.sin(elevation_rad)
    sin_el = np.clip(sin_el, 0.1, 1.0)  # 避免sin(el)过小导致数值问题
    return a + b / (sin_el ** 2)


def snr_model(snr, a, b):
    """信噪比模型：σ² = a + b × 10^(-SNR/10)"""
    return a + b * 10 ** (-snr / 10)


def linear_snr_model(snr, a, b):
    """线性信噪比模型：σ² = a × SNR + b"""
    return a * snr + b


def combined_model(params, elevation_rad, snr):
    """联合模型：σ² = a + b / sin²(el) + c × 10^(-SNR/10)"""
    a, b, c = params
    sin_el = np.sin(elevation_rad)
    sin_el = np.clip(sin_el, 0.1, 1.0)
    return a + b / (sin_el ** 2) + c * 10 ** (-snr / 10)


def error_function(params, model, x, y):
    """联合模型的误差函数"""
    return model(params, *x) - y ** 2


# 模型拟合函数
def fit_elevation_model(x, y):
    p0 = [1, 1]
    return curve_fit(elevation_model, x, y ** 2, p0=p0)[0]


def fit_snr_model(x, y):
    p0 = [1, 100]
    return curve_fit(snr_model, x, y ** 2, p0=p0)[0]


def fit_linear_snr_model(x, y):
    p0 = [0.1, 10]
    return curve_fit(linear_snr_model, x, y ** 2, p0=p0)[0]


def fit_combined_model(x, y):
    x_comb = (x[0], x[1])
    p0 = [1, 1, 100]
    return least_squares(error_function, p0, args=(combined_model, x_comb, y)).x


# 生成模型信息文本函数
def generate_model_info(system, freq, data_count,
                        elev_params, snr_params, linear_snr_params, combined_params):
    model_info = f"""
卫星系统: {system}, 频率: {freq}, 数据量: {data_count}

高度角模型：
σ² = {elev_params[0]:.4f}+{elev_params[1]:.4f}./sin(e).^2

指数信噪比模型：
σ² = {snr_params[0]:.4f}+{snr_params[1]:.4f}.*10.^(-SNR/10)

线性信噪比模型：
σ² = {linear_snr_params[0]:.4f}.*SNR+{linear_snr_params[1]:.4f}

联合模型：
σ² = {combined_params[0]:.4f}+{combined_params[1]:.4f}./sin(e).^2+{combined_params[2]:.4f}.*10.^(-SNR/10)
"""
    return model_info


# 可视化函数
def visualize_elevation_model(ax, subset, a_elev, b_elev, system, freq, data_count):
    ax.scatter(subset['ele'], subset['residual'].values ** 2, alpha=0.5, s=10, label='观测值')
    elev_range = np.linspace(min(subset['ele']), max(subset['ele']), 100)
    elev_range_rad = np.radians(elev_range)
    ax.plot(elev_range, elevation_model(elev_range_rad, a_elev, b_elev), 'r-', label='拟合曲线')
    ax.set_xlabel('高度角 (°)')
    ax.set_ylabel('伪距残差平方 (m²)')
    ax.set_title('高度角模型')
    ax.legend()


def visualize_snr_model(ax, subset, a_snr, b_snr, system, freq, data_count):
    ax.scatter(subset['snr'], subset['residual'].values ** 2, alpha=0.5, s=10, label='观测值')
    snr_range = np.linspace(min(subset['snr']), max(subset['snr']), 100)
    ax.plot(snr_range, snr_model(snr_range, a_snr, b_snr), 'g-', label='指数模型')
    ax.set_xlabel('信噪比 (dB-Hz)')
    ax.set_ylabel('伪距残差平方 (m²)')
    ax.set_title('指数信噪比模型')
    ax.legend()


def visualize_linear_snr_model(ax, subset, a_linear_snr, b_linear_snr, system, freq, data_count):
    ax.scatter(subset['snr'], subset['residual'].values ** 2, alpha=0.5, s=10, label='观测值')
    snr_range = np.linspace(min(subset['snr']), max(subset['snr']), 100)
    ax.plot(snr_range, linear_snr_model(snr_range, a_linear_snr, b_linear_snr), 'm-', label='线性模型')
    ax.set_xlabel('信噪比 (dB-Hz)')
    ax.set_ylabel('伪距残差平方 (m²)')
    ax.set_title('线性信噪比模型')
    ax.legend()


def visualize_combined_model(ax, subset, a_comb, b_comb, c_comb, system, freq, data_count):
    snr_mean = subset['snr'].mean()
    ax.scatter(subset['ele'], subset['residual'].values ** 2, alpha=0.5, s=10, label='观测值')
    elev_range = np.linspace(min(subset['ele']), max(subset['ele']), 100)
    elev_range_rad = np.radians(elev_range)
    ax.plot(elev_range, combined_model([a_comb, b_comb, c_comb], elev_range_rad, snr_mean), 'b-',
            label=f'拟合曲线(SNR={snr_mean:.1f}dB)')
    ax.set_xlabel('高度角 (°)')
    ax.set_ylabel('伪距残差平方 (m²)')
    ax.set_title('联合模型')
    ax.legend()


def fit_and_visualize(subset, system, freq, output_dir):
    y = subset['residual'].values
    x_elev = subset['elevation_rad'].values
    x_snr = subset['snr'].values

    # 模型拟合
    a_elev, b_elev = fit_elevation_model(x_elev, y)
    a_snr, b_snr = fit_snr_model(x_snr, y)
    a_linear_snr, b_linear_snr = fit_linear_snr_model(x_snr, y)
    a_comb, b_comb, c_comb = fit_combined_model((x_elev, x_snr), y)

    # 保存结果
    result = {
        'system': system,
        'frequency': freq,
        'elevation_model': (a_elev, b_elev),
        'snr_model': (a_snr, b_snr),
        'linear_snr_model': (a_linear_snr, b_linear_snr),
        'combined_model': (a_comb, b_comb, c_comb),
        'data_count': len(subset)
    }

    # 生成模型信息文本
    model_info = generate_model_info(system, freq, len(subset),
                                     (a_elev, b_elev), (a_snr, b_snr), (a_linear_snr, b_linear_snr),
                                     (a_comb, b_comb, c_comb))
    model_file = os.path.join(output_dir, f'model_{system}_{freq}.txt')
    with open(model_file, 'w', encoding='utf-8') as f:
        f.write(model_info)

    # 可视化拟合效果
    fig, axs = plt.subplots(1, 4, figsize=(20, 5))
    fig.suptitle(f'卫星系统: {system}, 频率: {freq}, 数据量: {len(subset)}')

    visualize_elevation_model(axs[0], subset, a_elev, b_elev, system, freq, len(subset))
    visualize_snr_model(axs[1], subset, a_snr, b_snr, system, freq, len(subset))
    visualize_linear_snr_model(axs[2], subset, a_linear_snr, b_linear_snr, system, freq, len(subset))
    visualize_combined_model(axs[3], subset, a_comb, b_comb, c_comb, system, freq, len(subset))

    plt.tight_layout()

    image_file = os.path.join(output_dir, f'fitting_{system}_{freq}.png')
    fig.savefig(image_file, dpi=300, bbox_inches='tight')
    plt.close(fig)

    return result


def main():
    # 输入文件路径
    input_file = 'results/K70128H/pseudorange_residuals.csv'

    # 创建输出文件夹
    output_dir = os.path.join(os.path.dirname(input_file), 'Weighting')
    os.makedirs(output_dir, exist_ok=True)

    # 配置 matplotlib
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
    plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题

    # 读取数据
    df = pd.read_csv(input_file, parse_dates=['epoch'])
    df['prn'] = df['prn'].astype(str)  # 确保PRN列为字符串类型
    df['system'] = df['prn'].str[0]  # 提取卫星系统信息

    # 数据预处理
    df = df.dropna()
    df = df[df['residual'] != 0]
    df['ele'] = pd.to_numeric(df['ele'], errors='coerce')
    df['snr'] = pd.to_numeric(df['snr'], errors='coerce')
    df['residual'] = pd.to_numeric(df['residual'], errors='coerce')
    df = df.dropna()
    df['elevation_rad'] = np.radians(df['ele'])

    # 按卫星系统和频率分组拟合
    systems = df['system'].unique()
    frequencies = df['frequency'].unique()
    results = []
    print("开始分组参数拟合...")
    for system in systems:
        for freq in frequencies:
            subset = df[(df['system'] == system) & (df['frequency'] == freq)]
            if len(subset) < 10:  # 数据量太少则不拟合
                continue
            result = fit_and_visualize(subset, system, freq, output_dir)
            results.append(result)

    # 所有卫星系统和频率一起进行全局拟合
    global_subset = df.copy()
    if len(global_subset) >= 50:  # 确保数据量足够
        print("开始全局参数拟合...")
        result = fit_and_visualize(global_subset, 'ALL', 'ALL', output_dir)
        result['system'] = 'ALL'
        result['frequency'] = 'ALL'
        results.append(result)
    else:
        print("数据量不足，无法进行全局拟合", file=sys.stderr)

    # 保存所有拟合结果到CSV
    if results:
        results_df = pd.DataFrame(results)
        results_file = os.path.join(output_dir, 'fitting_results.csv')
        results_df.to_csv(results_file, index=False)

    print(f"所有拟合已完成，结果已保存到 {output_dir}")


if __name__ == "__main__":
    main()
