import math


def calculate_phone_coordinates(lat_rtk, lon_rtk, alt_rtk,
                                ΔE_cm, ΔN_cm, ΔU_cm=0):
    """
    根据RTK坐标和手机相对偏移量计算手机的真值坐标。

    参数:
        lat_rtk, lon_rtk, alt_rtk: RTK天线的WGS84坐标（单位：度, 度, 米）
        ΔE_cm: 手机东方向偏移（厘米，东正西负）
        ΔN_cm: 手机北方向偏移（厘米，北正南负）
        ΔU_cm: 手机高程偏移（厘米，上正下负，默认0）

    返回:
        (lat_phone, lon_phone, alt_phone): 手机的WGS84坐标
    """
    # 厘米转米
    delta_east = ΔE_cm / 100
    delta_north = ΔN_cm / 100
    delta_up = ΔU_cm / 100

    # 计算纬度偏移（1米 ≈ 1/111320度）
    lat_phone = lat_rtk + (delta_north / 111320)

    # 计算经度偏移（考虑纬度圈收缩）
    lon_phone = lon_rtk + (delta_east / (111320 * math.cos(math.radians(lat_rtk))))

    # 高程直接相加
    alt_phone = alt_rtk + delta_up - 0.05
    return lat_phone, lon_phone, alt_phone


# 示例：手机在RTK西南8cm，高度一致
lat_phone, lon_phone, alt_phone = calculate_phone_coordinates(
    lat_rtk=30.702901527777, lon_rtk=104.05424769444, alt_rtk=458.684,
    ΔE_cm=-15 * math.cos(math.radians(45)),  # 东方向分量（西）
    ΔN_cm=-0 * math.sin(math.radians(45))  # 北方向分量（南）
)
print(f"Phone Coordinates: {lat_phone:.10f}, {lon_phone:.10f}, {alt_phone:.2f}")