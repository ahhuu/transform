#!/usr/bin/env python3
"""
GNSS 时间日期转换工具
========================
功能:
  1. 年月日 → DOY, MJD, JD, GPS周/周内日, BDS周/周内日, GST周/周内日, GLONASS日
     并自动显示对应的北京时间 (CST, UTC+8)
  2. UTC / CST / GPST / BDT / GST / GLONASST 多系统时间互转

时间系统说明:
  - UTC      : 协调世界时 (含跳秒)
  - CST      : 北京时间 (China Standard Time, UTC+8)
  - GPST     : GPS系统时, 历元 1980-01-06 00:00:00 UTC
  - BDT      : 北斗系统时, 历元 2006-01-01 00:00:00 UTC
  - GST      : 伽利略系统时, 历元 1999-08-22 00:00:00 UTC (≈GPST)
  - GLONASST : 格洛纳斯系统时, UTC + 3h (含跳秒修正)
"""

import tkinter as tk
from tkinter import ttk, messagebox
from datetime import datetime
import math

# =================================================================
# 1. 全局常量
# =================================================================

# 各GNSS系统历元的儒略日 (JD)
JD_GPS_EPOCH  = 2444244.5   # 1980-01-06 00:00:00 UTC
JD_BDS_EPOCH  = 2453736.5   # 2006-01-01 00:00:00 UTC
JD_GST_EPOCH  = 2451412.5   # 1999-08-22 00:00:00 UTC
JD_GLONASS_EPOCH = 2445157.5  # 实际上GLONASS没有独立的周计数，这里是UTC起点参考

# MJD 偏移量
MJD_OFFSET = 2400000.5

# 跳秒表: (UTC日期, 新增跳秒数)
# 说明: 从GPS历元(1980-01-06)以来的累积跳秒
# 格式: (年, 月, 日, 该日之后新增跳秒)
LEAP_SECONDS_TABLE = [
    (1980, 1,  1,  19),   # GPS历元时 TAI-UTC = 19s, GPS = TAI - 19s
    (1981, 7,  1,  20),
    (1982, 7,  1,  21),
    (1983, 7,  1,  22),
    (1985, 7,  1,  23),
    (1988, 1,  1,  24),
    (1990, 1,  1,  25),
    (1991, 1,  1,  26),
    (1992, 7,  1,  27),
    (1993, 7,  1,  28),
    (1994, 7,  1,  29),
    (1996, 1,  1,  30),
    (1997, 7,  1,  31),
    (1999, 1,  1,  32),
    (2006, 1,  1,  33),
    (2009, 1,  1,  34),
    (2012, 7,  1,  35),
    (2015, 7,  1,  36),
    (2017, 1,  1,  37),   # TAI-UTC = 37s
]

# 当前累积跳秒 (TAI - UTC)
# GPST = TAI - 19s => GPST - UTC = (TAI - 19) - UTC = (TAI - UTC) - 19
# 当前 TAI-UTC = 37, 所以 GPST-UTC = 37-19 = 18s

# 北京时间偏移 (CST, China Standard Time, UTC+8)
LOCAL_TZ_OFFSET = 8


# =================================================================
# 2. 核心计算函数
# =================================================================

def get_leap_seconds(year, month, day):
    """获取指定UTC日期时的 TAI - UTC 跳秒数"""
    tai_minus_utc = 10  # 1972年之前为10s
    for y, m, d, ls in LEAP_SECONDS_TABLE:
        if (year, month, day) >= (y, m, d):
            tai_minus_utc = ls
        else:
            break
    return tai_minus_utc


def get_gps_utc_offset(year, month, day):
    """
    计算指定UTC日期时 GPST - UTC 的差值(秒)。
    GPST = TAI - 19s, 所以 GPST - UTC = (TAI - UTC) - 19
    """
    tai_utc = get_leap_seconds(year, month, day)
    return tai_utc - 19


def get_bdt_utc_offset(year, month, day):
    """
    计算指定UTC日期时 BDT - UTC 的差值(秒)。
    BDT = TAI - 33s, 所以 BDT - UTC = (TAI - UTC) - 33
    """
    tai_utc = get_leap_seconds(year, month, day)
    return tai_utc - 33


def julian_day(year, month, day):
    """
    计算给定日期(UTC)的儒略日 JD (以世界时12h为日界)
    输入: year, month, day (均为整数或浮点数)
    返回: JD (float)
    """
    if month <= 2:
        year -= 1
        month += 12
    A = int(year / 100)
    B = 2 - A + int(A / 4)
    JD = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + B - 1524.5
    return JD


def jd_to_mjd(jd):
    """儒略日 → 简化儒略日"""
    return jd - MJD_OFFSET


def mjd_to_jd(mjd):
    """简化儒略日 → 儒略日"""
    return mjd + MJD_OFFSET


def jd_to_date(jd):
    """
    儒略日 → 年/月/日 (UTC)
    返回: (year, month, day, hour, minute, second)
    """
    jd = jd + 0.5
    Z = int(jd)
    F = jd - Z

    if Z < 2299161:
        A = Z
    else:
        alpha = int((Z - 1867216.25) / 36524.25)
        A = Z + 1 + alpha - int(alpha / 4)

    B = A + 1524
    C = int((B - 122.1) / 365.25)
    D = int(365.25 * C)
    E = int((B - D) / 30.6001)

    day = B - D - int(30.6001 * E) + F
    month = E - 1 if E < 14 else E - 13
    year = C - 4716 if month > 2 else C - 4715

    day_int = int(day)
    frac = day - day_int
    total_seconds = frac * 86400
    hour = int(total_seconds // 3600)
    minute = int((total_seconds % 3600) // 60)
    second = total_seconds % 60

    return year, month, day_int, hour, minute, second


def day_of_year(year, month, day):
    """年积日 DOY (1月1日为1)"""
    return datetime(year, month, day).timetuple().tm_yday


def _weeks_and_doy_from_jd(jd, epoch_jd):
    """
    通用函数: 根据儒略日和系统历元计算周数和周内日
    返回: (week_number, day_of_week)
    day_of_week: 0=星期日, 1=星期一, ..., 6=星期六
    """
    days_since_epoch = jd - epoch_jd
    week = int(math.floor(days_since_epoch / 7))
    dow = days_since_epoch - week * 7  # 0=星期日
    return week, dow


def compute_gps_time(jd):
    """计算GPS周和GPS周内日"""
    return _weeks_and_doy_from_jd(jd, JD_GPS_EPOCH)


def compute_bds_time(jd):
    """计算BDS周和BDS周内日"""
    return _weeks_and_doy_from_jd(jd, JD_BDS_EPOCH)


def compute_gst_time(jd):
    """计算GST周和GST周内日"""
    return _weeks_and_doy_from_jd(jd, JD_GST_EPOCH)


def compute_all_from_ymd(year, month, day, hour=0, minute=0, second=0.0):
    """
    根据年月日时分秒计算所有时间参数
    返回字典
    """
    # 将时分秒转为天的分数
    day_frac = day + hour / 24 + minute / 1440 + second / 86400
    jd = julian_day(year, month, day_frac)
    mjd = jd_to_mjd(jd)
    doy = day_of_year(year, month, day)

    gps_week, gps_dow = compute_gps_time(jd)
    bds_week, bds_dow = compute_bds_time(jd)
    gst_week, gst_dow = compute_gst_time(jd)

    # 时间系统间偏移
    leap_sec = get_leap_seconds(year, month, day)
    gps_utc_offset = get_gps_utc_offset(year, month, day)
    bdt_utc_offset = get_bdt_utc_offset(year, month, day)

    return {
        "jd": jd,
        "mjd": mjd,
        "doy": doy,
        "gps_week": int(gps_week),
        "gps_dow": int(gps_dow),
        "bds_week": int(bds_week),
        "bds_dow": int(bds_dow),
        "gst_week": int(gst_week),
        "gst_dow": int(gst_dow),
        "leap_seconds": leap_sec,
        "gps_utc_offset": gps_utc_offset,
        "bdt_utc_offset": bdt_utc_offset,
    }


def time_system_convert(source_sys, dest_sys, utc_jd):
    """
    不同时间系统之间的转换（基于UTC JD转换到目标系统）。
    
    参数:
        source_sys: 源时间系统名称 (仅供错误提示)
        dest_sys:   目标时间系统名称 ('UTC', 'CST', 'GPST', 'BDT', 'GST', 'GLONASST')
        utc_jd:    UTC 下的儒略日
    
    返回: 目标系统下的儒略日
    """
    if source_sys == dest_sys:
        return utc_jd

    # UTC → 目标系统
    if dest_sys == "UTC":
        return utc_jd
    elif dest_sys == "CST":
        return utc_jd + LOCAL_TZ_OFFSET / 24
    elif dest_sys in ("GPST", "GST"):
        y, m, d, h, mi, s = jd_to_date(utc_jd)
        offset = get_gps_utc_offset(int(y), int(m), int(d))
        return utc_jd + offset / 86400
    elif dest_sys == "BDT":
        y, m, d, h, mi, s = jd_to_date(utc_jd)
        offset = get_bdt_utc_offset(int(y), int(m), int(d))
        return utc_jd + offset / 86400
    elif dest_sys == "GLONASST":
        return utc_jd + 3 / 24
    else:
        raise ValueError(f"不支持的目标时间系统: {dest_sys}")


# =================================================================
# 3. 周内日名称映射
# =================================================================

DOW_NAMES = {
    0: "星期日 (Sunday)",
    1: "星期一 (Monday)",
    2: "星期二 (Tuesday)",
    3: "星期三 (Wednesday)",
    4: "星期四 (Thursday)",
    5: "星期五 (Friday)",
    6: "星期六 (Saturday)",
}

SYSTEM_NAMES = {
    "GPS":  "GPS系统时 (GPST)",
    "BDS":  "北斗系统时 (BDT)",
    "GST":  "伽利略系统时 (GST)",
    "GLONASS": "格洛纳斯系统时 (GLONASST)",
}

DOW_SHORT = {
    0: "Sun", 1: "Mon", 2: "Tue", 3: "Wed",
    4: "Thu", 5: "Fri", 6: "Sat",
}


# =================================================================
# 4. GUI 界面
# =================================================================

class TimeDateConverter(tk.Toplevel):
    """GNSS时间日期转换工具主窗口"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.title("GNSS 时间日期转换工具")
        self.geometry("780x880")
        self.resizable(True, True)

        # 如果不是独立窗口且有父窗口，设置模态
        if parent:
            self.transient(parent)
            self.grab_set()

        # 样式
        self.style = ttk.Style()
        self.style.configure("Title.TLabel", font=("微软雅黑", 12, "bold"))
        self.style.configure("Result.TLabel", font=("Consolas", 10))
        self.style.configure("Section.TLabelframe.Label", font=("微软雅黑", 10, "bold"))

        self._updating = False  # 双向同步防递归标志
        self._last_result1 = None  # 模块1结果缓存，用于联动模块2

        self._build_ui()

        # 初始化：填入当前时间
        self._fill_current_time()

        self.protocol("WM_DELETE_WINDOW", self._on_close)

    def _on_close(self):
        self.master.destroy()

    def _build_ui(self):
        """构建界面"""
        main_frame = ttk.Frame(self, padding="15")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # ==================== 模块1: 日期 → 多系统时间参数 ====================
        frame1 = ttk.LabelFrame(main_frame, text="模块1: 日期 → 多系统时间参数", 
                                style="Section.TLabelframe", padding="10")
        frame1.pack(fill=tk.X, pady=(0, 10))

        # 输入行
        input_frame = ttk.Frame(frame1)
        input_frame.pack(fill=tk.X, pady=5)

        ttk.Label(input_frame, text="UTC", font=("Consolas", 10, "bold")).grid(row=0, column=0, padx=(0, 4))
        ttk.Label(input_frame, text="年:").grid(row=0, column=1, padx=2)
        self.entry_year = ttk.Entry(input_frame, width=8)
        self.entry_year.grid(row=0, column=2, padx=2)
        self.entry_year.bind("<KeyRelease>", lambda e: self._update_time_display("utc"))

        ttk.Label(input_frame, text="月:").grid(row=0, column=3, padx=2)
        self.entry_month = ttk.Entry(input_frame, width=6)
        self.entry_month.grid(row=0, column=4, padx=2)
        self.entry_month.bind("<KeyRelease>", lambda e: self._update_time_display("utc"))

        ttk.Label(input_frame, text="日:").grid(row=0, column=5, padx=2)
        self.entry_day = ttk.Entry(input_frame, width=6)
        self.entry_day.grid(row=0, column=6, padx=2)
        self.entry_day.bind("<KeyRelease>", lambda e: self._update_time_display("utc"))

        ttk.Label(input_frame, text="时:").grid(row=0, column=7, padx=2)
        self.entry_hour = ttk.Entry(input_frame, width=6)
        self.entry_hour.insert(0, "0")
        self.entry_hour.grid(row=0, column=8, padx=2)
        self.entry_hour.bind("<KeyRelease>", lambda e: self._update_time_display("utc"))

        ttk.Label(input_frame, text="分:").grid(row=0, column=9, padx=2)
        self.entry_min = ttk.Entry(input_frame, width=6)
        self.entry_min.insert(0, "0")
        self.entry_min.grid(row=0, column=10, padx=2)
        self.entry_min.bind("<KeyRelease>", lambda e: self._update_time_display("utc"))

        ttk.Label(input_frame, text="秒:").grid(row=0, column=11, padx=2)
        self.entry_sec = ttk.Entry(input_frame, width=8)
        self.entry_sec.insert(0, "0")
        self.entry_sec.grid(row=0, column=12, padx=2)
        self.entry_sec.bind("<KeyRelease>", lambda e: self._update_time_display("utc"))

        # CST 输入行 (第二行)
        ttk.Label(input_frame, text="CST", font=("Consolas", 10, "bold")).grid(row=1, column=0, padx=(0, 4))
        ttk.Label(input_frame, text="年:").grid(row=1, column=1, padx=2)
        self.cst_year = ttk.Entry(input_frame, width=8)
        self.cst_year.grid(row=1, column=2, padx=2)
        self.cst_year.bind("<KeyRelease>", lambda e: self._update_time_display("cst"))

        ttk.Label(input_frame, text="月:").grid(row=1, column=3, padx=2)
        self.cst_month = ttk.Entry(input_frame, width=6)
        self.cst_month.grid(row=1, column=4, padx=2)
        self.cst_month.bind("<KeyRelease>", lambda e: self._update_time_display("cst"))

        ttk.Label(input_frame, text="日:").grid(row=1, column=5, padx=2)
        self.cst_day = ttk.Entry(input_frame, width=6)
        self.cst_day.grid(row=1, column=6, padx=2)
        self.cst_day.bind("<KeyRelease>", lambda e: self._update_time_display("cst"))

        ttk.Label(input_frame, text="时:").grid(row=1, column=7, padx=2)
        self.cst_hour = ttk.Entry(input_frame, width=6)
        self.cst_hour.insert(0, "0")
        self.cst_hour.grid(row=1, column=8, padx=2)
        self.cst_hour.bind("<KeyRelease>", lambda e: self._update_time_display("cst"))

        ttk.Label(input_frame, text="分:").grid(row=1, column=9, padx=2)
        self.cst_min = ttk.Entry(input_frame, width=6)
        self.cst_min.insert(0, "0")
        self.cst_min.grid(row=1, column=10, padx=2)
        self.cst_min.bind("<KeyRelease>", lambda e: self._update_time_display("cst"))

        ttk.Label(input_frame, text="秒:").grid(row=1, column=11, padx=2)
        self.cst_sec = ttk.Entry(input_frame, width=8)
        self.cst_sec.insert(0, "0")
        self.cst_sec.grid(row=1, column=12, padx=2)
        self.cst_sec.bind("<KeyRelease>", lambda e: self._update_time_display("cst"))

        # 按钮行
        btn_frame = ttk.Frame(frame1)
        btn_frame.pack(fill=tk.X, pady=5)
        ttk.Button(btn_frame, text="转换 →", command=self._convert_date).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="填入当前UTC时间", command=self._fill_current_time).pack(side=tk.LEFT, padx=5)

        # 结果显示
        self.result_text1 = tk.Text(frame1, height=14, width=85, font=("Consolas", 10),
                                    relief=tk.SUNKEN, borderwidth=1)
        self.result_text1.pack(fill=tk.X, pady=5)

        # 发送到模块2的按钮
        link_frame = ttk.Frame(frame1)
        link_frame.pack(fill=tk.X, pady=(0, 5))
        ttk.Label(link_frame, text="发送到模块2:", font=("微软雅黑", 9)).pack(side=tk.LEFT, padx=5)
        ttk.Button(link_frame, text="→ GPST", width=8,
                   command=lambda: self._send_to_module2("GPST")).pack(side=tk.LEFT, padx=2)
        ttk.Button(link_frame, text="→ BDT", width=8,
                   command=lambda: self._send_to_module2("BDT")).pack(side=tk.LEFT, padx=2)
        ttk.Button(link_frame, text="→ GST", width=8,
                   command=lambda: self._send_to_module2("GST")).pack(side=tk.LEFT, padx=2)

        # ==================== 模块2: 时间系统互转 ====================
        frame2 = ttk.LabelFrame(main_frame, text="模块2: 多系统时间互转",
                                style="Section.TLabelframe", padding="10")
        frame2.pack(fill=tk.X, pady=(0, 10))

        # 源时间系统
        sys_frame = ttk.Frame(frame2)
        sys_frame.pack(fill=tk.X, pady=5)

        ttk.Label(sys_frame, text="源系统:").grid(row=0, column=0, padx=5)
        self.src_sys = ttk.Combobox(sys_frame, values=["UTC", "CST", "GPST", "BDT", "GST", "GLONASST"],
                                    state="readonly", width=12)
        self.src_sys.current(0)
        self.src_sys.grid(row=0, column=1, padx=5)
        self.src_sys.bind("<<ComboboxSelected>>", self._on_src_sys_change)

        # 源系统时间输入 (根据选择显示不同的输入框)
        self.src_time_frame = ttk.Frame(sys_frame)
        self.src_time_frame.grid(row=0, column=2, columnspan=8, padx=5, sticky=tk.W)

        # 默认: 年/月/日/时/分/秒
        self.src_labels = []
        self.src_entries = []
        self._build_src_time_entries_ymd()

        ttk.Label(sys_frame, text="目标系统:").grid(row=1, column=0, padx=5, pady=5)
        self.dst_sys = ttk.Combobox(sys_frame, values=["UTC", "CST", "GPST", "BDT", "GST", "GLONASST"],
                                    state="readonly", width=12)
        self.dst_sys.current(1)  # 默认 GPST
        self.dst_sys.grid(row=1, column=1, padx=5, pady=5)

        # 转换按钮
        ttk.Button(sys_frame, text="转换 →", command=self._convert_system).grid(
            row=1, column=2, padx=10, pady=5)

        # 模块2结果显示
        self.result_text2 = tk.Text(frame2, height=10, width=85, font=("Consolas", 10),
                                    relief=tk.SUNKEN, borderwidth=1)
        self.result_text2.pack(fill=tk.X, pady=5)

        # ==================== 底部说明 ====================
        info_frame = ttk.LabelFrame(main_frame, text="时间系统说明", padding="5")
        info_frame.pack(fill=tk.X)

        info_text = (
            "• UTC      : 协调世界时 (含跳秒)\n"
            "• CST      : 北京时间, UTC+8 (China Standard Time)\n"
            "• GPST     : GPS系统时, 历元 1980-01-06 00:00:00 UTC, 无跳秒\n"
            "• BDT      : 北斗系统时, 历元 2006-01-01 00:00:00 UTC, 无跳秒\n"
            "• GST      : 伽利略系统时, 历元 1999-08-22 00:00:00 UTC, 与GPST同步\n"
            "• GLONASST : 格洛纳斯系统时, 比UTC快3小时 (UTC+3h)\n"
            "📌 跳秒数据更新至 TAI-UTC = 37s (2017年1月至今), GPST-UTC = 18s"
        )
        ttk.Label(info_frame, text=info_text, font=("微软雅黑", 9)).pack(anchor=tk.W, pady=3)

    def _build_src_time_entries_ymd(self):
        """构建源系统时间输入的 年/月/日/时/分/秒 输入框"""
        for w in self.src_time_frame.winfo_children():
            w.destroy()
        self.src_labels = []
        self.src_entries = []

        fields = [("年:", 8), ("月:", 6), ("日:", 6), ("时:", 6), ("分:", 6), ("秒:", 8)]
        for i, (label, width) in enumerate(fields):
            lb = ttk.Label(self.src_time_frame, text=label)
            lb.grid(row=0, column=i*2, padx=1)
            entry = ttk.Entry(self.src_time_frame, width=width)
            entry.grid(row=0, column=i*2+1, padx=1)
            self.src_labels.append(lb)
            self.src_entries.append(entry)

        # 填入默认值
        self.src_entries[0].insert(0, "2026")
        self.src_entries[1].insert(0, "1")
        self.src_entries[2].insert(0, "1")
        self.src_entries[3].insert(0, "0")
        self.src_entries[4].insert(0, "0")
        self.src_entries[5].insert(0, "0")

    def _build_src_time_entries_week(self):
        """构建源系统时间输入的 周/周内日/时/分/秒 输入框"""
        for w in self.src_time_frame.winfo_children():
            w.destroy()
        self.src_labels = []
        self.src_entries = []

        fields = [("周:", 10), ("周内日:", 8), ("时:", 6), ("分:", 6), ("秒:", 8)]
        for i, (label, width) in enumerate(fields):
            lb = ttk.Label(self.src_time_frame, text=label)
            lb.grid(row=0, column=i*2, padx=1)
            entry = ttk.Entry(self.src_time_frame, width=width)
            entry.grid(row=0, column=i*2+1, padx=1)
            self.src_labels.append(lb)
            self.src_entries.append(entry)

        self.src_entries[0].insert(0, "2400")
        self.src_entries[1].insert(0, "0")
        self.src_entries[2].insert(0, "0")
        self.src_entries[3].insert(0, "0")
        self.src_entries[4].insert(0, "0")

    def _on_src_sys_change(self, event=None):
        """源系统改变时切换输入方式"""
        sys_name = self.src_sys.get()
        if sys_name in ("UTC", "CST"):
            self._build_src_time_entries_ymd()
        else:
            self._build_src_time_entries_week()

    def _fill_current_time(self):
        """填入当前UTC时间，并同步CST"""
        now = datetime.utcnow()
        self.entry_year.delete(0, tk.END)
        self.entry_year.insert(0, str(now.year))
        self.entry_month.delete(0, tk.END)
        self.entry_month.insert(0, str(now.month))
        self.entry_day.delete(0, tk.END)
        self.entry_day.insert(0, str(now.day))
        self.entry_hour.delete(0, tk.END)
        self.entry_hour.insert(0, str(now.hour))
        self.entry_min.delete(0, tk.END)
        self.entry_min.insert(0, str(now.minute))
        self.entry_sec.delete(0, tk.END)
        self.entry_sec.insert(0, str(now.second))
        self._update_time_display("utc")

    def _update_time_display(self, source="utc"):
        """双向同步 UTC ↔ CST 输入
        source="utc" → 修改UTC，更新CST行
        source="cst" → 修改CST，更新UTC行
        """
        if self._updating:
            return
        self._updating = True
        try:
            if source == "utc":
                year   = int(self.entry_year.get().strip() or "0")
                month  = int(self.entry_month.get().strip() or "0")
                day    = int(self.entry_day.get().strip() or "0")
                hour   = int(self.entry_hour.get().strip() or "0")
                minute = int(self.entry_min.get().strip() or "0")
                second = float(self.entry_sec.get().strip() or "0")
                if year < 1 or month < 1 or month > 12 or day < 1 or day > 31:
                    raise ValueError
                day_frac = day + hour / 24 + minute / 1440 + second / 86400
                jd = julian_day(year, month, day_frac)
                local_jd = jd + LOCAL_TZ_OFFSET / 24
                local = jd_to_date(local_jd)

                self.cst_year.delete(0, tk.END); self.cst_year.insert(0, str(local[0]))
                self.cst_month.delete(0, tk.END); self.cst_month.insert(0, str(local[1]))
                self.cst_day.delete(0, tk.END); self.cst_day.insert(0, str(local[2]))
                self.cst_hour.delete(0, tk.END); self.cst_hour.insert(0, str(local[3]))
                self.cst_min.delete(0, tk.END); self.cst_min.insert(0, str(local[4]))
                self.cst_sec.delete(0, tk.END); self.cst_sec.insert(0, f"{local[5]:.1f}")
            else:
                year   = int(self.cst_year.get().strip() or "0")
                month  = int(self.cst_month.get().strip() or "0")
                day    = int(self.cst_day.get().strip() or "0")
                hour   = int(self.cst_hour.get().strip() or "0")
                minute = int(self.cst_min.get().strip() or "0")
                second = float(self.cst_sec.get().strip() or "0")
                if year < 1 or month < 1 or month > 12 or day < 1 or day > 31:
                    raise ValueError
                day_frac = day + hour / 24 + minute / 1440 + second / 86400
                cst_jd = julian_day(year, month, day_frac)
                utc_jd = cst_jd - LOCAL_TZ_OFFSET / 24
                utc = jd_to_date(utc_jd)

                self.entry_year.delete(0, tk.END); self.entry_year.insert(0, str(utc[0]))
                self.entry_month.delete(0, tk.END); self.entry_month.insert(0, str(utc[1]))
                self.entry_day.delete(0, tk.END); self.entry_day.insert(0, str(utc[2]))
                self.entry_hour.delete(0, tk.END); self.entry_hour.insert(0, str(utc[3]))
                self.entry_min.delete(0, tk.END); self.entry_min.insert(0, str(utc[4]))
                self.entry_sec.delete(0, tk.END); self.entry_sec.insert(0, f"{utc[5]:.1f}")
        except (ValueError, ZeroDivisionError):
            pass
        finally:
            self._updating = False

    def _format_result1(self, data):
        """格式化模块1的转换结果"""
        dow_gps = DOW_NAMES.get(data["gps_dow"], f"Day {data['gps_dow']}")
        dow_bds = DOW_NAMES.get(data["bds_dow"], f"Day {data['bds_dow']}")
        dow_gst = DOW_NAMES.get(data["gst_dow"], f"Day {data['gst_dow']}")

        lines = []
        lines.append(f"{'='*60}")
        lines.append(f"  日期时间转换结果")
        lines.append(f"{'='*60}")
        lines.append(f"")
        lines.append(f"  儒略日        JD  = {data['jd']:.5f}")
        lines.append(f"  简化儒略日    MJD = {data['mjd']:.5f}")
        lines.append(f"  年积日        DOY = {data['doy']}")
        lines.append(f"")

        lines.append(f"  {'─'*50}")
        lines.append(f"  【GPS系统时 - GPST】  (历元: 1980-01-06 UTC)")
        lines.append(f"  GPS 周        = {data['gps_week']}")
        lines.append(f"  GPS 周内日    = {data['gps_dow']}  ({dow_gps})")
        lines.append(f"")
        lines.append(f"  【北斗系统时 - BDT】  (历元: 2006-01-01 UTC)")
        lines.append(f"  BDS 周        = {data['bds_week']}")
        lines.append(f"  BDS 周内日    = {data['bds_dow']}  ({dow_bds})")
        lines.append(f"")
        lines.append(f"  【伽利略系统时 - GST】 (历元: 1999-08-22 UTC)")
        lines.append(f"  GST 周        = {data['gst_week']}")
        lines.append(f"  GST 周内日    = {data['gst_dow']}  ({dow_gst})")
        lines.append(f"")
        lines.append(f"  {'─'*50}")
        lines.append(f"  TAI - UTC = {data['leap_seconds']}s")
        lines.append(f"  GPST - UTC = {data['gps_utc_offset']}s")
        lines.append(f"  BDT  - UTC = {data['bdt_utc_offset']}s")
        lines.append(f"  GLONASST - UTC = +3h")
        lines.append(f"  GST - GPST ≈ 0s (同步精度<50ns)")
        lines.append(f"{'='*60}")

        return "\n".join(lines)

    def _convert_date(self):
        """模块1: 日期转换"""
        try:
            year = int(self.entry_year.get().strip())
            month = int(self.entry_month.get().strip())
            day = int(self.entry_day.get().strip())
            hour = int(self.entry_hour.get().strip() or "0")
            minute = int(self.entry_min.get().strip() or "0")
            second = float(self.entry_sec.get().strip() or "0")

            # 验证
            if month < 1 or month > 12:
                raise ValueError("月份必须在1-12之间")
            if day < 1 or day > 31:
                raise ValueError("日期必须在1-31之间")

            data = compute_all_from_ymd(year, month, day, hour, minute, second)
            self._last_result1 = data
            result = self._format_result1(data)

            self.result_text1.delete(1.0, tk.END)
            self.result_text1.insert(1.0, result)

        except ValueError as e:
            messagebox.showerror("输入错误", str(e))
        except Exception as e:
            messagebox.showerror("转换错误", f"日期转换失败: {e}")

    def _send_to_module2(self, target_sys):
        """将模块1的转换结果发送到模块2进行进一步转换"""
        if self._last_result1 is None:
            messagebox.showinfo("提示", "请先在模块1中点击「转换 →」")
            return

        data = self._last_result1

        # 设置模块2的源系统
        self.src_sys.set(target_sys)
        self._on_src_sys_change()

        # 填入周/周内日/时/分/秒
        sys_key = {"GPST": "gps", "BDT": "bds", "GST": "gst"}
        key = sys_key.get(target_sys, "gps")
        self.src_entries[0].delete(0, tk.END)
        self.src_entries[0].insert(0, str(data[f"{key}_week"]))
        self.src_entries[1].delete(0, tk.END)
        self.src_entries[1].insert(0, str(data[f"{key}_dow"]))
        # 时/分/秒清零
        for i in range(2, 5):
            self.src_entries[i].delete(0, tk.END)
            self.src_entries[i].insert(0, "0")

        # 自动执行转换 (目标系统设为 UTC)
        self.dst_sys.current(0)
        self._convert_system()

    def _convert_system(self):
        """模块2: 时间系统互转"""
        try:
            src = self.src_sys.get()
            dst = self.dst_sys.get()

            if src == dst:
                messagebox.showinfo("提示", "源系统与目标系统相同，无需转换")
                return

            if src in ("UTC", "CST"):
                # 读取年月日时分秒
                vals = [e.get().strip() for e in self.src_entries]
                year, month, day = int(vals[0]), int(vals[1]), int(vals[2])
                hour, minute = int(vals[3] or "0"), int(vals[4] or "0")
                second = float(vals[5] or "0")
                day_frac = day + hour / 24 + minute / 1440 + second / 86400
                src_jd = julian_day(year, month, day_frac)
            else:
                # 读取周/周内日/时/分/秒
                vals = [e.get().strip() for e in self.src_entries]
                week = int(vals[0])
                dow = int(vals[1])
                hour = int(vals[2] or "0")
                minute = int(vals[3] or "0")
                second = float(vals[4] or "0")

                # 确定历元JD
                epoch_map = {"GPST": JD_GPS_EPOCH, "BDT": JD_BDS_EPOCH, "GST": JD_GST_EPOCH}
                epoch_jd = epoch_map.get(src, JD_GPS_EPOCH)

                # 计算源系统JD: 历元JD + 周数*7 + 周内日 + 日内秒/86400
                src_jd = epoch_jd + week * 7 + dow + (hour / 24 + minute / 1440 + second / 86400)

            # 计算UTC JD (用于北京时间显示)
            if src == "UTC":
                utc_jd = src_jd
            elif src == "CST":
                utc_jd = src_jd - LOCAL_TZ_OFFSET / 24
            elif src in ("GPST", "GST"):
                y, m, d, *_ = jd_to_date(src_jd)
                offset = get_gps_utc_offset(int(y), int(m), int(d))
                utc_jd = src_jd - offset / 86400
            elif src == "BDT":
                y, m, d, *_ = jd_to_date(src_jd)
                offset = get_bdt_utc_offset(int(y), int(m), int(d))
                utc_jd = src_jd - offset / 86400
            elif src == "GLONASST":
                utc_jd = src_jd - 3 / 24
            else:
                utc_jd = src_jd

            # 执行转换 (从UTC到目标系统)
            result_jd = time_system_convert("UTC", dst, utc_jd)

            # 格式化结果
            result_date = jd_to_date(result_jd)

            # 计算对应的北京时间 (CST = UTC+8)
            local_jd = utc_jd + LOCAL_TZ_OFFSET / 24
            local_date = jd_to_date(local_jd)

            lines = []
            lines.append(f"{'='*55}")
            lines.append(f"  {src} → {dst}  转换结果")
            lines.append(f"{'='*55}")
            lines.append(f"")

            if dst in ("UTC", "CST"):
                # 显示年月日时分秒 (CST是民用时，同UTC格式)
                sys_label = "CST (北京时间)" if dst == "CST" else dst
                lines.append(f"  {sys_label} 日期时间:")
                lines.append(f"  {result_date[0]:04d}-{result_date[1]:02d}-{result_date[2]:02d} "
                             f"{result_date[3]:02d}:{result_date[4]:02d}:{result_date[5]:.3f}")
                lines.append(f"  JD  = {result_jd:.5f}")
                lines.append(f"  MJD = {jd_to_mjd(result_jd):.5f}")
            else:
                # 显示周/周内日
                lines.append(f"  {dst} 时间参数:")
                lines.append(f"  JD = {result_jd:.5f}")
                if dst == "GPST":
                    w, d = compute_gps_time(result_jd)
                elif dst == "BDT":
                    w, d = compute_bds_time(result_jd)
                elif dst == "GST":
                    w, d = compute_gst_time(result_jd)
                else:
                    w, d = 0, 0
                lines.append(f"  周(Week) = {int(w)}")
                dow_name = DOW_NAMES.get(int(d), f"Day {int(d)}")
                lines.append(f"  周内日(DOW) = {int(d)} ({dow_name})")
                # 日内秒
                total_sec = result_date[3] * 3600 + result_date[4] * 60 + result_date[5]
                lines.append(f"  日内秒(SOD) = {total_sec:.3f}s")

            # 北京时间 (CST) 显示 (目标不是 CST 时才显示)
            if dst != "CST":
                lines.append(f"")
                lines.append(f"  {'─'*35}")
                lines.append(f"  对应北京时间(CST):")
                lines.append(f"  {local_date[0]:04d}-{local_date[1]:02d}-{local_date[2]:02d} "
                             f"{local_date[3]:02d}:{local_date[4]:02d}:{local_date[5]:.1f}  (UTC+8)")

            lines.append(f"{'='*55}")

            self.result_text2.delete(1.0, tk.END)
            self.result_text2.insert(1.0, "\n".join(lines))

        except ValueError as e:
            messagebox.showerror("输入错误", str(e))
        except Exception as e:
            messagebox.showerror("转换错误", f"系统转换失败: {e}")


# =================================================================
# 5. 独立运行入口
# =================================================================

def main():
    root = tk.Tk()
    root.withdraw()  # 隐藏主窗口
    app = TimeDateConverter()
    app.mainloop()


if __name__ == "__main__":
    main()
