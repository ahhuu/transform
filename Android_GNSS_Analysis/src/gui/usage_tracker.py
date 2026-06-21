"""
工具使用频率追踪模块
=====================
记录每个工具的启动次数，用于主界面动态显示最常用工具。
数据存储在工作空间根目录下的 .tool_usage.json 文件中。
"""

import json
import os
import threading


class UsageTracker:
    """工具使用频率追踪器"""

    def __init__(self, filepath=None):
        """
        参数:
            filepath: 存储文件路径，默认在项目根目录下 .tool_usage.json
        """
        if filepath is None:
            # 默认放在工程根目录
            gui_dir = os.path.dirname(os.path.abspath(__file__))
            project_root = os.path.abspath(os.path.join(gui_dir, "..", ".."))
            filepath = os.path.join(project_root, ".tool_usage.json")
        self._filepath = filepath
        self._lock = threading.Lock()
        self._data = self._load()

    def _load(self):
        """从文件加载使用数据"""
        try:
            if os.path.exists(self._filepath):
                with open(self._filepath, "r", encoding="utf-8") as f:
                    return json.load(f)
        except (json.JSONDecodeError, OSError):
            pass
        return {}

    def _save(self):
        """保存使用数据到文件"""
        try:
            with open(self._filepath, "w", encoding="utf-8") as f:
                json.dump(self._data, f, ensure_ascii=False, indent=2)
        except OSError:
            pass  # 保存失败静默处理

    def record_usage(self, label):
        """记录一次工具使用

        参数:
            label: 工具显示名称（与菜单中的标签一致）
        """
        with self._lock:
            self._data[label] = self._data.get(label, 0) + 1
            self._save()

    def get_top(self, n=4):
        """获取使用次数最多的前 n 个工具

        返回:
            list[tuple[label, count]] 按使用次数降序排列
        """
        with self._lock:
            sorted_items = sorted(
                self._data.items(), key=lambda x: x[1], reverse=True
            )
            return sorted_items[:n]

    def get_count(self, label):
        """获取指定工具的使用次数"""
        with self._lock:
            return self._data.get(label, 0)

    def get_all_sorted(self):
        """获取所有工具的使用数据（按次数降序）"""
        with self._lock:
            return sorted(
                self._data.items(), key=lambda x: x[1], reverse=True
            )

    def reset(self):
        """重置所有统计数据"""
        with self._lock:
            self._data.clear()
            self._save()

    @property
    def filepath(self):
        return self._filepath
