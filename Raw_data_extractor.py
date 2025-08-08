import os
from tkinter import Tk
from tkinter.filedialog import askopenfilename


def extract_and_save_data():
    # 初始化 Tkinter 根窗口并隐藏它
    root = Tk()
    root.withdraw()

    # 弹出文件选择对话框，仅允许选择 txt 文件
    input_file_path = askopenfilename(title="选择输入文件", filetypes=[("Text files", "*.txt")])

    if not input_file_path:
        print("未选择文件。")
        return

    version = ""
    data_lines = []

    # 读取文件
    with open(input_file_path, 'r', encoding='utf-8') as file:
        for line in file:
            if line.startswith('#Version'):
                version = line.strip()
            elif not line.startswith('#'):
                data_lines.append(line.strip())

    # 生成新文件名
    base_name, extension = os.path.splitext(os.path.basename(input_file_path))
    directory = os.path.dirname(input_file_path)

    if "gnss_log" in base_name:
        new_base_name = base_name.replace("gnss_log", "Raw")
    else:
        new_base_name = "Raw" + base_name

    new_file_name = os.path.join(directory, new_base_name + extension)

    # 写入新文件
    with open(new_file_name, 'w', encoding='utf-8') as new_file:
        if version:
            new_file.write(version + '\n')
        for line in data_lines:
            new_file.write(line + '\n')

    print(f"数据已成功提取并保存到 {new_file_name}")


if __name__ == "__main__":
    extract_and_save_data()
