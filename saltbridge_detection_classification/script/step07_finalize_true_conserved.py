#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import os
import argparse
from glob import glob
from collections import defaultdict

# 组头行，例如：3RGB-I280_414 与 3CHX-I295_426 : pi1-pi2:...
HEADER_RE = re.compile(r'^(\S+)\s+与\s+(\S+)\s*:\s*(.*)$')

# 盐桥明细行，例如：
# 3RGB-I280_414 pi_chain:I pi_resn:GLU pi_resi:316 cation_chain:I cation_resn:ARG cation_resi:323 distance:...
BRIDGE_RE = re.compile(
    r'^(?P<struct>\S+)\s+'
    r'pi_chain:(?P<pi_chain>\S+)\s+pi_resn:(?P<pi_resn>\S+)\s+pi_resi:(?P<pi_resi>\S+)\s+'
    r'cation_chain:(?P<cat_chain>\S+)\s+cation_resn:(?P<cat_resn>\S+)\s+cation_resi:(?P<cat_resi>\S+)'
)

def normalize_spaces(s: str) -> str:
    """压缩多余空白，便于文本一致性判断。"""
    return ' '.join(s.strip().split())

def bridge_key_line(line: str) -> str:
    """以整行（压缩空白后）作为盐桥的唯一键。最严格、和你示例一致。"""
    return normalize_spaces(line)

def bridge_key_residue(line: str) -> str:
    """
    以残基身份作为键（忽略结构ID、距离、坐标等）。
    若解析失败，回退为整行键。
    """
    m = BRIDGE_RE.match(normalize_spaces(line))
    if not m:
        return bridge_key_line(line)
    d = m.groupdict()
    # 只用链+残基名+编号来定义“同一盐桥”
    return f"{d['pi_chain']}:{d['pi_resn']}{d['pi_resi']}--{d['cat_chain']}:{d['cat_resn']}{d['cat_resi']}"

def parse_groups(lines):
    """
    将原始文本解析为“组”的列表。
    每组 = {'header': header_line, 'bridge_lines': [line1, line2, ...]}
    """
    groups = []
    current = None
    for raw in lines:
        line = raw.rstrip('\n')
        if not line.strip():
            continue
        if HEADER_RE.match(line):
            if current:
                groups.append(current)
            current = {'header': line, 'bridge_lines': []}
        else:
            if current is None:
                # 文件可能从明细开始，给它建一个无头组
                current = {'header': '(NO_HEADER)', 'bridge_lines': []}
            current['bridge_lines'].append(line)
    if current:
        groups.append(current)
    return groups

def process_one_file(in_path: str, out_path: str, k: int, key_mode: str) -> dict:
    """
    处理单个txt：
    - 找到出现次数 ≥ k 的盐桥
    - 输出包含这些盐桥的所有组（去重、按原顺序）到 out_path
    返回统计信息
    """
    with open(in_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    groups = parse_groups(lines)
    if not groups:
        # 空文件或无有效组
        open(out_path, 'w', encoding='utf-8').close()
        return {'groups': 0, 'selected_groups': 0, 'bridges': 0, 'selected_bridges': 0}

    # 选择键函数
    key_fn = bridge_key_line if key_mode == 'line' else bridge_key_residue

    # 统计：每条“盐桥键” -> 出现在哪些组（用组索引）
    bridge_to_groups = defaultdict(set)
    bridge_count = 0
    for gi, grp in enumerate(groups):
        for bl in grp['bridge_lines']:
            key = key_fn(bl)
            bridge_to_groups[key].add(gi)
            bridge_count += 1

    # 筛选：出现次数 ≥ k 的盐桥
    qualified_keys = {key for key, gset in bridge_to_groups.items() if len(gset) >= k}

    # 收集需要输出的组（出现过至少一条“达阈值盐桥”）
    selected_group_ids = []
    seen = set()
    for gi, grp in enumerate(groups):
        # 如果本组中存在任一达阈值盐桥，则选入
        has_qualified = any(key_fn(bl) in qualified_keys for bl in grp['bridge_lines'])
        if has_qualified and gi not in seen:
            seen.add(gi)
            selected_group_ids.append(gi)

    # 写出
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as wf:
        for gi in selected_group_ids:
            grp = groups[gi]
            wf.write(grp['header'] + "\n")
            for bl in grp['bridge_lines']:
                wf.write(bl + "\n")
            wf.write("\n")

    return {
        'groups': len(groups),
        'selected_groups': len(selected_group_ids),
        'bridges': bridge_count,
        'selected_bridges': len(qualified_keys)
    }

def main():
    parser = argparse.ArgumentParser(
        description="批量：对目录中每个txt，保留“包含出现次数≥k的盐桥”的所有组，并写入新目录下同名txt。"
    )
    parser.add_argument("indir", help="输入目录（包含多个 .txt）")
    parser.add_argument("-o", "--outdir", required=True, help="输出目录（将写入同名 .txt）")
    parser.add_argument("-k", "--min_occurrence", type=int, default=2, help="出现次数阈值（≥k），默认 2")
    parser.add_argument("--key-mode", choices=["line", "residue"], default="line",
                        help="判定同一盐桥的方式：line=整行一致；residue=按残基身份（链+残基名+编号）。默认 line")
    args = parser.parse_args()

    in_files = sorted(glob(os.path.join(args.indir, "*.txt")))
    if not in_files:
        raise SystemExit("输入目录下未找到 .txt 文件。")

    os.makedirs(args.outdir, exist_ok=True)

    total = 0
    kept_any = 0
    for in_path in in_files:
        base = os.path.basename(in_path)
        out_path = os.path.join(args.outdir, base)

        stats = process_one_file(in_path, out_path, k=args.min_occurrence, key_mode=args.key_mode)
        total += 1
        if stats['selected_groups'] > 0:
            kept_any += 1

        print(f"[{base}] groups={stats['groups']}, selected_groups={stats['selected_groups']}, "
              f"unique_bridges={stats['selected_bridges']} (k={args.min_occurrence}, key={args.key_mode})")

    print(f"\n完成：共处理 {total} 个文件；有输出内容的文件 {kept_any} 个。输出目录：{args.outdir}")

if __name__ == "__main__":
    main()

