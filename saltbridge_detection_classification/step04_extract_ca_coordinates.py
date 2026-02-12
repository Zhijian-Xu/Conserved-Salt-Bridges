#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
盐桥保守性分析脚本 - 简化版
功能：
1. 读取盐桥信息，叠合蛋白，计算CA原子距离
2. 处理多构象残基（自动尝试 alt A/B/C）
3. 根据距离分类保存结果（cla/nocla/other/error）
4. 可选生成PSE文件用于可视化
"""

import os
import sys
import re
import argparse
from pymol import cmd
import traceback

def load_reference_list(ref_file):
    """读取参考蛋白列表"""
    ref_dict = {}
    try:
        with open(ref_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    parts = line.split()
                    if len(parts) >= 2:
                        family_id = parts[0]
                        ref_protein = parts[1]
                        ref_dict[family_id] = ref_protein
        print(f"成功加载 {len(ref_dict)} 个参考蛋白")
        return ref_dict
    except Exception as e:
        print(f"读取参考列表失败: {e}")
        return {}

def parse_saltbridge_line(line):
    """解析盐桥信息行"""
    try:
        pi_chain_match = re.search(r'pi_chain:(\w+)', line)
        pi_resn_match = re.search(r'pi_resn:(\w+)', line)
        pi_resi_match = re.search(r'pi_resi:(-?\d+[A-Za-z]?)', line)
        
        cation_chain_match = re.search(r'cation_chain:(\w+)', line)
        cation_resn_match = re.search(r'cation_resn:(\w+)', line)
        cation_resi_match = re.search(r'cation_resi:(-?\d+[A-Za-z]?)', line)
        
        if all([pi_chain_match, pi_resn_match, pi_resi_match,
                cation_chain_match, cation_resn_match, cation_resi_match]):
            return {
                'pi_chain': pi_chain_match.group(1),
                'pi_resn': pi_resn_match.group(1),
                'pi_resi': pi_resi_match.group(1),
                'cation_chain': cation_chain_match.group(1),
                'cation_resn': cation_resn_match.group(1),
                'cation_resi': cation_resi_match.group(1)
            }
        else:
            print(f"  解析失败的行: {line}")
            print(f"    pi_chain: {'✓' if pi_chain_match else '✗'}")
            print(f"    pi_resn: {'✓' if pi_resn_match else '✗'}")
            print(f"    pi_resi: {'✓' if pi_resi_match else '✗'}")
            print(f"    cation_chain: {'✓' if cation_chain_match else '✗'}")
            print(f"    cation_resn: {'✓' if cation_resn_match else '✗'}")
            print(f"    cation_resi: {'✓' if cation_resi_match else '✗'}")
    except Exception as e:
        print(f"解析盐桥信息异常: {e}")
        print(f"  问题行: {line}")
    return None

def try_select_with_alt(protein_name, chain, resi, resn):
    """
    尝试选择残基CA原子，如果有多个构象则尝试alt A/B/C
    返回：(选择字符串, 原子数)
    """
    # 基础选择（不带alt）
    base_sel = f"{protein_name} and chain {chain} and resi {resi} and resn {resn} and name CA"
    count = cmd.count_atoms(base_sel)
    
    if count == 1:
        # 只有1个原子，直接返回
        return (base_sel, count)
    elif count == 0:
        # 找不到原子
        return (base_sel, 0)
    else:
        # 多个原子，尝试alt
        print(f"      残基 {chain}:{resi}:{resn} 选中 {count} 个原子，尝试 alt 过滤...")
        
        # 尝试常见的alt标识：A, B, C, D, 空
        for alt_id in ['A', 'B', 'C', 'D', '']:
            if alt_id:
                alt_sel = f"{protein_name} and chain {chain} and resi {resi} and resn {resn} and name CA and (alt ''+{alt_id})"
            else:
                alt_sel = f"{protein_name} and chain {chain} and resi {resi} and resn {resn} and name CA and alt ''"
            
            alt_count = cmd.count_atoms(alt_sel)
            
            if alt_count == 1:
                print(f"        ✓ 使用 alt '{alt_id}' 成功选中 1 个原子")
                return (alt_sel, 1)
            elif alt_count > 0:
                print(f"        - alt '{alt_id}' 选中 {alt_count} 个原子")
        
        # 所有alt都试过了，还是选不到单个原子
        print(f"        ✗ 尝试所有 alt 标识后仍无法选中单个原子")
        return (base_sel, count)

def get_protein_path(protein_name, family_id, base_dir):
    """构建蛋白质文件路径"""
    protein_file = f"{protein_name}.cif"
    protein_path = os.path.join(base_dir, family_id, protein_file)
    
    if not os.path.exists(protein_path):
        protein_path = os.path.join(base_dir, family_id, f"{protein_name}.pdb")
    
    return protein_path

def calculate_ca_distances(protein1_name, sb1_info, protein2_name, sb2_info):
    """
    计算两个蛋白质盐桥残基CA原子之间的距离
    自动处理多构象残基
    """
    try:
        print(f"    选择残基...")
        
        # 尝试选择蛋白1的残基
        p1_pi_sel, p1_pi_count = try_select_with_alt(
            protein1_name, sb1_info['pi_chain'], 
            sb1_info['pi_resi'], sb1_info['pi_resn']
        )
        
        p1_cation_sel, p1_cation_count = try_select_with_alt(
            protein1_name, sb1_info['cation_chain'],
            sb1_info['cation_resi'], sb1_info['cation_resn']
        )
        
        # 尝试选择蛋白2的残基
        p2_pi_sel, p2_pi_count = try_select_with_alt(
            protein2_name, sb2_info['pi_chain'],
            sb2_info['pi_resi'], sb2_info['pi_resn']
        )
        
        p2_cation_sel, p2_cation_count = try_select_with_alt(
            protein2_name, sb2_info['cation_chain'],
            sb2_info['cation_resi'], sb2_info['cation_resn']
        )
        
        # 打印PyMOL选择命令
        print(f"    PyMOL选择命令:")
        print(f"      p1_pi:      {p1_pi_sel}")
        print(f"      p1_cation:  {p1_cation_sel}")
        print(f"      p2_pi:      {p2_pi_sel}")
        print(f"      p2_cation:  {p2_cation_sel}")
        
        print(f"    选中原子数:")
        print(f"      p1_pi:      {p1_pi_count} 个原子")
        print(f"      p1_cation:  {p1_cation_count} 个原子")
        print(f"      p2_pi:      {p2_pi_count} 个原子")
        print(f"      p2_cation:  {p2_cation_count} 个原子")
        
        # 检查是否选中0个原子
        if p1_pi_count == 0:
            print(f"    ✗ 错误：找不到 {protein1_name} 的 pi 残基 CA 原子")
            return None
        if p1_cation_count == 0:
            print(f"    ✗ 错误：找不到 {protein1_name} 的 cation 残基 CA 原子")
            return None
        if p2_pi_count == 0:
            print(f"    ✗ 错误：找不到 {protein2_name} 的 pi 残基 CA 原子")
            return None
        if p2_cation_count == 0:
            print(f"    ✗ 错误：找不到 {protein2_name} 的 cation 残基 CA 原子")
            return None
        
        # 检查是否仍选中多个原子
        if any(c > 1 for c in [p1_pi_count, p1_cation_count, p2_pi_count, p2_cation_count]):
            print(f"    ✗ 错误：仍选中多个原子")
            return {'error': 'multiple_atoms'}
        
        # 计算距离
        print(f"    计算距离...")
        dist_pi1_pi2 = cmd.distance("tmp_dist", p1_pi_sel, p2_pi_sel)
        cmd.delete("tmp_dist")
        
        dist_pi1_cation2 = cmd.distance("tmp_dist", p1_pi_sel, p2_cation_sel)
        cmd.delete("tmp_dist")
        
        dist_cation1_pi2 = cmd.distance("tmp_dist", p1_cation_sel, p2_pi_sel)
        cmd.delete("tmp_dist")
        
        dist_cation1_cation2 = cmd.distance("tmp_dist", p1_cation_sel, p2_cation_sel)
        cmd.delete("tmp_dist")
        
        return {
            'pi1-pi2': round(dist_pi1_pi2, 2),
            'pi1-cation2': round(dist_pi1_cation2, 2),
            'cation1-pi2': round(dist_cation1_pi2, 2),
            'cation1-cation2': round(dist_cation1_cation2, 2)
        }
    except Exception as e:
        print(f"    ✗ 计算距离失败: {e}")
        traceback.print_exc()
        return None

def classify_saltbridge(distances):
    """根据距离分类盐桥"""
    pi1_pi2 = distances['pi1-pi2']
    pi1_cation2 = distances['pi1-cation2']
    cation1_pi2 = distances['cation1-pi2']
    cation1_cation2 = distances['cation1-cation2']
    
    if pi1_pi2 < pi1_cation2 and cation1_cation2 < cation1_pi2:
        return 'cla'
    elif pi1_pi2 > pi1_cation2 and cation1_cation2 > cation1_pi2:
        return 'nocla'
    else:
        return 'other'

def visualize_saltbridge(protein1_name, sb1_info, protein2_name, sb2_info, 
                        distances, category, pse_file):
    """生成PSE文件用于可视化"""
    try:
        cmd.hide("everything")
        cmd.show("cartoon")
        cmd.color("gray80", "all")
        
        # 选择残基（不带name CA，显示整个残基）
        p1_pi_sel, _ = try_select_with_alt(
            protein1_name, sb1_info['pi_chain'],
            sb1_info['pi_resi'], sb1_info['pi_resn']
        )
        p1_cation_sel, _ = try_select_with_alt(
            protein1_name, sb1_info['cation_chain'],
            sb1_info['cation_resi'], sb1_info['cation_resn']
        )
        p2_pi_sel, _ = try_select_with_alt(
            protein2_name, sb2_info['pi_chain'],
            sb2_info['pi_resi'], sb2_info['pi_resn']
        )
        p2_cation_sel, _ = try_select_with_alt(
            protein2_name, sb2_info['cation_chain'],
            sb2_info['cation_resi'], sb2_info['cation_resn']
        )
        
        # 去掉 "and name CA" 用于显示
        p1_pi_sel = p1_pi_sel.replace(" and name CA", "")
        p1_cation_sel = p1_cation_sel.replace(" and name CA", "")
        p2_pi_sel = p2_pi_sel.replace(" and name CA", "")
        p2_cation_sel = p2_cation_sel.replace(" and name CA", "")
        
        # 显示残基
        cmd.show("sticks", f"{p1_pi_sel} or {p1_cation_sel}")
        cmd.show("sticks", f"{p2_pi_sel} or {p2_cation_sel}")
        
        # 着色
        cmd.color("red", p1_pi_sel)
        cmd.color("blue", p1_cation_sel)
        cmd.color("orange", p2_pi_sel)
        cmd.color("cyan", p2_cation_sel)
        
        # 显示CA原子之间的距离
        cmd.distance("dist_pi1_pi2", 
                    f"{p1_pi_sel} and name CA", 
                    f"{p2_pi_sel} and name CA")
        cmd.distance("dist_pi1_cation2", 
                    f"{p1_pi_sel} and name CA", 
                    f"{p2_cation_sel} and name CA")
        cmd.distance("dist_cation1_pi2", 
                    f"{p1_cation_sel} and name CA", 
                    f"{p2_pi_sel} and name CA")
        cmd.distance("dist_cation1_cation2", 
                    f"{p1_cation_sel} and name CA", 
                    f"{p2_cation_sel} and name CA")
        
        # 标签
        cmd.label(f"{p1_pi_sel} and name CA", "'P1-pi'")
        cmd.label(f"{p1_cation_sel} and name CA", "'P1-cat'")
        cmd.label(f"{p2_pi_sel} and name CA", "'P2-pi'")
        cmd.label(f"{p2_cation_sel} and name CA", "'P2-cat'")
        
        cmd.zoom("all", 5)
        cmd.save(pse_file)
        print(f"    PSE文件已保存: {pse_file}")
        
    except Exception as e:
        print(f"生成PSE文件失败: {e}")

def process_family_file(family_id, input_file, ref_dict, protein_base_dir, 
                       output_dirs, generate_pse=False, pse_dir=None):
    """处理单个家族的盐桥文件"""
    print(f"\n{'='*60}")
    print(f"处理家族 {family_id}")
    print(f"{'='*60}")
    
    if family_id not in ref_dict:
        print(f"警告：家族 {family_id} 没有参考蛋白，跳过")
        return
    
    ref_protein = ref_dict[family_id]
    ref_protein_path = get_protein_path(ref_protein, family_id, protein_base_dir)
    
    if not os.path.exists(ref_protein_path):
        print(f"错误：参考蛋白文件不存在: {ref_protein_path}")
        return
    
    # 初始化PyMOL
    cmd.reinitialize()
    cmd.feedback("disable", "all", "everything")
    
    # 加载参考蛋白
    try:
        print(f"正在加载参考蛋白: {ref_protein}")
        cmd.load(ref_protein_path, "reference")
        print(f"✓ 参考蛋白加载成功")
    except Exception as e:
        print(f"✗ 加载参考蛋白失败: {e}")
        return
    
    # 读取输入文件
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            content = f.read()
    except Exception as e:
        print(f"✗ 读取文件失败: {e}")
        return
    
    # 准备输出内容
    results = {
        'all': [],
        'cla': [],
        'nocla': [],
        'other': [],
        'error': []
    }
    
    # 按组处理盐桥对
    groups = re.split(r'\n(?=\S+\s+与\s+\S+\s*:)', content)
    
    group_count = 0
    success_count = 0
    failed_count = 0
    error_count = 0
    
    for group in groups:
        group = group.strip()
        if not group:
            continue
        
        lines = group.split('\n')
        if len(lines) < 3:
            continue
        
        group_count += 1
        header = lines[0].strip()
        
        match = re.match(r'(\S+)\s+与\s+(\S+)\s*:', header)
        if not match:
            continue
        
        protein1_name = match.group(1)
        protein2_name = match.group(2)
        
        print(f"\n[{group_count}] 处理: {protein1_name} <-> {protein2_name}")
        
        if len(lines) < 3:
            continue
        
        sb1_line = lines[1].strip()
        sb2_line = lines[2].strip()
        
        sb1_info = parse_saltbridge_line(sb1_line)
        sb2_info = parse_saltbridge_line(sb2_line)
        
        if not sb1_info or not sb2_info:
            print(f"  ✗ 解析盐桥信息失败")
            failed_count += 1
            continue
        
        protein1_path = get_protein_path(protein1_name, family_id, protein_base_dir)
        protein2_path = get_protein_path(protein2_name, family_id, protein_base_dir)
        
        if not os.path.exists(protein1_path):
            print(f"  ✗ 蛋白质文件不存在: {protein1_path}")
            failed_count += 1
            continue
        if not os.path.exists(protein2_path):
            print(f"  ✗ 蛋白质文件不存在: {protein2_path}")
            failed_count += 1
            continue
        
        try:
            print(f"  加载蛋白质...")
            cmd.load(protein1_path, protein1_name)
            cmd.load(protein2_path, protein2_name)
            
            print(f"  叠合到参考蛋白...")
            cmd.super(protein1_name, "reference")
            cmd.super(protein2_name, "reference")
            
            print(f"  计算CA距离...")
            distances = calculate_ca_distances(protein1_name, sb1_info, 
                                              protein2_name, sb2_info)
            
            if distances:
                if 'error' in distances:
                    print(f"  ✗ 错误：选中多个原子")
                    result_text = f"{header} ERROR: multiple atoms selected\n{sb1_line}\n{sb2_line}\n"
                    results['error'].append(result_text)
                    error_count += 1
                else:
                    dist_str = (f"pi1-pi2:{distances['pi1-pi2']},"
                               f"pi1-cation2:{distances['pi1-cation2']},"
                               f"cation1-pi2:{distances['cation1-pi2']},"
                               f"cation1-cation2:{distances['cation1-cation2']}")
                    new_header = f"{header} {dist_str}"
                    
                    category = classify_saltbridge(distances)
                    
                    result_text = f"{new_header}\n{sb1_line}\n{sb2_line}\n"
                    results['all'].append(result_text)
                    results[category].append(result_text)
                    
                    print(f"  ✓ 距离: {dist_str}")
                    print(f"  ✓ 分类: {category}")
                    
                    success_count += 1
                    
                    if generate_pse and pse_dir:
                        pse_subdir = os.path.join(pse_dir, category)
                        os.makedirs(pse_subdir, exist_ok=True)
                        pse_file = os.path.join(pse_subdir, 
                                               f"{family_id}_{protein1_name}_{protein2_name}.pse")
                        visualize_saltbridge(protein1_name, sb1_info, 
                                           protein2_name, sb2_info, 
                                           distances, category, pse_file)
            else:
                failed_count += 1
            
            # 立即删除蛋白，避免内存占用
            print(f"  清理蛋白质对象...")
            cmd.delete(protein1_name)
            cmd.delete(protein2_name)
            
        except Exception as e:
            print(f"  ✗ 处理失败: {e}")
            traceback.print_exc()
            failed_count += 1
            try:
                cmd.delete(protein1_name)
            except:
                pass
            try:
                cmd.delete(protein2_name)
            except:
                pass
            continue
    
    # 写入输出文件
    print(f"\n{'='*60}")
    print(f"写入结果文件...")
    for category in ['all', 'cla', 'nocla', 'other', 'error']:
        if results[category]:
            output_dir = output_dirs[category]
            os.makedirs(output_dir, exist_ok=True)
            output_file = os.path.join(output_dir, f"{family_id}.txt")
            
            try:
                with open(output_file, 'w', encoding='utf-8') as f:
                    f.write('\n'.join(results[category]))
                print(f"  ✓ {category:6s}: {output_file} ({len(results[category])} 组)")
            except Exception as e:
                print(f"  ✗ 写入失败 {output_file}: {e}")
    
    cmd.delete("all")
    cmd.reinitialize()
    
    print(f"\n{'='*60}")
    print(f"家族 {family_id} 处理完成")
    print(f"总组数:   {group_count}")
    print(f"成功:     {success_count}")
    print(f"错误:     {error_count} (选中多个原子)")
    print(f"失败:     {failed_count}")
    print(f"{'='*60}")

def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='盐桥保守性分析工具',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例:
  python %(prog)s -f 4007548 -i /path/to/input -o /path/to/output -r reference.txt -p /path/to/proteins
  python %(prog)s -f 4007548 -i /path/to/input -o /path/to/output -r reference.txt -p /path/to/proteins --pse --pse-dir /path/to/pse
        '''
    )
    
    parser.add_argument('-f', '--family', required=True, help='家族ID')
    parser.add_argument('-i', '--input-dir', required=True, help='输入目录路径')
    parser.add_argument('-o', '--output-dir', required=True, help='输出目录路径')
    parser.add_argument('-r', '--reference', required=True, help='参考蛋白列表文件路径')
    parser.add_argument('-p', '--protein-dir', required=True, help='蛋白质文件基础目录')
    parser.add_argument('--pse', action='store_true', help='生成PSE可视化文件')
    parser.add_argument('--pse-dir', default=None, help='PSE文件保存目录')
    
    args = parser.parse_args()
    
    if args.pse and not args.pse_dir:
        parser.error("使用 --pse 时必须指定 --pse-dir")
    
    output_dirs = {
        'all': os.path.join(args.output_dir, 'all'),
        'cla': os.path.join(args.output_dir, 'cla'),
        'nocla': os.path.join(args.output_dir, 'nocla'),
        'other': os.path.join(args.output_dir, 'other'),
        'error': os.path.join(args.output_dir, 'error')
    }
    
    for dir_path in output_dirs.values():
        os.makedirs(dir_path, exist_ok=True)
    
    if args.pse and args.pse_dir:
        os.makedirs(args.pse_dir, exist_ok=True)
    
    ref_dict = load_reference_list(args.reference)
    if not ref_dict:
        print("错误：无法加载参考列表")
        sys.exit(1)
    
    input_file = os.path.join(args.input_dir, f"{args.family}.txt")
    
    if not os.path.exists(input_file):
        print(f"错误：输入文件不存在: {input_file}")
        sys.exit(1)
    
    try:
        process_family_file(args.family, input_file, ref_dict, 
                           args.protein_dir, output_dirs, 
                           args.pse, args.pse_dir)
        print(f"\n✓✓✓ 全部处理完成！✓✓✓\n")
    except Exception as e:
        print(f"\n✗✗✗ 处理出错: {e} ✗✗✗\n")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()