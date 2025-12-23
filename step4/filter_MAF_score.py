import re

# 设置你想要的最低分数阈值
score_threshold = 20.0

# 打开MAF文件
with open('B_ZM4_6_arabidopsis_aradu_chickpea_pea_phase_populus_soybean_vigna_vitis.maf', 'r') as maf_file:
    # 为了存储过滤后的结果创建一个新文件
    with open('filtered.maf', 'w') as filtered_file:
        # 初始化一个标记，指示当前是否在一个比对块内部
        in_block = False
        # 初始化一个标记，指示当前比对块是否满足条件
        block_accepted = False

        # 逐行遍历MAF文件
        for line in maf_file:
            # 如果行以#开头，说明是注释，直接写入新文件
            if line.startswith('#'):
                filtered_file.write(line)
                continue

            # 检查是否是一个新的比对块的开始
            if line.startswith('a score='):
                # 提取分数并转换为浮点数
                score = float(re.search(r'score=([\d\.\-]+)', line).group(1))
                # 重置比对块的接受状态为假
                block_accepted = False
                # 如果分数未达到阈值，则跳过当前比对块
                if score < score_threshold:
                    in_block = False
                    continue
                else:
                    in_block = True
                    # 预阅读下一行以检查序列名
                    next_line = next(maf_file)
                    # 如果下一行以`s`开头且紧接着的是Chr6(参考基因组开头的序列)，那么比对块被接受
                    if next_line.startswith('s') and re.match(r's Chr6', next_line):
                        block_accepted = True
                        # 将比对块开始行和序列行写入新文件
                        filtered_file.write(line)
                        filtered_file.write(next_line)
                    else:
                        # 如果不是Chr6，跳过当前比对块
                        continue
            
            # 如果当前处于比对块内，且比对块被接受，写入余下的行
            elif in_block and block_accepted:
                filtered_file.write(line)
            # 如果行不是以#开头和不是比对块的一部分，则不做任何事
