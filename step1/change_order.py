def sort_gff3_chr1_to_chr8(input_gff3, output_gff3):
    # 存储所有行的列表
    gff3_lines = []

    with open(input_gff3, 'r') as infile:
        for line in infile:
            line = line.strip()
            if line:  # 只添加非空行
                gff3_lines.append(line)

    # 过滤出注释行和数据行
    header_lines = [line for line in gff3_lines if line.startswith('#')]
    data_lines = [line for line in gff3_lines if not line.startswith('#')]

    # 分离 Chr1 到 Chr8 和其他染色体
    chr1_to_chr8 = []
    others = []

    for line in data_lines:
        seq_id = line.split('\t')[0]
        if seq_id in [f'Chr{i}' for i in range(1, 9)]:  # 如果是 Chr1 到 Chr8
            chr1_to_chr8.append(line)
        else:
            others.append(line)

    # 对 Chr1 到 Chr8 进行排序
    chr1_to_chr8_sorted = sorted(chr1_to_chr8, key=lambda x: int(x.split('\t')[0][3:]))

    # 写入排序后的GFF3文件
    with open(output_gff3, 'w') as outfile:
        # 写入注释行
        for header in header_lines:
            outfile.write(header + '\n')
        # 写入排序后的 Chr1 到 Chr8 的数据行
        for sorted_line in chr1_to_chr8_sorted:
            outfile.write(sorted_line + '\n')
        # 写入其他染色体的数据行
        for other_line in others:
            outfile.write(other_line + '\n')

# 使用示例
sort_gff3_chr1_to_chr8('Medicago_arabica_rename_reorder.gff3', 'Medicago_arabica.gff3')

