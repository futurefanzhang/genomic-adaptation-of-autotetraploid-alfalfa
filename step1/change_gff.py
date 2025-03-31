def flip_gff3(input_gff3, output_gff3, mapping_file):
    # 读取染色体长度映射文件并创建字典
    length_dict = {}
    with open(mapping_file, 'r') as map_file:
        for line in map_file:
            fields = line.strip().split('\t')
            if len(fields) >= 2:
                chromosome_name = fields[0]
                chromosome_length = int(fields[1])
                length_dict[chromosome_name] = chromosome_length

    # 打开GFF3文件进行处理
    with open(input_gff3, 'r') as infile, open(output_gff3, 'w') as outfile:
        for line in infile:
            # 去掉两端空白字符
            line = line.strip()
            if not line:  # 跳过空行
                continue
            
            if line.startswith('#'):  # 跳过注释行
                outfile.write(line + '\n')
                continue
            
            fields = line.split('\t')
            
            # 检查行是否包含至少8个字段
            if len(fields) < 8:
                print(f"警告：行格式不正确，跳过此行: {line}")
                continue
            
            seq_id = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]

            # 检查该染色体是否在长度字典中
            if seq_id in length_dict:
                new_length = length_dict[seq_id]
                # 计算新的坐标
                new_start = new_length - end + 1
                new_end = new_length - start + 1

                # 更新方向
                new_strand = '-' if strand == '+' else '+' if strand == '-' else strand

                # 写入新的GFF3文件
                fields[3] = str(new_start)
                fields[4] = str(new_end)
                fields[6] = new_strand
            
            # 写入输出文件
            outfile.write('\t'.join(fields) + '\n')

# 使用示例
flip_gff3('Medicago_arabica_rename.gff3', 'Medicago_arabica_rename_reorder.gff3', 'mapping.txt')

