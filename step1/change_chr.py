def replace_chromosomes(gff_file, mapping_file, output_file):
    # 读取映射文件并创建字典
    mapping_dict = {}
    with open(mapping_file, 'r') as map_file:
        for line in map_file:
            original, new_name = line.strip().split('\t')
            mapping_dict[original] = new_name

    # 替换GFF文件中的染色体名称
    with open(gff_file, 'r') as gff_in, open(output_file, 'w') as gff_out:
        for line in gff_in:
            if line.startswith('#'):  # 如果是注释行，直接写入
                gff_out.write(line)
            else:
                fields = line.strip().split('\t')
                original_chromosome = fields[0]
                # 替换染色体名称
                new_chromosome = mapping_dict.get(original_chromosome, original_chromosome)
                fields[0] = new_chromosome
                gff_out.write('\t'.join(fields) + '\n')

# 使用示例
replace_chromosomes('Medicago_arabica.gff3', 'mapping.txt', 'Medicago_arabica_rename.gff3')

