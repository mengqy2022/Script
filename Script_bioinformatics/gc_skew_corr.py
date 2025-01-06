import argparse

def calculate_gc_skew(sequence):
    # 去除换行符并简化字符串
    sequence = sequence.replace("\n", "").replace(">", "")
    
    G = 0
    C = 0
    gc_skew_values = []
    
    # 计算 G 和 C 的数量
    for char in sequence:
        if char == 'G':
            G += 1
        elif char == 'C':
            C += 1
        
        # 防止除以零
        if G + C > 0:
            gc_skew = (G - C) / (G + C)
        else:
            gc_skew = 0
        
        gc_skew_values.append(gc_skew)
    
    return gc_skew_values

def find_intersection_coordinates(gc_skew_values):
    coordinates = []
    
    for idx in range(1, len(gc_skew_values)):
        # 检查每一段从正变为负或从负变为正
        if (gc_skew_values[idx-1] >= 0 and gc_skew_values[idx] < 0) or (gc_skew_values[idx-1] < 0 and gc_skew_values[idx] >= 0):
            coordinates.append(idx)
    
    return coordinates

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', type=str, help='Input file containing the DNA sequence')
    
    args = parser.parse_args()

    # 读取序列文件
    with open(args.input_file, 'r') as file:
        sequence = file.read()
    
    # 计算 GC-skew
    gc_skew_values = calculate_gc_skew(sequence)

    # 找到相交坐标
    intersection_coordinates = find_intersection_coordinates(gc_skew_values)

    print("Intersection Coordinates:", intersection_coordinates)

if __name__ == "__main__":
    main()
