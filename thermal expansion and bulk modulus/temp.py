from pathlib import Path
import numpy as np
import re
import os

script_dir = os.path.dirname(__file__)

def inject_temperature_block_with_transform(temperature, source_path, target_path):
    with open(source_path, 'r') as file:
        lines = file.readlines()

    pattern = re.compile(r"# Temperature:\s+([0-9.]+)")
    blocks = []
    current_block = []
    current_temp = None

    for line in lines:
        match = pattern.match(line)
        if match:
            if current_block and current_temp is not None:
                blocks.append((float(current_temp), current_block))
            current_temp = match.group(1)
            current_block = [line]
        else:
            current_block.append(line)

    if current_block and current_temp is not None:
        blocks.append((float(current_temp), current_block))

    selected_block = min(blocks, key=lambda x: abs(x[0] - temperature))
    selected_lines = selected_block[1]

    transformed_lines = []
    for line in selected_lines:
        if line.strip().startswith("#") or line.strip() == "":
            transformed_lines.append(line)
        else:
            try:
                parts = line.strip().split()
                x = float(parts[0])
                y = float(parts[1])
                x_new = x**(1/3)
                transformed_lines.append(f"{x_new:.12f}      {y:.12f}\n")
            except:
                transformed_lines.append(line)

    with open(target_path, 'w') as out_file:
        out_file.writelines(transformed_lines)

    return selected_block[0]