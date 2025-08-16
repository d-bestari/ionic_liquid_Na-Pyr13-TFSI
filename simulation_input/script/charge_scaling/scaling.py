#!/usr/bin/python3

import argparse
import os

def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def scale_charges(lines, scaling_factor): 
    scaled_lines = [] 
    in_atoms = False 
    in_atomtypes = False
    atom_count = 0 
    atomtypes_count = 0
    current_molecule = '' 
    skip_molecules = set() 

    for line in lines: 
        stripped = line.strip() 

        if stripped.startswith('['): 
            in_atoms = '[ atoms ]' in stripped 
            in_atomtypes = '[ atomtypes ]' in stripped
            current_molecule = ''
            scaled_lines.append(line) 
            continue 

        if in_atoms and not stripped.startswith(';') and stripped: 
            parts = line.split(';', 1) 
            data = parts[0].split() 
            comment = f";{parts[1]}" if len(parts) > 1 else ''

  
            try:
                if len(data) >= 7 and is_float(data[6]): 
                    original_charge = float(data[6])
                    scaled_charge = original_charge * scaling_factor
                    data[6] = f"{scaled_charge:.6f}"
                    scaled_line = '\t'.join(data) 
                if comment:
                        scaled_line += f"   {comment}"
                scaled_lines.append(scaled_line + "\n") 
                atom_count += 1
                continue
            except ValueError:
                pass
    
        elif in_atomtypes and not stripped.startswith(';') and stripped:
            parts = line.split(';', 1)
            data = parts[0].split()
            comment = f";{parts[1]}" if len(parts) > 1 else ''

            try:
                if len(data) >= 4 and is_float(data[3]):
                    original_charge = float(data[3])
                    scaled_charge = original_charge * scaling_factor
                    data[3] = f"{scaled_charge:.4f}"

                    scaled_line = '\t'.join(data)
                    if comment:
                        scaled_line += f"   {comment}"

                    scaled_lines.append(scaled_line + "\n")
                    atomtypes_count += 1
                    continue
            except ValueError:
                pass

        scaled_lines.append(line)
    print("Scaling charges completed:")
    print(f"- Total atoms scaled in [ atoms ]: {atom_count}")
    print(f"- Total atomtypes scaled in [ atomtypes ]: {atomtypes_count}")
    return scaled_lines 

# --------------------------
# MAIN SCRIPT
# --------------------------

if __name__ == "__main__": 
    parser = argparse.ArgumentParser(description="Scale atomic charges in a .top/.itp file.") 
    parser.add_argument("input_file", help="Input .top or .itp file") 
    parser.add_argument("scaling_factor", type=float, help="Scaling factor for charges") 

    args = parser.parse_args() 
    input_file = args.input_file
    scaling_factor = args.scaling_factor

    base, ext = os.path.splitext(input_file) 
    output_file = f"{base}_scaled_{scaling_factor:.4f}{ext}" #{base}
    
    print(f"Reading input file: {input_file}")
    with open(input_file, 'r') as f: 
        lines = f.readlines()
    
    print(f"Scaling charges with factor: {scaling_factor}")
    scaled_lines = scale_charges(lines, scaling_factor) 
    
    print(f"Writing scaled output to: {output_file}")
    with open(output_file, 'w') as f:
        f.writelines(scaled_lines)

    print("Done!")
