import sys

def fix_gff(input_file, output_file):
    with open(input_file) as fin, open(output_file, 'w') as fout:
        # Write GFF3 header
        fout.write("##gff-version 3\n")
        
        # Process each line
        for line in fin:
            if line.startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            if len(fields) != 9:
                continue
                
            # Ensure proper phase for CDS features
            if fields[2] == 'CDS' and fields[7] == '.':
                fields[7] = '0'
                
            # Add missing required attributes
            attributes = fields[8].split(';')
            attr_dict = dict(attr.split('=') if '=' in attr else (attr, '') for attr in attributes)
            
            # Ensure ID and Parent attributes are properly formatted
            if 'ID' not in attr_dict:
                attr_dict['ID'] = f"{fields[2]}_{fields[3]}_{fields[4]}"
            if fields[2] in ['CDS', 'exon'] and 'Parent' not in attr_dict:
                attr_dict['Parent'] = attr_dict['ID']
                
            # Reconstruct attributes field
            fields[8] = ';'.join(f"{k}={v}" for k, v in attr_dict.items())
            
            # Write the modified line
            fout.write('\t'.join(fields) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python fix_gff.py input.gff output.gff")
        sys.exit(1)
        
    fix_gff(sys.argv[1], sys.argv[2])
