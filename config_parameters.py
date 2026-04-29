import re
import csv
import os
from collections.abc import Hashable
import ast

# Read the config file
config_path = 'config.py'
with open(config_path, 'r', encoding='utf-8') as f:
    config_content = f.read()

# Extract all top-level assignments
assignments = re.findall(r'^([A-Z_]+)\s*=\s*(.+)$', config_content, re.MULTILINE)
output = []

for var_name, value_str in assignments:
    # Safely evaluate value
    try:
        if value_str.strip().startswith('['):
            value = ast.literal_eval(value_str)
        else:
            value = eval(value_str)
        output.append((var_name, str(value)))
    except:
        # Fallback for strings or complex cases
        output.append((var_name, value_str.strip()))

# Write CSV
with open('config_vars.csv', 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerow(['Variable', 'Value'])
    writer.writerows(output)

print("CSV exported: config_vars.csv")
