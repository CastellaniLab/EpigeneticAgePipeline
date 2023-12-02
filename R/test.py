def calculate_indent_size(line):
    if not line.strip():  # Check if the line is empty or contains only whitespace
        return 0
    else:
        return len(line) - len(line.lstrip())

def output_indent_sizes(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

        print("Indentation sizes for each non-empty line:")
        for i, line in enumerate(lines, start=1):
            indent_size = calculate_indent_size(line)
            if indent_size > 0:
                print(f"Line {i}: {indent_size}")
                if indent_size % 4 != 0:
                    print(f"   - Indentation size is not a multiple of 4")

# Replace 'your_file.py' with the path to your Python source file
output_indent_sizes('file1.R')
