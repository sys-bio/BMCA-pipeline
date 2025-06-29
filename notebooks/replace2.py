import os
import nbformat

# Set parameters
target_directory = "./notebooks"
match_lines = ["def run_prior_predictive(BMCA_obj):", "random_range = np.random.default_rng(SEED)"]
replacement_lines = ["def run_prior_predictive(BMCA_obj):", ""]  # Change if different

def process_notebook(filepath):
    print(f"\nChecking: {filepath}")
    with open(filepath, "r", encoding="utf-8") as f:
        nb = nbformat.read(f, as_version=4)

    modified = False

    for i, cell in enumerate(nb.cells):
        if cell["cell_type"] != "code":
            continue

        lines = cell["source"].splitlines()
        new_lines = []
        j = 0
        while j < len(lines):
            # Check for 2-line match
            if (
                j + 1 < len(lines) and
                lines[j].strip() == match_lines[0] and
                lines[j + 1].strip() == match_lines[1]
            ):
                print(f"  ➤ Replacing lines {j} and {j+1} in cell {i}")
                new_lines.extend(replacement_lines)
                j += 2  # skip both matched lines
                modified = True
                break  # stop after first replacement
            else:
                new_lines.append(lines[j])
                j += 1

        if modified:
            new_lines.extend(lines[j:])  # include remaining lines if any
            cell["source"] = "\n".join(new_lines)
            break

    if modified:
        with open(filepath, "w", encoding="utf-8") as f:
            nbformat.write(nb, f)
        print(f"  ✔ Notebook updated")
    else:
        print(f"  ✖ No matching pair of lines found")

# Walk through all notebooks
for root, dirs, files in os.walk(target_directory):
    for file in files:
        if file.endswith(".ipynb"):
            full_path = os.path.join(root, file)
            process_notebook(full_path)
