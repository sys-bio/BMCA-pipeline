import os
import nbformat

# Set parameters
target_directory = "./notebooks"
match_line = "ppc_vi = pm.sample_posterior_predictive(trace, random_seed=random_range)"
replacement_lines = ["        ppc_vi = pm.sample_posterior_predictive(trace, random_seed=SEED)"]

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
        replaced = False

        for line in lines:
            if not replaced and line.strip() == match_line:
                print(f"  ➤ Replacing first match in cell {i}")
                new_lines.extend(replacement_lines)
                replaced = True
                modified = True
            else:
                new_lines.append(line)

        if replaced:
            cell["source"] = "\n".join(new_lines)
            break  # stop after first match in notebook

    if modified:
        with open(filepath, "w", encoding="utf-8") as f:
            nbformat.write(nb, f)
        print(f"  ✔ Notebook updated")
    else:
        print(f"  ✖ No matching line found")

# Walk through all notebooks
for root, dirs, files in os.walk(target_directory):
    for file in files:
        if file.endswith(".ipynb"):
            full_path = os.path.join(root, file)
            process_notebook(full_path)
