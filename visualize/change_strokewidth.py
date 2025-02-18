import sys
from lxml import etree

def change_stroke_width(input_svg, output_svg):
    tree = etree.parse(input_svg)
    root = tree.getroot()
    ns = {"svg": "http://www.w3.org/2000/svg"}  # Namespace for proper search

    for path in root.findall(".//svg:path", namespaces=ns):
        if "stroke-width" in path.attrib:
            path.set("stroke-width", "0.353")

        if "style" in path.attrib:
            styles = dict(item.split(":") for item in path.attrib["style"].split(";") if ":" in item)
            styles["stroke-width"] = "0.353"
            path.set("style", ";".join(f"{k}:{v}" for k, v in styles.items()))

    tree.write(output_svg, pretty_print=True, xml_declaration=True, encoding="utf-8")
    print(f"✅ Updated SVG saved as {output_svg}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python change_path_thickness.py input.svg output.svg")
    else:
        change_stroke_width(sys.argv[1], sys.argv[2])
