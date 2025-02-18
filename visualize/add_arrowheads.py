import sys
from lxml import etree

SVG_NS = "http://www.w3.org/2000/svg"
NSMAP = {"svg": SVG_NS}  # Namespace mapping for XPath

def add_arrowheads_to_svg(input_svg, output_svg):
    # Parse the SVG file
    tree = etree.parse(input_svg)
    root = tree.getroot()

    # Find <path> elements (now using proper namespaces)
    paths = root.findall(".//svg:path", namespaces=NSMAP)

    if not paths:
        print("⚠️ No <path> elements found. Check your SVG file!")
        return

    print(f"✅ Found {len(paths)} <path> elements. Applying arrowheads...")

    # Ensure there's a <defs> section for markers
    defs = root.find("svg:defs", namespaces=NSMAP)
    if defs is None:
        defs = etree.Element(f"{{{SVG_NS}}}defs")
        root.insert(0, defs)  # Insert at the beginning

    # Define a sharp arrowhead marker
    marker_id = "ConcaveTriangle"
    marker = etree.Element(f"{{{SVG_NS}}}marker", {
        "id": marker_id,
        "viewBox": "0 0 10 10",
        "refX": "5",  # Moves arrowhead forward
        "refY": "0",
        "markerWidth": "1",
        "markerHeight": "1",
        "orient": "auto"
    })
    
    arrow_path = etree.Element(f"{{{SVG_NS}}}path", {
        "d": "M 0 0 L 10 5 L 0 10 z",  # Sharp arrow shape
        "fill": "black"
    })
    
    marker.append(arrow_path)
    defs.append(marker)

    # Apply arrowheads to each <path>
    for path in paths:
        path.set("marker-end", f"url(#{marker_id})")
        print(f"🖋 Applied arrowhead to path: {etree.tostring(path).decode()}")

    # Save the modified SVG
    tree.write(output_svg, pretty_print=True, xml_declaration=True, encoding="utf-8")
    print(f"✅ SVG successfully updated and saved to {output_svg}")

# Example usage
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python add_arrowheads.py input.svg output.svg")
    else:
        input_svg = sys.argv[1]
        output_svg = sys.argv[2]
        add_arrowheads_to_svg(input_svg, output_svg)
