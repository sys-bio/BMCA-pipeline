import sys
from lxml import etree

def flip_rects_with_text(input_svg, output_svg):
    # Parse the input SVG
    tree = etree.parse(input_svg)
    root = tree.getroot()
    
    # Namespace for proper XPath search
    ns = {"svg": "http://www.w3.org/2000/svg"}  # SVG namespace
    
    # Iterate through all <g> elements, which contain both <rect> and <text>
    for group in root.findall(".//svg:g", namespaces=ns):
        rect = group.find(".//svg:rect", namespaces=ns)
        text = group.find(".//svg:text", namespaces=ns)

        if rect is not None and text is not None:
            # Get the position (x, y) of the rect
            x = float(rect.attrib['x'])
            y = float(rect.attrib['y'])
            width = float(rect.attrib['width'])
            height = float(rect.attrib['height'])

            # Compute the center of the rectangle
            center_x = x + width / 2
            center_y = y + height / 2

            # Rotate the rect around its center
            rect_transform = f"rotate(90 {center_x} {center_y})"
            rect.attrib['transform'] = rect_transform

            # Rotate the text around its own center
            text_transform = f"rotate(90 {center_x} {center_y})"
            text.attrib['transform'] = text_transform

            print(f"Rotated rect and text at ({center_x}, {center_y}) by 90 degrees.")

    # Save the modified SVG to the output file
    tree.write(output_svg, pretty_print=True, xml_declaration=True, encoding="utf-8")
    print(f"✅ Updated SVG saved as {output_svg}")

# Execute the function when running the script
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python flip_rects_with_text.py input.svg output.svg")
    else:
        flip_rects_with_text(sys.argv[1], sys.argv[2])
