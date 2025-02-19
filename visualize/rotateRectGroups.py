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

        # Only consider groups that contain both a <rect> and <text>
        if rect is not None and text is not None:
            # Get the position (x, y) of the rect
            x = float(rect.attrib['x'])
            y = float(rect.attrib['y'])
            width = float(rect.attrib['width'])
            height = float(rect.attrib['height'])

            # Apply the transform (rotate the rect 90 degrees around the top-left corner)
            transform = f"rotate(90 {x} {y})"
            rect.attrib['transform'] = transform

            # Rotate the text (apply the same rotation)
            text.attrib['transform'] = transform

            # After rotating, adjust the text position to keep it inside the rotated rectangle
            new_text_x = x + height  # New x position for text after rotation
            new_text_y = y - width   # New y position for text after rotation

            text.attrib['x'] = str(new_text_x)
            text.attrib['y'] = str(new_text_y)

            # Debugging output
            print(f"Rotated rect at ({x}, {y}) by 90 degrees and adjusted text position.")

    # Save the modified SVG to the output file
    tree.write(output_svg, pretty_print=True, xml_declaration=True, encoding="utf-8")
    print(f"✅ Updated SVG saved as {output_svg}")

# Execute the function when running the script
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python flip_rects_with_text.py input.svg output.svg")
    else:
        flip_rects_with_text(sys.argv[1], sys.argv[2])
