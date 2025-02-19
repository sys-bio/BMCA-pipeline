import inkex
from lxml import etree

print("Metabolic Map Generator script is running...")

class MetabolicMap(inkex.EffectExtension):
    def add_arguments(self, pars):
        pars.add_argument("--width", type=float, default=8.5*96, help="Document width in inches")
        pars.add_argument("--height", type=float, default=9*96, help="Document height in inches")
    
    def effect(self):
        # Check if an input file was provided; if not, create a new blank document
        if not self.options.input_file:
            self.document = inkex.ElementTree(inkex.elements.Svg())
            print("No input file detected. Using a blank SVG document.")
        else:
            print(f"Using input file: {self.options.input_file}")
        
        layer = self.svg.get_current_layer()
        if layer is None:
            layer = inkex.Layer.new("Metabolic Map")
            self.svg.add(layer)
        
        print(f"DEBUG: Using SVG layer: {layer.get_id()}")  # Debugging

        # Ensure correct document dimensions
        svg_width = self.options.width  # Should be 816px (8.5 inches * 96 DPI)
        svg_height = self.options.height  # Should be 864px (9 inches * 96 DPI)
        
        print(f"DEBUG: Setting SVG size to {svg_width} x {svg_height} px")

        # Explicitly set width, height, and viewBox
        self.svg.set("width", f"{svg_width}px")  # Force width in pixels
        self.svg.set("height", f"{svg_height}px")  # Force height in pixels
        # self.svg.set("viewBox", f"0 0 {svg_width} {svg_height}")  # Ensure correct scaling
        self.svg.set("viewBox", "0 0 612 648")  # Ensure correct scaling
        
        # Elements
        metabolites = {
            "A": (100, 100), "B": (200, 100), "C": (300, 100), "D": (400, 100), "E": (500, 100),
            "F": (100, 200), "G": (200, 200), "H": (300, 200), "I": (400, 200), "J": (500, 200),
            "L": (100, 300), "M": (200, 300), "N": (300, 300), "O": (400, 300), "P": (500, 300),
            "Q": (200, 400), "R": (300, 400)
        }

        reactions = {
            "v1": (100, 50), "v2": (200, 50), "v3": (300, 50), "v4": (400, 50), "v5": (500, 50),
            "v6": (100, 150), "v7": (200, 150), "v8": (300, 150), "v9": (400, 150), "v10": (500, 150),
            "v11": (100, 250), "v12": (200, 250), "v13": (300, 250), "v14": (400, 250), "v15": (500, 250),
            "v16": (100, 350), "v17": (200, 350), "v18": (300, 350), "v19": (400, 350)
        }


        # 🚀 Ensure all metabolites are drawn
        for metabolite, (x, y) in metabolites.items():
            print(f"DEBUG: Adding metabolite {metabolite} at ({x}, {y})")  # Debugging
            
            # Estimate text width based on character count
            char_width = 6.5  # Approximate width of a character in pixels
            text_width = len(metabolite) * char_width  # Adjust box width based on text length
            box_padding = 5  # Padding around the text
            box_width = text_width + box_padding  # Total box width
            box_height = 12  # Fixed height for consistency

            # Create rectangle (box) around text
            rect = etree.Element('rect', {
                'x': str(x - box_width / 2),  # Center the box around text
                'y': str(y - box_height / 2),
                'width': str(box_width),
                'height': str(box_height),
                'fill': 'none',
                'stroke': 'black',
                'ry': '2'  # Rounded corners for better aesthetics
            })

            # Create text element
            text = etree.Element('text', {
                'x': str(x),
                'y': str(y + 3),
                'font-size': '10',
                'font-family': 'Arial',
                'text-anchor': 'middle'
            })
            text.text = metabolite

            # Group the rect and text together
            group = etree.Element('g')  # Create a group element
            group.append(rect)  # Add the rect to the group
            group.append(text)  # Add the text to the group

            # Add the group to the layer
            layer.append(group)

        # 🚀 Ensure all reactions are drawn
        for reaction, (x, y) in reactions.items():
            # Create text element
            rxn_label = etree.Element('text', {
                'x': str(x),
                'y': str(y + 3),
                'font-size': '8',
                'font-family': 'Arial',
                'text-anchor': 'middle',
                'fill': 'red',
                'font-weight': 'bold',
                'text-anchor': 'middle'
            })
            rxn_label.text = reaction
            layer.append(rxn_label)

        # save svg to new file
        output_filename = "output_metabolic_map.svg"
        self.document.write(output_filename)
        print(f"SVG successfully saved to {output_filename}")
        
if __name__ == '__main__':
    MetabolicMap().run()
