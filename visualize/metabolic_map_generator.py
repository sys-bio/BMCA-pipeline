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
            "G6P": (100, 100), "F6P": (200, 100), "FDP": (300, 100), "GAP": (400, 100), "PYR": (500, 100),
            "ATP": (100, 200), "ADP": (200, 200), "NADH": (300, 200), "NAD": (400, 200), "OAA": (500, 200),
            "PEP": (600, 200), "NADP": (700, 200), "NADPH": (800, 200), "P": (100, 300), "BPG": (200, 300),
            "FAD": (300, 300), "FADH2": (400, 300), "AMP": (500, 300), "CAMP": (600, 300), "DAP": (700, 300),
            "PGA3": (800, 300), "PGA2": (100, 400), "GL6P": (200, 400), "PGN": (300, 400), "RU5P": (400, 400),
            "X5P": (500, 400), "R5P": (600, 400), "E4P": (700, 400), "S7P": (800, 400), "MAL": (100, 500),
            "Q": (200, 500), "QH2": (300, 500), "ICIT": (400, 500), "AKG": (500, 500), "SUCCOA": (600, 500),
            "SUC": (700, 500), "FUM": (800, 500), "GLX": (100, 600), "ACCOA": (200, 600), "ACEx": (300, 600),
            "ACP": (400, 600)
        }

        reactions = [
            ("vPGI", "G6P", "F6P"), ("vPFK", "ATP", "F6P"), ("vFBA", "FDP", "GAP"),
            ("vTPI", "DAP", "GAP"), ("vGDH", "GAP", "NADH"), ("vPGK", "ADP", "BPG"),
            ("vGPM", "PGA3", "PGA2"), ("vENO", "PGA2", "PEP"), ("vPYK", "ADP", "PEP"),
            ("vZWF", "G6P", "NADP"), ("vPGL", "GL6P", "PGN"), ("vGND", "NADP", "PGN"),
            ("vRPE", "RU5P", "X5P"), ("vRPI", "RU5P", "R5P"), ("vF6P_GAP_TAL", "GAP", "F6P"),
            ("vPPC", "PEP", "OAA"), ("vPCK", "ATP", "OAA"), ("vMAD", "MAL", "NADH"),
            ("vPDH", "NAD", "PYR"), ("vGLT", "ACCOA", "OAA"), ("vICD", "ICIT", "NADP"),
            ("vLPD", "AKG", "NADH"), ("vSK", "ADP", "SUCCOA"), ("vSDH", "FAD", "SUC"),
            ("vFUMA", "FUM", "MAL"), ("vMQO", "MAL", "Q"), ("vMDH", "QH2", "OAA"),
            ("vATP_syn", "ADP", "P"), ("vCYA", "ATP", "CAMP"), ("vDOS", "CAMP", "AMP"),
            ("vACK", "ACP", "ADP"), ("vACS", "ACEx", "ATP"), ("vPTA", "ACCOA", "P")
        ]
        
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
        for reaction, src, dst in reactions:
            if src in metabolites and dst in metabolites:
                src_x, src_y = metabolites[src]
                dst_x, dst_y = metabolites[dst]
                print(f"DEBUG: Drawing reaction {reaction} from {src} to {dst}")  # Debugging

                arrow = etree.Element('path', {
                    'd': f'M {src_x+4.3},{src_y+2.4} L {dst_x+4.3},{dst_y+2.4}',
                    'stroke': 'black', 'fill': 'none', 'marker-end': 'url(#Arrow1Lend)'
                })
                label = etree.Element('text', {
                    'x': str((src_x + dst_x) / 2), 'y': str((src_y + dst_y) / 2 - 5),
                    'font-size': '10', 'font-family': 'Arial', 'font-weight': 'bold',
                    'fill': 'red', 'text-anchor': 'middle'
                })
                label.text = reaction
                # layer.append(arrow)
                layer.append(label)
            else:
                print(f"WARNING: Reaction {reaction} has missing metabolites: {src}, {dst}")

        # save svg to new file
        output_filename = "output_metabolic_map.svg"
        self.document.write(output_filename)
        print(f"SVG successfully saved to {output_filename}")
        
if __name__ == '__main__':
    MetabolicMap().run()
