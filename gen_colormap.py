import colorpy.thinfilm as thinfilm
import colorpy.illuminants as illuminants
import math
import xml.etree.ElementTree as ET

illuminant = illuminants.get_illuminant_D65()
illuminants.scale_illuminant (illuminant, 9.50)
n1 = 1.000293
n2 = 1.333
n3 = 1.000293
tf = thinfilm.thin_film(n1, n2, n3, 0)

def get_color_single(h):
    tf.thickness_nm = h
    tf.phase_factor = -2.0 * tf.thickness_nm * 2.0 * math.pi * n2
    return tf.illuminated_color(illuminant)

colormaps = ET.Element("ColorMaps")

colormap = ET.SubElement(colormaps, "ColorMap")
colormap.set("name", "ThinFilmThickness")
colormap.set("space", "RGB")
colormap.set("indexed", "0")

for i in range(1000):
    val = i * 1e-9
    r, g, b = get_color_single(i)
    point = ET.SubElement(colormap, "Point")
    point.set("x", f"{val:.9e}")
    point.set("r", f"{r:.6f}")
    point.set("g", f"{g:.6f}")
    point.set("b", f"{b:.6f}")

tree = ET.ElementTree(colormaps)
tree.write("thinfilm_colormap.xml", encoding="utf-8", xml_declaration=True)

print("✅ Color map saved as 'thinfilm_colormap.xml'")
