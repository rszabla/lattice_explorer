![Lattice Explorer Overview](images/logo_chatgpt.png)

Lattice Explorer is a plugin for generating and visualizing crystal lattices directly within PyMOL. It automates the generation of symmetry mates and allows for easy exploration of crystal packing, unit cell boundaries, and lattice interactions. It is particularly useful for visually comparing different crystal lattices and generating high-quality, consistent figures for publication.

---
## Prerequisites
- [PyMOL v2.0+](https://pymol.org/)

---
## Basic Usage

1. Download [`lattice_explorer.py`](lattice_explorer.py).
2. Open PyMOL.
3. `File` > `Run Script...` and select `lattice_explorer.py`.
4. The Lattice Explorer window should open automatically.
   
   ![Lattice Explorer Overview](images/main_screen.png)

5. **Input:**
    * **Crystal structure from file:** Browse for a local PDB/CIF file.
    * **PDB code:** Enter a 4-character PDB code to fetch from the RCSB.
6. **Parameters:**
    * **Unit cell expansion:** Number of unit cells to expand in each direction (e.g., 1 creates a 3x3x3 block).
    * **Wrap atoms:** Moves atoms that are outside the unit cell box back into it (based on fractional coordinates).
7. Click **Generate lattice!**
    * The script will generate the unit cell and symmetry mates.
    * Progress is shown in the status bar.

---
## Visualization & Management

Once the lattice is generated, you can control the visualization:
![Lattice Explorer Overview](images/main_screen_demo.png)

* **Select lattice:** Choose which loaded model to control.
* **Delete Lattice:** Removes all generated lattice and cell objects for the selected model.
* **Set lattice view:** Quickly align the view to the `ab`, `ac`, or `bc` faces of the unit cell.
* **Color by symmetry operator:** Colors each ASU symmetry mate according to its symmetry operator.
* **Show grid points:** Toggles spheres at the unit cell vertices.
* **Grid point size / Cell thickness:** Fine-tune the visual representation of the unit cell boundaries.

### Exporting Views
Click **Export views...** to open the advanced export dialog. The tool remembers your last used directory and DPI settings for a faster workflow.
![Lattice Explorer Overview](images/export_screen.png)

1. **Output folder:** Select where to save your images.
2. **Dimensions & Units:** Set the image size in `px`, `in`, or `cm`.
3. **DPI:** Choose a standard resolution or select **custom...** to input your own.
4. **Background:** Choose white, black, or **custom...** (supports color names or packed RGB hex values like `0xffffff`).
5. **Render faces:** Select which unit cell faces (`ab`, `ac`, `bc`) to batch render.
6. **Distance Width/Height:** Define the physical span in Angstroms for consistent scaling across different structures.
7. **Ray trace:** Toggle high-quality ray tracing for the final output.
8. Click **Render!** to generate the PNG images.

---
## Example Outputs

### 6A27
![6A27](images/example_6a27_600dpi.png)

### 6A28
![6A28](images/example_6a28_600dpi.png)

### 6A29
![6A29](images/example_6a29_600dpi.png)

### 6BDU
![6BDU](images/example_6bdu_600dpi.png)

### 6MC6
![6MC6](images/example_6mc6_600dpi.png)

### 6MC8
![6MC8](images/example_6mc8_600dpi.png)

### 6NEO
![6NEO](images/example_6neo_600dpi.png)

### 6O5L
![6O5L](images/example_6o5l_600dpi.png)

### 9OM8
![9OM8](images/example_9om8_600dpi.png)

### 9OR6
![9OR6](images/example_9or6_600dpi.png)

### 9Y6E
![9Y6E](images/example_9y6e_600dpi.png)

### 9YI3
![9YI3](images/example_9yi3_600dpi.png)

### 9YL4
![9YL4](images/example_9yl4_600dpi.png)

### 9YUP
![9YUP](images/example_9yup_600dpi.png)

---
## Credits

Lattice Explorer was developed to facilitate the analysis of crystal packing interfaces.

[Lattice Explorer](lattice_explorer.py) was written by Robert Szabla of the [Junop Lab](https://junoplab.wordpress.com/) at Western Univerity. 

[PyMOL](https://pymol.org/) is maintained and distributed by Schrödinger.

---
## License
Lattice Explorer is licensed under the GNU General Public License V3.0.
