#####################
# Import Statements #
#####################
from pymol import cmd, util, stored
from string import ascii_uppercase
from tkinter import Tk, Toplevel, Label, Entry, Button, Checkbutton, Radiobutton, StringVar, IntVar, OptionMenu, W, E, DISABLED, NORMAL, Frame
from tkinter.ttk import Progressbar, Separator
from tkinter import filedialog
from tkinter.messagebox import showinfo, askyesno
import math
import itertools
import os
import numpy as np



###################
# Script Defaults #
###################
script_dir = os.path.dirname(os.path.abspath(__file__))
input_dir = os.path.join(script_dir, "input_structures")
output_dir = os.path.join(script_dir, "outputs")

default_unit_cell_expansion = 1
default_wrap_atoms = True
default_color_asu_by_operator = True
default_show_grid = False
default_face = 'ab'

######################
# Rendering Settings #
######################
face = default_face   # either 'ab', 'ac', or 'bc'
w = '1.55in'  # include un as either 'in' or 'cm'
h = '1.55in'  # include un as either 'in' or 'cm'
dpi = 300
min_span_w = 0
min_span_h = 0
framing_mode = True
show_atoms = True
export_png = False
ray_trace=0 # 0: no ray tracing, 1: ray tracing
wrap_atoms = default_wrap_atoms  # wrap atoms into the unit cell when generating the lattice
color_asu_by_operator = default_color_asu_by_operator
show_grid = default_show_grid

ori_color = 'yellow'
unit_cell_edge_color = 'red'
unit_cell_cent_color = 'orange'
unit_cell_corner_size = 4.0
unit_cell_thickness_view = 5  # set thickness of unit cell in viewport
unit_cell_thickness_ray = 2.5 # set thickness of unit cell when ray tracing








# Convert unit cell parameters to vectors in cartesian coordinates
def unit_cell_vectors(unit_cell_params, zero_threshold=1e-10):
    a, b, c, alpha, beta, gamma, _ = unit_cell_params

    # Convert angles to radians
    alpha_rad = math.radians(alpha)
    beta_rad = math.radians(beta)
    gamma_rad = math.radians(gamma)

    # Calculate components
    ax = a
    ay = 0.0
    az = 0.0

    bx = b * math.cos(gamma_rad)
    by = b * math.sin(gamma_rad)
    bz = 0.0

    cx = c * math.cos(beta_rad)
    cy = c * (math.cos(alpha_rad) - math.cos(beta_rad) * math.cos(gamma_rad)) / math.sin(gamma_rad)
    
    cz_square = c**2 - cx**2 - cy**2
    cz = math.sqrt(cz_square) if cz_square > 0 else 0.0

    # Helper to clean small values
    def clean(val):
        return 0.0 if abs(val) < zero_threshold else val

    # Define and clean the vectors
    vec_a = [clean(ax), clean(ay), clean(az)]
    vec_b = [clean(bx), clean(by), clean(bz)]
    vec_c = [clean(cx), clean(cy), clean(cz)]

    return vec_a, vec_b, vec_c



def unit_cell_to_cartesian(unit_cell_params):
    a, b, c, alpha, beta, gamma, _ = unit_cell_params

    # Convert angles from degrees to radians
    alpha_r = math.radians(alpha)
    beta_r = math.radians(beta)
    gamma_r = math.radians(gamma)

    # Compute Cartesian coordinates of the unit cell basis vectors
    ax, ay, az = a, 0.0, 0.0
    bx = b * math.cos(gamma_r)
    by = b * math.sin(gamma_r)
    bz = 0.0

    cx = c * math.cos(beta_r)
    cy = c * (math.cos(alpha_r) - math.cos(beta_r) * math.cos(gamma_r)) / math.sin(gamma_r)
    cz = math.sqrt(c**2 - cx**2 - cy**2)

    # Define origin and 7 other vertices from combinations of vectors a, b, and c
    origin = (0.0, 0.0, 0.0)
    a_vec = (ax, ay, az)
    b_vec = (bx, by, bz)
    c_vec = (cx, cy, cz)

    def add_vectors(*vecs):
        return tuple(sum(coords) for coords in zip(*vecs))

    vertices = [
        origin,
        a_vec,
        b_vec,
        c_vec,
        add_vectors(a_vec, b_vec),
        add_vectors(a_vec, c_vec),
        add_vectors(b_vec, c_vec),
        add_vectors(a_vec, b_vec, c_vec)
    ]

    # Calculate the geometric center as the average of all 8 vertices
    center = tuple(sum(coord[i] for coord in vertices) / 8.0 for i in range(3))

    return vertices + [center]



def show_cell_bounds(model):
    print(f'Generating unit cell for {model}...')
    cell = model+'_cell'
    cmd.copy(cell, model)
    cmd.remove(cell) # remove all atoms from the cell object
    
    # place pseudoatom at the vertices and center of the unit cell
    unit_cell_params = cmd.get_symmetry(model)
    vertices = unit_cell_to_cartesian(unit_cell_params)
    resi = 0
    resns = ['ORI', 'A', 'B', 'C', 'AB', 'AC', 'BC', 'ABC']
    for vertex in vertices[:-1]:
        cmd.pseudoatom(cell, pos=vertex, resi=resi, resn=resns[resi], name='VTX')
        resi+=1
    cmd.pseudoatom(cell, pos=vertices[-1], resi=9, resn='CEN', name='CEN')
    
    # Set cell representation
    cmd.show_as('cell', cell)
    cmd.color(unit_cell_edge_color, cell)
    cmd.color(unit_cell_cent_color, f"/{cell}/PSDO/P/CEN/CEN")
    #cmd.spectrum('resi', selection=cell)
    cmd.set('cgo_line_width', unit_cell_thickness_view, cell) 
    cmd.set('cgo_line_radius', unit_cell_thickness_ray, cell) 
    print('\tDone!')



def find_stretched_bonds(selection='all', state=1, cutoff=5.0,
                         selname='stretched_bond_atoms', do_select=True,
                         do_unbond=False, quiet=0):
    """
    Find (and optionally break) bonds whose current length exceeds `cutoff` Å.
    - selection: atoms/objects to consider (default: all)
    - state: coordinate state to use
    - cutoff: length threshold in Å
    - selname: name for the atom selection containing all atoms in stretched bonds
    - do_select: create/replace `selname` selection
    - do_unbond: if True, call cmd.unbond() for each stretched pair
    """
    m = cmd.get_model(selection, state)                # atoms & bonds in a chempy model
    atoms = m.atom
    if do_select:
        cmd.select(selname, 'none')

    n = 0
    for b in m.bond:                                   # iterate existing bonds
        i, j = b.index                                 # zero-based indices into m.atom
        ai, aj = atoms[i], atoms[j]
        dx = ai.coord[0] - aj.coord[0]
        dy = ai.coord[1] - aj.coord[1]
        dz = ai.coord[2] - aj.coord[2]
        d = math.sqrt(dx*dx + dy*dy + dz*dz)
        if d > cutoff:
            # use (object AND index N) to target the exact atoms
            sel1 = f"({ai.model} and index {ai.index})"
            sel2 = f"({aj.model} and index {aj.index})"
            if do_select:
                cmd.select(selname, f"{sel1} or {sel2}", merge=1)
            if do_unbond:
                cmd.unbond(sel1, sel2)
            if not quiet:
                print(f"> {ai.model}`{ai.index} -- {aj.model}`{aj.index}: {d:.2f} Å")
            n += 1

    if not quiet:
        action = "unbonded" if do_unbond else "found"
        print(f"{action} {n} stretched bonds (> {cutoff:.2f} Å).")



def _cell_vectors(a, b, c, alpha_deg, beta_deg, gamma_deg):
    """
    Build cell vectors (columns of A) in Cartesian space.
    """
    alpha = math.radians(alpha_deg)
    beta  = math.radians(beta_deg)
    gamma = math.radians(gamma_deg)

    ax, ay, az = a, 0.0, 0.0
    bx = b * math.cos(gamma)
    by = b * math.sin(gamma)
    bz = 0.0

    # See e.g. Giacovazzo (Fundamentals of Crystallography) for this construction
    cx = c * math.cos(beta)
    if abs(math.sin(gamma)) < 1e-12:
        raise ValueError("Unit cell gamma too close to 0/180°, cannot build cell.")
    cy = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
    # Ensure cz is positive (right-handed basis)
    cz_sq = c*c - cx*cx - cy*cy
    if cz_sq < 0 and cz_sq > -1e-8:  # numerical noise
        cz_sq = 0.0
    if cz_sq < 0:
        raise ValueError("Computed negative cz^2; check unit cell parameters.")
    cz = math.sqrt(cz_sq)

    # A has columns [a_vec, b_vec, c_vec]
    A = [[ax, bx, cx],
         [ay, by, cy],
         [az, bz, cz]]
    return A

def _mat_inv_3x3(M):
    """
    Invert a 3x3 matrix M (list of lists). Returns Minv.
    """
    (a,b,c), (d,e,f), (g,h,i) = M
    A =  (e*i - f*h)
    B = -(d*i - f*g)
    C =  (d*h - e*g)
    D = -(b*i - c*h)
    E =  (a*i - c*g)
    F = -(a*h - b*g)
    G =  (b*f - c*e)
    H = -(a*f - c*d)
    I =  (a*e - b*d)
    det = a*A + b*B + c*C
    if abs(det) < 1e-18:
        raise ValueError("Singular matrix (det≈0); invalid unit cell.")
    inv_det = 1.0/det
    return [[A*inv_det, D*inv_det, G*inv_det],
            [B*inv_det, E*inv_det, H*inv_det],
            [C*inv_det, F*inv_det, I*inv_det]]

def _mat_vec(M, v):
    return [M[0][0]*v[0] + M[0][1]*v[1] + M[0][2]*v[2],
            M[1][0]*v[0] + M[1][1]*v[1] + M[1][2]*v[2],
            M[2][0]*v[0] + M[2][1]*v[1] + M[2][2]*v[2]]

def _inside_frac(frac_xyz, tol):
    # In [0,1) with tolerance; atoms exactly on 1.0 - tol are considered inside.
    for u in frac_xyz:
        if u < -tol or u >= 1.0 - tol:
            return False
    return True

def _wrap01(x):
    """
    Wrap a float into [0,1). Uses floor; works for negatives too.
    """
    return x - math.floor(x)



def wrap_to_unit_cell(obj_name, selection=None, state=0, quiet=0):
    """
    Bring atoms of a model back into the crystallographic unit cell by
    translating along integer multiples of a/b/c.

    Parameters
    ----------
    obj_name : str
        The PyMOL object name.
    selection : str or None
        Selection string to restrict which atoms to wrap (default: object only).
    state : int
        Model state to operate on. 0 means all states.
    quiet : int
        If non-zero, suppress messages.

    Notes
    -----
    - Uses cmd.get_symmetry(obj) for (a,b,c,alpha,beta,gamma,spacegroup).
    - Works per-atom; atoms exactly on a boundary are left as-is (wrapped with floor).
    - Pure Python math, no NumPy required.
    """
    if selection is None:
        selection = f"model {obj_name}"
    else:
        selection = f"({selection}) and model {obj_name}"

    sym = cmd.get_symmetry(obj_name)
    if not sym or len(sym) < 6:
        raise ValueError("No unit cell found for object. Was this structure solved by crystallography?")
    a, b, c, alpha, beta, gamma = sym[:6]

    # Build cell matrix and its inverse (Cartesian <-> fractional)
    A = _cell_vectors(a, b, c, alpha, beta, gamma)
    Ainv = _mat_inv_3x3(A)

    # Determine which states to operate on
    if state == 0:
        n_states = cmd.count_states(obj_name)
        states = range(1, n_states + 1)
    else:
        states = [state]

    total_moved = 0

    for st in states:
        # Read atom coords & indexes
        m = cmd.get_model(selection, state=st)
        if not m.atom:
            continue

        # Build new coordinate map keyed by atom index
        coord_map = {}
        moved_here = 0

        for at in m.atom:
            r = [at.coord[0], at.coord[1], at.coord[2]]
            # Cartesian -> fractional
            f = _mat_vec(Ainv, r)

            # Wrap into [0,1)
            fw = [_wrap01(f[0]), _wrap01(f[1]), _wrap01(f[2])]

            # If unchanged, skip (small speedup)
            if (abs(fw[0]-f[0]) < 1e-12 and
                abs(fw[1]-f[1]) < 1e-12 and
                abs(fw[2]-f[2]) < 1e-12):
                # Already inside [0,1)^3
                continue

            # fractional -> Cartesian (wrapped)
            r_new = _mat_vec(A, fw)

            coord_map[at.index] = (r_new[0], r_new[1], r_new[2])
            moved_here += 1

        # Apply updates for this state (if any)
        if coord_map:
            # Use alter_state so we can set x,y,z directly; pass our map in 'space'
            cmd.alter_state(
                st,
                selection,
                "x,y,z = coord_map.get(index, (x,y,z))",
                space={'coord_map': coord_map}
            )
            # force a small rebuild to refresh geometry
            #cmd.rebuild(selection, 1)

        total_moved += moved_here
        if not quiet:
            print(f"[wrap_to_unit_cell] State {st}: moved {moved_here} atoms.")

    if not quiet:
        print(f"[wrap_to_unit_cell] Done. Total atoms moved: {total_moved}.")



def generate_unit_cell(model, basename='', wrap_atoms=True, progress_cb=None):
    print(f"Generating a full unit cell for '{model}'...")
    if not basename:
        basename=model
    # generate symmates within at least 1 full unit cell
    sym_prefix = basename + '_sym_'
    identity_sym = sym_prefix + '00000000'
    unit_cell_params = cmd.get_symmetry(model)
    sym_search_cuttoff = unit_cell_params[0] + unit_cell_params[1] + unit_cell_params[2]
    cmd.copy(identity_sym, model)
    cmd.symexp(sym_prefix, model, model, cutoff=sym_search_cuttoff)

    # delete all symmates that are not in the 00 unit cell
    address_values = ['00', '01', '-1']
    all_cells = itertools.product(address_values, repeat=3)
    for cell_address in all_cells:
        cell_address_str = ''.join(cell_address)
        if cell_address_str != '000000':
            cmd.delete(f"{sym_prefix}*{cell_address_str}")
    print("\tDone!")

    # Wrap atoms that lie outside the unit cell back into the unit cell
    if progress_cb:
        progress_cb(0)

    if wrap_atoms:
        symmate_list = cmd.get_object_list(sym_prefix+'*')
        total_symmates = max(len(symmate_list), 1)
        if not symmate_list and progress_cb:
            progress_cb(25)
        for idx, symmate in enumerate(symmate_list, start=1):
            print(f"Wrapping atoms inside unit cell for '{symmate}'...")
            wrap_to_unit_cell(symmate, quiet=1)
            find_stretched_bonds(symmate, cutoff=5.0, do_select=False, do_unbond=True, quiet=1)
            print("\tDone!")
            if progress_cb:
                progress_cb(25 * (idx / total_symmates))
    elif progress_cb:
        progress_cb(25)


def generate_expansion_range(unit_cell_expansion: int) -> list[str]:
    if not (0 <= unit_cell_expansion <= 9):
        raise ValueError("unit_cell_expansion must be an integer between 0 and 9")

    result = []
    for i in range(-unit_cell_expansion, unit_cell_expansion + 1):
        if i < 0:
            result.append(f"{i}")
        else:
            result.append(f"{i:02d}")
    return result


def get_lattice_translation_vectors(model, lattice_size=1):
    # calculate the translation vectors for each unit cell in a crystal lattice
    # The crystal size is defined by the parameter "lattice_size"
    # "lattice_size" is the number of unit cells to expand in each direction from the reference cell 
    # eg: 0 = 1x1x1, 1 = 3x3x3, 2 = 5x5x5, etc...
    print(f"Calculating lattice translation vectors for '{model}'....")
    unit_cell_params = cmd.get_symmetry(model)
    vec_a, vec_b, vec_c = unit_cell_to_cartesian(unit_cell_params)[1:4]
    address_values = generate_expansion_range(lattice_size)
    all_cells = list(itertools.product(address_values, repeat=3))
    translation_vectors = {}
    # move the reference cell address to the front of the list
    all_cells.sort(key=lambda x: x != ('00', '00', '00'))
    for cell_address in all_cells:
        a, b, c = map(int, cell_address)
        cell_address_str = ''.join(cell_address)
        # Scale each vector
        vec_a_scaled = [a * component for component in vec_a]
        vec_b_scaled = [b * component for component in vec_b]
        vec_c_scaled = [c * component for component in vec_c]
        # combine the scaled vectors into a single translation vector
        vec_sum = [x + y + z for x, y, z in zip(vec_a_scaled, vec_b_scaled, vec_c_scaled)]
        translation_vectors[cell_address_str] = vec_sum
    print("\tDone!")
    return translation_vectors



def expand_cell_bounds(base_cell, translation_vectors):
    state=1
    lattice_obj = f"{base_cell}_lattice"
    for cell_address_str, vector in translation_vectors.items():
        # Apply the cell translations - store translated coordinates as a separate object state
        cmd.create(lattice_obj, base_cell, target_state=state)
        cmd.translate(vector, lattice_obj, state=state, camera=0)
        # move to next unit cell/state
        state=state+1



def expand_model_lattice(base_asu, translation_vectors, basename='', split_by_symop=False, split_by_group=False, groupings=None, lattice_size=1, wrap_atoms=True, progress_cb=None):
    # sanity checks
    if split_by_group and not groupings:
        raise ValueError("Must provide a groupings dict if split_by_group=True!")
    if not basename:
        basename=base_asu

    # Generate a full unit cell
    generate_unit_cell(base_asu, basename, wrap_atoms=wrap_atoms, progress_cb=progress_cb)

    print(f"Generating lattice by unit cell translation of '{basename}' within +/-{lattice_size} unit cells...")
    identity_sym = f"{basename}_sym_00000000"
    
    # If not splitting by symmetry operator, merge all symmates into one base_symmate
    symmate_list = cmd.get_object_list(f"{basename}_sym*")
    if not split_by_symop:
        for symmate in symmate_list:
            symdef = symmate.split(f"{basename}_sym_")[1]
            symop = symdef[-8:-6]
            # use segi identifier to tag each symmate with its symmetry operator to avoid overwriting atoms when merging objects
            cmd.alter(symmate, f"segi=segi+'{symop}'" )
        # merge all symmetry mates to one object
        cmd.copy_to(identity_sym, f"{basename}_sym_* and not /{identity_sym}", rename='')
        symmate_list.remove(identity_sym)
        cmd.delete(' '.join(symmate_list))

    # if splitting by custom grouping, split the unit cell into strucutre groups
    symmate_list = cmd.get_object_list(f"{basename}_sym*")
    if split_by_group:
        for symmate in symmate_list:
            for group_name, (group_sele, color) in groupings.items():
                symdef = symmate.split(f"{basename}_sym_")[1]
                atom_selection = f"{symmate} and ({group_sele})"
                split_symmate = f"{basename}_sym_{group_name}_{symdef}"
                cmd.create(split_symmate, atom_selection)
        cmd.delete(' '.join(symmate_list))

    # Applt the unit cell translateion to all symmates
    symmate_list = cmd.get_object_list(f"{basename}_sym*")
    state=1
    total_states = max(len(translation_vectors), 1)
    for state_index, (cell_address_str, vector) in enumerate(translation_vectors.items(), start=1):
        # Apply the cell translations - store translated coordinates as a separate object state
        print(f"\tState #{state:02d}: Translating cell {cell_address_str} to {vector}")
        for symmate in symmate_list:
            symop = symmate[-9:-6]
            if not split_by_symop:
                symop=''
            if split_by_group:
                group = symmate.split(f"{basename}_sym")[1][:-9]
            else:
                group = ''
            lattice_frag = f"{basename}_lattice{group}_model{symop}"
            cmd.create(lattice_frag, symmate, target_state=state)
            cmd.translate(vector, lattice_frag, state=state, camera=0)
        # move to next unit cell/state
        state=state+1
        if progress_cb:
            progress_cb(25 + 50 * (state_index / total_states))
    if progress_cb and not translation_vectors:
        progress_cb(75)
    
    # Cleanup
    cmd.delete(f"{basename}_sym*")
    if split_by_symop and split_by_group:
        for group_name in list(groupings.keys()):
            cmd.group(f"{basename}_lattice_{group_name}_model", f"{basename}_lattice_{group_name}_model_*")            
    return


def generate_lattice(model, lattice_size=1):
    print(f'Generating lattice for {model} within +/-{lattice_size} unit cells...')
    
    # generate symmates within at least 1 full unit cell
    lattice_prefix = model + '_lattice_'
    sym_prefix = model + '_sym_'
    identity_sym = sym_prefix + '00000000'
    unit_cell_params = cmd.get_symmetry(model)
    sym_search_cuttoff = unit_cell_params[0] + unit_cell_params[1] + unit_cell_params[2]
    cmd.symexp(sym_prefix, model, model, cutoff=sym_search_cuttoff)
    cmd.copy(identity_sym, model)

    # delete all symmates that are not in the 00 unit cell
    address_values = ['00', '01', '-1']
    all_cells = itertools.product(address_values, repeat=3)
    for cell_address in all_cells:
        cell_address_str = ''.join(cell_address)
        if cell_address_str != '000000':
            cmd.delete(f"{sym_prefix}*{cell_address_str}")
    
    # expand all symmates beyond the 00 unit cell by translation
    vec_a, vec_b, vec_c = unit_cell_to_cartesian(unit_cell_params)[1:4]
    sym_mate_list = cmd.get_object_list(sym_prefix+'*')
    address_values = generate_expansion_range(lattice_size)
    all_cells = list(itertools.product(address_values, repeat=3))
    # move the reference cell address to the front of the list
    all_cells.sort(key=lambda x: x != ('00', '00', '00'))
    state=1
    for cell_address in all_cells:
        a, b, c = map(int, cell_address)
        cell_address_str = ''.join(cell_address)
        sym_name = f"{sym_prefix}{cell_address_str}"
        # Scale each vector
        vec_a_scaled = [a * component for component in vec_a]
        vec_b_scaled = [b * component for component in vec_b]
        vec_c_scaled = [c * component for component in vec_c]
        # combine the scaled vectors into a single translation vector
        vec_sum = [x + y + z for x, y, z in zip(vec_a_scaled, vec_b_scaled, vec_c_scaled)]
        # Apply the translation to each symmetry operator - store translated coordinates as a separate object state
        print(f"Translating {model} cell {cell_address_str} to {vec_sum}")
        for base_symmate in sym_mate_list:
            symop = base_symmate[-8:-6]
            sym_name = f"{lattice_prefix}{symop}"
            cmd.create(sym_name, base_symmate, target_state=state)
            cmd.translate(vec_sum, sym_name, state=state, camera=0)
        state=state+1
    cmd.delete(sym_prefix+'*')


    print('\tDone!')
    return 



def color_by_operator(model, spectrum='lightblue palecyan palegreen paleyellow lightorange lightpink'):
    print(f'Coloring {model} by symmetry operator...')
    #cmd.spectrum( 'p.op', selection=f'{model}_lattice*' )
    #cmd.spectrum( "int(model[-2:])", selection=f'{model}_lattice*' )
    selection = f'{model}_lattice*_model*'
    cmd.spectrum( "int(model[-2:])", selection=selection, palette=spectrum )
    print('\tDone!')


def organize_groups(model):
    print(f'Organizing groups for {model}...')
    group_name = f'{model}_lattice'
    members = f'{model}_lattice_* + {model}_cell_lattice'
    cmd.group(group_name, members)
    group_name = f'{model}_group'
    members = f'{model} {model}_cell {model}_lattice'
    cmd.group(group_name, members)
    print('\tDone!')



def get_rotation_matrix(vec1, vec2, vec3):
    def compute_axes(z_axis):
        z_axis = z_axis / np.linalg.norm(z_axis)

        # Project vec1 onto plane orthogonal to z_axis
        v1_proj = vec1 - np.dot(vec1, z_axis) * z_axis
        x_axis = v1_proj / np.linalg.norm(v1_proj)

        # Ensure x_axis points in the direction of vec1
        if np.dot(x_axis, vec1) < 0:
            x_axis = -x_axis

        # Compute y_axis as right-handed system
        y_axis = np.cross(z_axis, x_axis)
        y_axis /= np.linalg.norm(y_axis)

        # Recompute x_axis for strict orthogonality
        x_axis = np.cross(y_axis, z_axis)
        x_axis /= np.linalg.norm(x_axis)

        return x_axis, y_axis, z_axis

    # Normalize input vectors
    vec1 = np.array(vec1, dtype=float)
    vec2 = np.array(vec2, dtype=float)
    vec3 = np.array(vec3, dtype=float)

    # Try both directions of z_axis
    z_options = [vec3 / np.linalg.norm(vec3), -vec3 / np.linalg.norm(vec3)]
    best_R = None
    best_score = -np.inf

    for z_axis in z_options:
        x_axis, y_axis, z_axis = compute_axes(z_axis)

        # Project vec2 onto the plane orthogonal to z
        v2_proj = vec2 - np.dot(vec2, z_axis) * z_axis
        v2_proj_norm = np.linalg.norm(v2_proj)

        # If vec2 is nearly parallel to z (no projection), skip
        if v2_proj_norm < 1e-6:
            continue

        v2_proj /= v2_proj_norm

        # Score based on alignment of y_axis with projected vec2
        score = np.dot(y_axis, v2_proj)

        if score > best_score:
            best_score = score
            best_R = np.stack([x_axis, y_axis, z_axis], axis=1)

    if best_R is None:
        raise ValueError("Could not compute a valid rotation matrix with given vectors.")

    return best_R.flatten().tolist()



def set_views(view_name:str='default', rotation:list=None, camera_xyz:list=None, origin:list=None, clip:list=None, orthoscopic:int=None, xtal:str=None, face:str=None, min_span_w=0, min_span_h=0):
    """
    DESCRIPTION:
        Move camera to a specific view. 

    ARGUMENTS:
        `view`:        (required `str`)
                       Name of a pre-defined view from the `camera_view_library`. 
        `rotation`:    (optional `list`)
                       A 3x3 rotaion matrix to position the camera. 
                         If provided, overrides the rotation component of the given view
        `camera_xyz`:  (optional `list`)
                       xyz coordinates of the camera translation (relative to origin of rotation)
                         If provided, overrides the camer_xyz component of the given view
        `origin`:      (optional `list`)
                       xyz coordinates in model space of the origin of rotation. 
                         If provided, overrides the origin component of the given view
        `clip`:        (optional `list`)
                       front and rear clipping distances. 
                         If provided, overrides the clip component of the given view
        `orthoscopic`: (optional `int`)
                        20: orthoscopic view
                       -20: projection view
                       If provided, overrides the orthoscopic setting of the given view
        """
    #print(f'Setting camera view to {view_name}...')
    camera_view_library = {
        'default' : [
                    1.0,  0.0,  0.0,    # rotation   
                    0.0,  1.0,  0.0,    # rotation   
                    0.0,  0.0,  1.0,    # rotation   
                    0.0,  0.0, -1316, #-875,  # camera_xyz 
                    0.0,  0.0,  0.0,    # origin     
                    -1820.0, 2540.0,     # clip        
                    20  ],              # orthoscopic            
    }
    if view_name in camera_view_library.keys():
        camera_view = camera_view_library[view_name]
    else:
        print(f"'{view_name}' is not a valid view! Available views:")
        for available_view in camera_view_library.keys(): print(available_view)
    if rotation:
        camera_view[0:9]   = rotation
    if camera_xyz:
        camera_view[9:12]  = camera_xyz
    if origin:
        camera_view[12:15] = origin
    if clip:
        camera_view[15:17] = clip
    if orthoscopic:
        camera_view[17]    = orthoscopic
    
    if xtal:
        ori_selector = f'({xtal}*_cell) and resn CEN'
        stored.ori = []
        cmd.iterate_state(1, ori_selector, 'stored.ori.append((x,y,z))')
        camera_view[12:15] = stored.ori[0]
    if face:
        unit_cell_params = cmd.get_symmetry(f'({xtal}*_cell)')
        vec_a, vec_b, vec_c = unit_cell_to_cartesian(unit_cell_params)[1:4]
        vector_order = ['a', 'b', 'c']
        vector_reorder = []
        vector_lookup = { 'a': vec_a, 'b': vec_b, 'c': vec_c}

        for vec in face:
            if vec not in vector_order:
                raise ValueError(f"'{vec}' is not a valid face! Possible values: a, b, c")
            vector_order.remove(vec)
            vector_reorder.append(vec)
        vector_reorder.append(vector_order[0])
        rot = get_rotation_matrix(vector_lookup[vector_reorder[0]], vector_lookup[vector_reorder[1]], vector_lookup[vector_order[0]])
        camera_view[0:9] = rot

    cmd.set_view ( camera_view )

    if min_span_w or min_span_h:
        w,h = cmd.get_viewport()
        min_span = max(min_span_w, min_span_h)
        if min_span_w and (h < w):
            min_span = min_span*h/w
        if min_span_h and (w < h):
            min_span = min_span*h/w
        cmd.pseudoatom('min_span_cent', pos=camera_view[12:15])
        cmd.zoom('min_span_cent', min_span/2)     
        cmd.delete('min_span_cent')


def get_view_span():
    view = cmd.get_view()
    vp_w, vp_h = cmd.get_viewport()
    if vp_h <= 0:
        raise ValueError("Viewport height must be greater than zero.")
    fov = abs(view[17])
    cam_dist = abs(view[11])
    if fov <= 0:
        raise ValueError("Field of view must be greater than zero.")
    span_h = 2.0 * cam_dist * math.tan(math.radians(fov / 2.0))
    span_w = span_h * (vp_w / vp_h)
    return span_w, span_h


def find_smallest_face(model):
    unit_cell_params = cmd.get_symmetry(model)
    a, b, c, alpha, beta, gamma, _ = unit_cell_params
    face_areas = {
        'ab': a * b * math.sin(math.radians(gamma)),
        'ac': a * c * math.sin(math.radians(beta)),
        'bc': b * c * math.sin(math.radians(alpha))
    }
    smallest_face = min(face_areas, key=face_areas.get)
    return smallest_face


def _length_to_inches(value):
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip().lower()
    if text.endswith('cm'):
        return float(text[:-2]) / 2.54
    if text.endswith('in'):
        return float(text[:-2])
    return float(text)


def _parse_length_with_unit(value):
    if isinstance(value, (int, float)):
        return float(value), 'px'
    text = str(value).strip().lower()
    for unit in ('cm', 'in', 'px'):
        if text.endswith(unit):
            return float(text[:-len(unit)]), unit
    return float(text), 'px'


def _normalize_units(units):
    text = str(units or '').strip().lower()
    if text in ('px', 'pixel', 'pixels'):
        return 'px'
    if text in ('in', 'inch', 'inches'):
        return 'in'
    if text in ('cm', 'centimeter', 'centimeters'):
        return 'cm'
    return text


def _output_pixels(width_val, height_val, units, dpi_value):
    units = _normalize_units(units)
    if units == 'px':
        return float(width_val), float(height_val)
    if units == 'cm':
        scale = dpi_value / 2.54
        return float(width_val) * scale, float(height_val) * scale
    return float(width_val) * dpi_value, float(height_val) * dpi_value


def _format_png_size(value, units):
    units = _normalize_units(units)
    if units == 'px':
        return int(round(float(value)))
    return f"{value}{units}"


def render_views(
    model:str,
    ray=0,
    atoms=True,
    cell=True,
    output_dir_override=None,
    progress_cb=None,
    image_width=None,
    image_height=None,
    image_units=None,
    dpi_override=None,
    min_span_w_override=None,
    min_span_h_override=None,
    faces_to_render=None,
):
        print(f'Rendering views for {model}...')
        target_dir = output_dir_override or output_dir
        if target_dir and not os.path.exists(target_dir):
            os.makedirs(target_dir, exist_ok=True)
    
        if image_width is None or image_height is None:
            width_val, width_unit = _parse_length_with_unit(w)
            height_val, height_unit = _parse_length_with_unit(h)
            if width_unit != height_unit:
                height_unit = width_unit
            image_units = _normalize_units(image_units or width_unit)
        else:
            width_val = float(image_width)
            height_val = float(image_height)
            image_units = _normalize_units(image_units or 'px')
    
        dpi_value = float(dpi_override) if dpi_override is not None else dpi
        out_w_px, out_h_px = _output_pixels(width_val, height_val, image_units, dpi_value)
        if out_w_px <= 0 or out_h_px <= 0:
            raise ValueError("Output size must be greater than zero.")
        desired_ratio = out_w_px / out_h_px
        png_width = _format_png_size(width_val, image_units)
        png_height = _format_png_size(height_val, image_units)
    
        min_span_w_value = min_span_w if min_span_w_override is None else float(min_span_w_override)
        min_span_h_value = min_span_h if min_span_h_override is None else float(min_span_h_override)
    
        def match_viewport_aspect():
            vp_w, vp_h = cmd.get_viewport()
            if vp_h <= 0:
                return
            current_ratio = vp_w / vp_h
            if abs(current_ratio - desired_ratio) > 1e-3:
                new_width = max(1, int(round(vp_h * desired_ratio)))
                cmd.viewport(new_width, vp_h)
        
        faces = faces_to_render if faces_to_render is not None else ['ab', 'ac', 'bc']
        total_faces = len(faces)
        for idx, face in enumerate(faces, start=1):
            if progress_cb:
                progress_cb(face, idx, total_faces, started=True)
            print(f'\tRendering {face} images for {model}...')
            atoms_filepath = os.path.join(target_dir, f'{model}_{face}_atoms.png')
            cell_filepath  = os.path.join(target_dir, f'{model}_{face}_cell.png')
            
            if atoms:
                cmd.disable(f'all')
                cmd.enable(f'{model}_*lattice*',parents=1)
                match_viewport_aspect()
                clip = cmd.get_view()[15:17]
                set_views(xtal=model, face=face, clip=clip, min_span_w=min_span_w_value, min_span_h=min_span_h_value)
                cmd.center(f'{model}_cell')
                cmd.png(atoms_filepath, width=png_width, height=png_height, dpi=dpi_value, ray=ray)
            
            if cell:
                cmd.disable(f'all')
                cmd.enable(f'{model}',parents=1)
                cmd.enable(f'{model}_cell', parents=1)
                match_viewport_aspect()
                clip = cmd.get_view()[15:17]
                set_views(xtal=model, face=face, clip=clip, min_span_w=min_span_w_value, min_span_h=min_span_h_value)
                cmd.center(f'{model}_cell')
                cmd.png(cell_filepath, width=png_width, height=png_height, dpi=dpi_value, ray=1)
                #cmd.disable(f'{model}_cell')
                #cmd.enable(f'{model}')
                #cmd.png(asu_filepath, width=w, height=h, dpi=dpi, ray=1)
            if progress_cb:
                progress_cb(face, idx, total_faces, started=False)
        print('Done!')

    

def global_settings():
    print('Applying global settings...')
    cmd.set('orthoscopic', 1)
    cmd.set('all_states', 1)
    cmd.set('sphere_scale', unit_cell_corner_size)
    print('\tDone!')










def lattice_exists(model):
    lattice_name = f"{model}_lattice"
    try:
        object_names = set(cmd.get_names('objects') or [])
        group_names = set(cmd.get_names('groups') or [])
        if lattice_name not in object_names and lattice_name not in group_names:
            return False
        return bool(cmd.get_object_list(lattice_name))
    except Exception:
        return False


def get_existing_lattice_models():
    names = set()
    try:
        names.update(cmd.get_names('objects') or [])
    except Exception:
        pass
    try:
        names.update(cmd.get_names('groups') or [])
    except Exception:
        pass

    models = set()
    for name in names:
        if name.endswith('_lattice'):
            models.add(name[:-len('_lattice')])
        elif name.endswith('_cell_lattice'):
            models.add(name[:-len('_cell_lattice')])
        elif '_lattice_model' in name:
            models.add(name.split('_lattice_model')[0])
    models = {model for model in models if not model.endswith('_cell')}
    return sorted(models)


def generate_lattice_for_model(model, lattice_size, wrap_atoms, progress_cb=None):
    if progress_cb is None:
        progress_cb = lambda _value: None
    global_settings()
    if not cmd.get_symmetry(model):
        raise ValueError("No unit cell found for object. Is this a crystal structure?")

    progress_cb(0)
    cmd.remove(f'{model} and elem H')
    show_cell_bounds(model)
    translation_vectors = get_lattice_translation_vectors(model=model, lattice_size=lattice_size)

    expand_cell_bounds(base_cell=f"{model}_cell", translation_vectors=translation_vectors)
    if show_atoms:
        expand_model_lattice(
            base_asu=model,
            translation_vectors=translation_vectors,
            basename=model,
            split_by_symop=True,
            split_by_group=False,
            lattice_size=lattice_size,
            wrap_atoms=wrap_atoms,
            progress_cb=progress_cb,
        )

    progress_cb(75)
    final_steps = 5
    final_step = 0
    def bump_final():
        nonlocal final_step
        final_step += 1
        progress_cb(75 + (final_step / final_steps) * 25)

    cmd.color(color=unit_cell_edge_color, selection=f'/{model}_cell*////VTX')
    cmd.color(color=unit_cell_cent_color, selection=f'/{model}_cell*////CEN')
    cmd.color(color=ori_color, selection=f"/{model}_cell///ORI")
    cmd.set('sphere_scale', unit_cell_corner_size * 1.5, f"/{model}_cell///ORI")
    bump_final()
    util.cbc(model)
    cmd.show_as('cartoon', model)
    cmd.show_as('lines', f'{model}_lattice*')
    bump_final()

    organize_groups(model)
    bump_final()

    cmd.orient(f"{model}_group")
    camera_xyz = cmd.get_view()[9:12]
    clip = cmd.get_view()[15:17]
    set_views(xtal=model, face=default_face, camera_xyz=camera_xyz, clip=clip)
    cmd.center(f'{model}_cell')
    #cmd.viewport(700, 700)
    bump_final()


def launch_lattice_viewer():
    window = Tk()
    window.title('PyMOL Lattice Viewer')
    try:
        window.deiconify()
        window.lift()
        window.attributes('-topmost', True)
        window.after(200, lambda: window.attributes('-topmost', False))
    except Exception as exc:
        print(f"[pymol_lattice_viewer] window show hint failed: {exc}")

    input_path_var = StringVar()
    pdb_code_var = StringVar()
    input_mode_var = StringVar(value='file')
    unit_cell_var = StringVar(value=str(default_unit_cell_expansion))
    wrap_atoms_var = IntVar(value=1 if default_wrap_atoms else 0)
    show_model_var = IntVar(value=1)
    color_asu_var = IntVar(value=1 if default_color_asu_by_operator else 0)
    color_scheme_var = StringVar(value='Pastel')
    show_grid_var = IntVar(value=1 if default_show_grid else 0)
    grid_point_size_var = StringVar(value=str(unit_cell_corner_size))
    cell_thickness_var = StringVar(value=str(unit_cell_thickness_view))
    unit_cell_info_var = StringVar(value="No model selected")
    selected_model_var = StringVar()
    view_ab = None
    view_ac = None
    view_bc = None
    export_button = None
    delete_button = None
    scheme_menu = None
    export_window_ref = {'window': None}
    model_states = {}
    current_model = {'name': ''}
    
    last_export_settings = {
        'output_dir': output_dir,
        'dpi': dpi
    }
    color_scheme_map = {
        'Pastel': 'lightblue palecyan palegreen paleyellow lightorange lightpink',
        'Rainbow': 'rainbow',
        'Warm': 'red orange yellow lightorange lightpink magenta',
        'Cool': 'blue cyan lightblue palecyan',
    }

    def browse_input_file():
        filepath = filedialog.askopenfilename(
            initialdir="/",
            title="Select a File",
            filetypes=(("Atomic coordinates", ".pdb .cif"), ("All files", "*.*")),
        )
        if filepath:
            input_mode_var.set('file')
            input_path_var.set(filepath)
            pdb_code_var.set('')
            update_input_mode()
            input_entry.focus_set()

    def only_numbers(char):
        return char.isdigit()

    def update_view_buttons():
        if view_ab is None or view_ac is None or view_bc is None:
            return
        state = NORMAL if selected_model_var.get().strip() else DISABLED
        view_ab.configure(state=state)
        view_ac.configure(state=state)
        view_bc.configure(state=state)
        if export_button is not None:
            export_button.configure(state=state)
        if delete_button is not None:
            delete_button.configure(state=state)

    def save_model_state(model):
        if not model:
            return
        model_states[model] = {
            'color_asu': bool(color_asu_var.get()),
            'show_grid': bool(show_grid_var.get()),
            'wrap_atoms': bool(wrap_atoms_var.get()),
            'show_model': bool(show_model_var.get()),
            'color_scheme': color_scheme_var.get(),
        }

    def apply_display_options(*_args):
        model = selected_model_var.get().strip()
        if not model:
            return

        try:
            if show_model_var.get():
                cmd.enable(f"{model}_group")
            else:
                cmd.disable(f"{model}_group")
        except Exception as e:
            print(f"[pymol_lattice_viewer] Error setting model visibility: {e}")
        
        try:
            if color_asu_var.get():
                palette = color_scheme_map.get(color_scheme_var.get(), color_scheme_map['Pastel'])
                color_by_operator(model, spectrum=palette)
                if scheme_menu:
                    scheme_menu.configure(state=NORMAL)
            else:
                cmd.color(color='gray50', selection=f'{model}_lattice*_model*')
                if scheme_menu:
                    scheme_menu.configure(state=DISABLED)
        except Exception as e:
            print(f"[pymol_lattice_viewer] Error setting colors: {e}")

        try:
            cmd.color(color=unit_cell_edge_color, selection=f'/{model}_cell*////VTX')
            cmd.color(color=unit_cell_cent_color, selection=f'/{model}_cell*////CEN')
            cmd.color(color=ori_color, selection=f"/{model}_cell///ORI")

            # Apply grid point size
            try:
                g_size = float(grid_point_size_var.get())
                cmd.set('sphere_scale', g_size * 1.5, f"/{model}_cell///ORI")
                cmd.set('sphere_scale', g_size)
            except ValueError:
                pass

            # Apply cell thickness
            try:
                c_thick = float(cell_thickness_var.get())
                cmd.set('cgo_line_width', c_thick, f"{model}_cell")
                cmd.set('cgo_line_radius', c_thick / 2.0, f"{model}_cell")
            except ValueError:
                pass

            if show_grid_var.get():
                cmd.show('spheres', f"{model}_cell*")
            else:
                cmd.hide('spheres', f"{model}_cell*")
        except Exception as e:
            print(f"[pymol_lattice_viewer] Error setting grid: {e}")

        save_model_state(model)

    def delete_lattice_objects(model):
        cmd.delete(f"{model}_cell*")
        cmd.delete(f"{model}_lattice*")
        cmd.delete(f"{model}_cell_lattice")
        cmd.delete(f"{model}_group")
        cmd.delete(f"{model}_lattice")

    def on_model_change(*_args):
        new_model = selected_model_var.get().strip()
        
        if new_model:
            try:
                sym = cmd.get_symmetry(new_model)
                if sym:
                    a, b, c, alpha, beta, gamma, sg = sym
                    
                    # Count symmetry operators based on generated lattice objects
                    # Objects are named like {model}_lattice_model{symop}
                    sym_op_count = len(cmd.get_object_list(f"{new_model}_lattice_model*"))
                    
                    info_text = (
                        f"Space Group: {sg}\n"
                        f"Unit Cell: {a:.2f}, {b:.2f}, {c:.2f}\n"
                        f"Angles: {alpha:.2f}, {beta:.2f}, {gamma:.2f}\n"
                        f"Symmetry Operators: {sym_op_count}"
                    )
                    unit_cell_info_var.set(info_text)
                else:
                    unit_cell_info_var.set("No symmetry data found")
            except Exception:
                unit_cell_info_var.set("Error retrieving symmetry")
        else:
            unit_cell_info_var.set("No model selected")

        if current_model['name'] == new_model:
            update_view_buttons()
            apply_display_options()
            return
        save_model_state(current_model['name'])
        current_model['name'] = new_model
        state = model_states.get(new_model)
        if state:
            color_asu_var.set(1 if state.get('color_asu', True) else 0)
            show_grid_var.set(1 if state.get('show_grid', False) else 0)
            wrap_atoms_var.set(1 if state.get('wrap_atoms', True) else 0)
            show_model_var.set(1 if state.get('show_model', True) else 0)
            if state.get('color_scheme') in color_scheme_map:
                color_scheme_var.set(state.get('color_scheme'))
        update_view_buttons()
        apply_display_options()

    def refresh_model_menu():
        models = get_existing_lattice_models()
        menu = model_menu['menu']
        menu.delete(0, 'end')
        for model in models:
            menu.add_command(label=model, command=lambda value=model: selected_model_var.set(value))
        if models:
            if selected_model_var.get() not in models:
                selected_model_var.set(models[0])
        else:
            selected_model_var.set('')
        update_view_buttons()

    def run_button_press():
        input_path = input_path_var.get().strip()
        pdb_code = pdb_code_var.get().strip()
        lattice_text = unit_cell_var.get().strip()
        input_mode = input_mode_var.get()

        if input_mode == 'file':
            if not input_path:
                showinfo(title='Error!', message='Please select an input file.')
                return
            if not os.path.exists(input_path):
                showinfo(title='Error!', message='Input file not found. Please select a valid file.')
                return
        elif input_mode == 'pdb':
            if not pdb_code:
                showinfo(title='Error!', message='Please enter a PDB code.')
                return
            if not pdb_code.isalnum():
                showinfo(title='Error!', message='PDB code must be alphanumeric.')
                return
        else:
            showinfo(title='Error!', message='Please select an input mode.')
            return
        if not lattice_text.isdigit():
            showinfo(title='Error!', message='Unit cell expansion must be an integer.')
            return
        lattice_size = int(lattice_text)
        if not (0 <= lattice_size <= 9):
            showinfo(title='Error!', message='Unit cell expansion must be between 0 and 9.')
            return

        model_name = os.path.splitext(os.path.basename(input_path))[0] if input_mode == 'file' else pdb_code.lower()

        object_names = set(cmd.get_names('objects') or [])
        lattice_present = lattice_exists(model_name)
        model_present = model_name in object_names
        rebuild_existing = False
        if lattice_present or model_present:
            wrap_text = 'Yes' if wrap_atoms_var.get() else 'No'
            message = (
                f'Lattice for "{model_name}" already exists. Rebuild lattice with specified parameters:\n'
                f'Unit Cell expansion: {lattice_size}\n'
                f'Wrap Atoms in Unit Cell: {wrap_text}'
            )
            if not askyesno(title='Rebuild lattice?', message=message):
                return
            rebuild_existing = True
            delete_lattice_objects(model_name)
            if model_present:
                cmd.delete(model_name)

        needs_load = rebuild_existing or model_name not in set(cmd.get_names('objects') or [])
        if needs_load:
            try:
                if input_mode == 'file':
                    cmd.load(input_path, model_name)
                else:
                    cmd.fetch(model_name, model_name)
            except Exception as exc:
                showinfo(title='Error!', message=f'Failed to load model: {exc}')
                return

        progress_bar['value'] = 0
        status_label.configure(text='Working...')
        window.update_idletasks()
        try:
            generate_lattice_for_model(
                model=model_name,
                lattice_size=lattice_size,
                wrap_atoms=bool(wrap_atoms_var.get()),
                progress_cb=update_progress,
            )
        except Exception as exc:
            showinfo(title='Error!', message=f'Failed to generate lattice: {exc}')
            progress_bar['value'] = 0
            status_label.configure(text='Ready!')
            window.update_idletasks()
            return

        refresh_model_menu()
        save_model_state(model_name)
        selected_model_var.set(model_name)
        apply_display_options()
        progress_bar['value'] = 0
        status_label.configure(text='Ready!')
        window.update_idletasks()

    def delete_button_press():
        model = selected_model_var.get().strip()
        if not model:
            return
        if askyesno(title='Delete Lattice', message=f'Are you sure you want to delete the lattice for "{model}"?'):
            cmd.delete(f"{model}_group")
            refresh_model_menu()

    def set_lattice_view(face_value):
        model = selected_model_var.get().strip()
        if not model:
            showinfo(title='Error!', message='No lattice model selected.')
            return
        try:
            camera_xyz = cmd.get_view()[9:12]
            clip = cmd.get_view()[15:17]
            set_views(xtal=model, face=face_value, camera_xyz=camera_xyz, clip=clip)
            cmd.center(f'{model}_cell')
        except Exception as exc:
            showinfo(title='Error!', message=f'Failed to set lattice view: {exc}')

    def open_export_window():
        if export_window_ref['window'] and export_window_ref['window'].winfo_exists():
            export_window_ref['window'].lift()
            export_window_ref['window'].focus_force()
            return

        export_window = Toplevel(window)
        export_window.title('Export Views')
        export_window_ref['window'] = export_window
        models = get_existing_lattice_models()
        model_var = StringVar()
        if selected_model_var.get().strip() in models:
            model_var.set(selected_model_var.get().strip())
        elif models:
            model_var.set(models[0])
        else:
            model_var.set('')

        default_width, default_unit = _parse_length_with_unit(w)
        default_height, default_height_unit = _parse_length_with_unit(h)
        if default_height_unit != default_unit:
            default_height_unit = default_unit

        current_output_dir = last_export_settings.get('output_dir', output_dir)
        output_dir_var = StringVar(value=current_output_dir)
        image_units_var = StringVar(value=default_unit)
        
        current_dpi = last_export_settings.get('dpi', dpi)
        standard_dpis = ['96', '150', '300', '600', '1200']
        
        # Check if current_dpi matches a standard value
        is_standard = False
        try:
            if str(int(current_dpi)) in standard_dpis:
                is_standard = True
        except ValueError:
            pass

        if is_standard:
            dpi_var = StringVar(value=str(int(current_dpi)))
            custom_dpi_var = StringVar(value='')
        else:
            dpi_var = StringVar(value='custom...')
            custom_dpi_var = StringVar(value=str(current_dpi))

        viewport_w, viewport_h = cmd.get_viewport()
        units_default = _normalize_units(default_unit)
        
        # Use a temporary DPI value for calculation if custom is invalid/empty
        try:
            calc_dpi = float(custom_dpi_var.get()) if dpi_var.get() == 'custom...' else float(dpi_var.get())
        except ValueError:
            calc_dpi = 300.0

        if units_default == 'px':
            width_default = viewport_w
            height_default = viewport_h
        elif units_default == 'cm':
            width_default = (viewport_w / calc_dpi) * 2.54
            height_default = (viewport_h / calc_dpi) * 2.54
        else:
            width_default = viewport_w / calc_dpi
            height_default = viewport_h / calc_dpi

        if units_default == 'px':
            image_width_var = StringVar(value=str(int(round(width_default))))
            image_height_var = StringVar(value=str(int(round(height_default))))
        else:
            image_width_var = StringVar(value=f"{width_default:.3f}")
            image_height_var = StringVar(value=f"{height_default:.3f}")
        width_span_var = StringVar()
        height_span_var = StringVar()
        ray_trace_var = IntVar(value=0)
        
        face_ab_var = IntVar(value=1)
        face_ac_var = IntVar(value=1)
        face_bc_var = IntVar(value=1)
        
        current_bg = cmd.get('bg_rgb')

        def to_float_tuple(val):
            # Handle string representation if necessary
            if isinstance(val, str):
                import ast
                try:
                    val = ast.literal_eval(val)
                except (ValueError, SyntaxError):
                    pass

            if isinstance(val, (int, float)):
                # If it's a number, assume it's a packed 24-bit RGB integer
                ival = int(val)
                return (
                    ((ival >> 16) & 0xFF) / 255.0,
                    ((ival >> 8) & 0xFF) / 255.0,
                    (ival & 0xFF) / 255.0
                )
            
            if isinstance(val, (list, tuple)) and len(val) >= 3:
                try:
                    return tuple(float(v) for v in val[:3])
                except (ValueError, TypeError):
                    pass
            return None

        bg_tuple = to_float_tuple(current_bg)

        # Check if current background is close to white (1,1,1) or black (0,0,0)
        def is_color_match(c1, c2, tol=0.01):
            if c1 is None or c2 is None: return False
            if len(c1) != len(c2): return False
            return sum(abs(x-y) for x,y in zip(c1, c2)) < tol

        if is_color_match(bg_tuple, (1.0, 1.0, 1.0)):
            bg_color_var = StringVar(value='white')
            custom_bg_var = StringVar(value='')
        elif is_color_match(bg_tuple, (0.0, 0.0, 0.0)):
            bg_color_var = StringVar(value='black')
            custom_bg_var = StringVar(value='')
        else:
            bg_color_var = StringVar(value='custom...')
            # Use original value for display if possible, else the tuple
            display_val = str(current_bg) if bg_tuple is None else str(list(bg_tuple))
            custom_bg_var = StringVar(value=display_val)

        progress_text = StringVar(value='Ready.')
        span_update_lock = {'busy': False}

        def on_bg_change(*args):
            if bg_color_var.get() == 'custom...':
                custom_bg_entry.configure(state=NORMAL)
            else:
                custom_bg_entry.configure(state=DISABLED)

        bg_color_var.trace_add('write', on_bg_change)

        def on_dpi_change(*args):
            if dpi_var.get() == 'custom...':
                custom_dpi_entry.configure(state=NORMAL)
            else:
                custom_dpi_entry.configure(state=DISABLED)
        
        dpi_var.trace_add('write', on_dpi_change)

        def browse_output_dir():
            folderpath = filedialog.askdirectory(initialdir=output_dir or "/")
            if folderpath:
                output_dir_var.set(folderpath)

        def get_output_ratio():
            try:
                width_val = float(image_width_var.get())
                height_val = float(image_height_var.get())
            except ValueError:
                return None
            if width_val <= 0 or height_val <= 0:
                return None
            return width_val / height_val

        def update_span_from_image_size(*_args):
            if span_update_lock['busy']:
                return
            ratio = get_output_ratio()
            if ratio is None:
                return
            try:
                height_span = float(height_span_var.get())
            except ValueError:
                return
            span_update_lock['busy'] = True
            width_span_var.set(f"{height_span * ratio:.3f}")
            span_update_lock['busy'] = False

        def update_height_span_from_width(*_args):
            if span_update_lock['busy']:
                return
            ratio = get_output_ratio()
            if ratio is None:
                return
            try:
                width_span = float(width_span_var.get())
            except ValueError:
                return
            span_update_lock['busy'] = True
            height_span_var.set(f"{width_span / ratio:.3f}")
            span_update_lock['busy'] = False

        def update_width_span_from_height(*_args):
            if span_update_lock['busy']:
                return
            ratio = get_output_ratio()
            if ratio is None:
                return
            try:
                height_span = float(height_span_var.get())
            except ValueError:
                return
            span_update_lock['busy'] = True
            width_span_var.set(f"{height_span * ratio:.3f}")
            span_update_lock['busy'] = False

        try:
            span_w_default, span_h_default = get_view_span()
        except Exception:
            span_w_default, span_h_default = 0.0, 0.0

        ratio_default = get_output_ratio()
        if ratio_default and span_h_default:
            current_ratio = span_w_default / span_h_default if span_h_default else ratio_default
            if abs(current_ratio - ratio_default) > 1e-3:
                span_w_default = span_h_default * ratio_default

        width_span_var.set(f"{span_w_default:.3f}")
        height_span_var.set(f"{span_h_default:.3f}")

        def update_export_progress(face, idx, total, started):
            if started:
                progress_text.set(f"Rendering {face} images...")
            export_progress['maximum'] = total
            export_progress['value'] = idx if not started else max(idx - 1, 0)
            export_window.update_idletasks()

        def apply_bg_color():
            selection = bg_color_var.get()
            if selection == 'white':
                cmd.set('bg_rgb', [1.0, 1.0, 1.0])
            elif selection == 'black':
                cmd.set('bg_rgb', [0.0, 0.0, 0.0])
            elif selection == 'custom...':
                custom_val = custom_bg_var.get().strip()
                import ast
                try:
                    # Try to parse as list/tuple
                    bg_val = ast.literal_eval(custom_val)
                except Exception:
                    # Fallback: pass as string (e.g. 'grey')
                    bg_val = custom_val
                try:
                    cmd.set('bg_rgb', bg_val)
                except Exception as e:
                    print(f"Error setting custom background: {e}")

        def render_button_press():
            model = model_var.get().strip()
            if not model:
                showinfo(title='Error!', message='Please select a model to render.')
                return
            out_dir = output_dir_var.get().strip()
            if not out_dir:
                showinfo(title='Error!', message='Please select an output folder.')
                return
            
            # Save output directory
            last_export_settings['output_dir'] = out_dir

            try:
                width_val = float(image_width_var.get())
                height_val = float(image_height_var.get())
            except ValueError:
                showinfo(title='Error!', message='Image width and height must be numbers.')
                return
            units_val = _normalize_units(image_units_var.get())
            if units_val not in ('in', 'cm', 'px'):
                showinfo(title='Error!', message='Image units must be in, cm, or px.')
                return
            
            # DPI handling
            dpi_selection = dpi_var.get()
            if dpi_selection == 'custom...':
                try:
                    dpi_value = float(custom_dpi_var.get())
                except ValueError:
                    showinfo(title='Error!', message='Custom DPI must be a number.')
                    return
            else:
                try:
                    dpi_value = float(dpi_selection)
                except ValueError:
                    showinfo(title='Error!', message='DPI must be a number.')
                    return
            
            # Save DPI
            last_export_settings['dpi'] = dpi_value

            if width_val <= 0 or height_val <= 0:
                showinfo(title='Error!', message='Image width and height must be greater than zero.')
                return
            if dpi_value <= 0:
                showinfo(title='Error!', message='DPI must be greater than zero.')
                return
            try:
                span_w = float(width_span_var.get())
                span_h = float(height_span_var.get())
            except ValueError:
                showinfo(title='Error!', message='Span values must be numbers.')
                return
            
            selected_faces = []
            if face_ab_var.get(): selected_faces.append('ab')
            if face_ac_var.get(): selected_faces.append('ac')
            if face_bc_var.get(): selected_faces.append('bc')
            
            if not selected_faces:
                showinfo(title='Error!', message='Please select at least one face to render.')
                return

            try:
                apply_bg_color()
                export_progress['value'] = 0
                progress_text.set('Starting render...')
                export_window.update_idletasks()
                render_views(
                    model=model,
                    ray=int(ray_trace_var.get()),
                    output_dir_override=out_dir,
                    progress_cb=update_export_progress,
                    image_width=width_val,
                    image_height=height_val,
                    image_units=units_val,
                    dpi_override=dpi_value,
                    min_span_w_override=span_w,
                    min_span_h_override=span_h,
                    faces_to_render=selected_faces,
                )
            except Exception as exc:
                showinfo(title='Error!', message=f'Failed to render views: {exc}')
                return
            export_progress['value'] = 0
            progress_text.set('Ready.')
            export_window.update_idletasks()
            showinfo(title='Done', message='Render complete.')

        def update_viewport_from_inputs():
            ratio = get_output_ratio()
            if ratio is None:
                showinfo(title='Error!', message='Image width and height must be valid numbers.')
                return
            try:
                span_w = float(width_span_var.get())
                span_h = float(height_span_var.get())
            except ValueError:
                showinfo(title='Error!', message='Span values must be numbers.')
                return
            vp_w, vp_h = cmd.get_viewport()
            if vp_h <= 0:
                return
            new_width = max(1, int(round(vp_h * ratio)))
            cmd.viewport(new_width, vp_h)
            min_span = max(span_w, span_h)
            if span_w and (vp_h < new_width):
                min_span = min_span * vp_h / new_width
            if span_h and (new_width < vp_h):
                min_span = min_span * vp_h / new_width
            view_old = list(cmd.get_view())
            cmd.pseudoatom('min_span_cent', pos=view_old[12:15])
            cmd.zoom('min_span_cent', min_span / 2)
            cmd.delete('min_span_cent')
            view_new = list(cmd.get_view())
            view_new[15:17] = view_old[15:17]  # Restore clipping value
            cmd.set_view(view_new)
            try:
                apply_bg_color()
            except Exception:
                pass

        model_label = Label(export_window, text="Select lattice:", width=16, height=2, anchor='e')
        model_menu = OptionMenu(export_window, model_var, *(models or ['']))

        output_label = Label(export_window, text="Output folder:", width=16, height=2, anchor='e')
        output_entry = Entry(export_window, textvariable=output_dir_var, width=40)
        output_browse = Button(export_window, text="Browse Folder...", command=browse_output_dir)

        image_width_label = Label(export_window, text="Image width:", width=16, height=2, anchor='e')
        image_width_entry = Entry(export_window, textvariable=image_width_var, width=12)
        image_height_label = Label(export_window, text="Image height:", width=16, height=2, anchor='e')
        image_height_entry = Entry(export_window, textvariable=image_height_var, width=12)
        units_label = Label(export_window, text="Image Units:", width=16, height=2, anchor='e')
        image_units_menu = OptionMenu(export_window, image_units_var, 'in', 'cm', 'px')

        dpi_label = Label(export_window, text="DPI:", width=16, height=2, anchor='e')
        dpi_menu = OptionMenu(export_window, dpi_var, '96', '150', '300', '600', '1200', 'custom...')
        custom_dpi_entry = Entry(export_window, textvariable=custom_dpi_var, width=12)

        bg_label = Label(export_window, text="Background:", width=12, height=2, anchor='e')
        bg_menu = OptionMenu(export_window, bg_color_var, 'white', 'black', 'custom...')
        custom_bg_entry = Entry(export_window, textvariable=custom_bg_var, width=15)
        
        faces_label = Label(export_window, text="Render faces:", width=16, height=2, anchor='e')
        faces_frame = Frame(export_window)
        face_ab_cb = Checkbutton(faces_frame, text="ab", variable=face_ab_var)
        face_ac_cb = Checkbutton(faces_frame, text="ac", variable=face_ac_var)
        face_bc_cb = Checkbutton(faces_frame, text="bc", variable=face_bc_var)
        
        face_ab_cb.pack(side="left", padx=5)
        face_ac_cb.pack(side="left", padx=5)
        face_bc_cb.pack(side="left", padx=5)

        width_span_label = Label(export_window, text="Distance Width (Angstrom):", width=24, height=2, anchor='e')
        width_span_entry = Entry(export_window, textvariable=width_span_var, width=12)

        height_span_label = Label(export_window, text="Distance Height (Angstrom):", width=24, height=2, anchor='e')
        height_span_entry = Entry(export_window, textvariable=height_span_var, width=12)

        ray_cb = Checkbutton(export_window, text="Ray trace", variable=ray_trace_var)
        progress_label = Label(export_window, textvariable=progress_text, width=24, height=2, anchor='e')
        export_progress = Progressbar(export_window, orient='horizontal', mode='determinate')
        render_button = Button(export_window, text="Render!", command=render_button_press)
        update_viewport_button = Button(export_window, text="Update viewport", command=update_viewport_from_inputs)

        model_label.grid(column=1, row=1)
        model_menu.grid(column=2, row=1, sticky=W)

        output_label.grid(column=1, row=2)
        output_entry.grid(column=2, row=2, columnspan=3, sticky=W+E)
        output_browse.grid(column=5, row=2, padx=10)

        units_label.grid(column=1, row=3)
        image_units_menu.grid(column=2, row=3, sticky=W)

        image_width_label.grid(column=1, row=4)
        image_width_entry.grid(column=2, row=4, sticky=W)
        width_span_label.grid(column=3, row=4)
        width_span_entry.grid(column=4, row=4, sticky=W)

        image_height_label.grid(column=1, row=5)
        image_height_entry.grid(column=2, row=5, sticky=W)
        height_span_label.grid(column=3, row=5)
        height_span_entry.grid(column=4, row=5, sticky=W)

        dpi_label.grid(column=1, row=6)
        dpi_menu.grid(column=2, row=6, sticky=W)
        bg_label.grid(column=3, row=6)
        bg_menu.grid(column=4, row=6, sticky=W)
        
        custom_dpi_entry.grid(column=2, row=7, sticky=W)
        custom_bg_entry.grid(column=4, row=7, sticky=W)
        on_bg_change() # Set initial state
        on_dpi_change() # Set initial state
        
        faces_label.grid(column=1, row=8)
        faces_frame.grid(column=2, row=8, columnspan=2, sticky=W)
        ray_cb.grid(column=4, row=8, sticky=W) # Moved to column 4 to align with custom bg entry

        progress_label.grid(column=1, row=9)
        export_progress.grid(column=2, row=9, columnspan=3, sticky=W+E, padx=10)

        update_viewport_button.grid(column=2, row=10, pady=10, sticky=E, padx=6)
        render_button.grid(column=3, row=10, pady=10, sticky=W, padx=6)

        image_width_var.trace_add('write', update_span_from_image_size)
        image_height_var.trace_add('write', update_span_from_image_size)
        width_span_var.trace_add('write', update_height_span_from_width)
        height_span_var.trace_add('write', update_width_span_from_height)

        def close_export_window():
            export_window_ref['window'] = None
            export_window.destroy()

        export_window.protocol("WM_DELETE_WINDOW", close_export_window)

    def update_input_mode():
        use_file = input_mode_var.get() == 'file'
        entry_state = NORMAL if use_file else DISABLED
        pdb_state = NORMAL if not use_file else DISABLED
        input_entry.configure(state=entry_state)
        browse_button.configure(state=entry_state)
        pdb_entry.configure(state=pdb_state)
        input_label.configure(fg='black' if use_file else 'gray50')
        pdb_label.configure(fg='black' if not use_file else 'gray50')

    def update_progress(value):
        progress_bar['value'] = value
        window.update_idletasks()

    input_label = Label(window, text="Crystal structure from file:", width=24, height=2, anchor='e')
    input_entry = Entry(window, textvariable=input_path_var, width=40)
    browse_button = Button(window, text="Browse File...", command=browse_input_file)

    pdb_label = Label(window, text="PDB code:", width=24, height=2, anchor='e')
    pdb_entry = Entry(window, textvariable=pdb_code_var, width=20)
    input_mode_file = Radiobutton(window, variable=input_mode_var, value='file', command=update_input_mode)
    input_mode_pdb = Radiobutton(window, variable=input_mode_var, value='pdb', command=update_input_mode)

    unit_cell_label = Label(window, text="Unit cell expansion:", width=24, height=2, anchor='e')
    validation = window.register(only_numbers)
    # unit_cell_entry created in layout

    # wrap_cb, color_cb, grid_cb created in layout

    run_button = Button(window, text="Generate lattice!", command=run_button_press)

    view_label = Label(window, text="Set lattice view:", width=24, height=2, anchor='e')
    model_label = Label(window, text="Select lattice:", width=24, height=2, anchor='e')
    model_menu = OptionMenu(window, selected_model_var, '')
    view_ab = Button(window, text="ab", width=6, command=lambda: set_lattice_view('ab'))
    view_ac = Button(window, text="ac", width=6, command=lambda: set_lattice_view('ac'))
    view_bc = Button(window, text="bc", width=6, command=lambda: set_lattice_view('bc'))
    refresh_model_menu()
    selected_model_var.trace_add('write', on_model_change)
    update_view_buttons()

    export_button = Button(window, text="Export views...", command=open_export_window)
    delete_button = Button(window, text="Delete Lattice", command=delete_button_press)
    update_view_buttons()

    input_mode_file.grid(column=0, row=1, sticky=E)
    input_label.grid(column=1, row=1)
    input_entry.grid(column=2, row=1, sticky=W+E)
    browse_button.grid(column=3, row=1, padx=10)

    input_mode_pdb.grid(column=0, row=2, sticky=E)
    pdb_label.grid(column=1, row=2)
    pdb_entry.grid(column=2, row=2, sticky=W+E)

    unit_cell_label.grid(column=1, row=3)

    # Frame for Unit Cell Entry and Wrap Checkbutton
    unit_cell_frame = Frame(window)
    unit_cell_frame.grid(column=2, row=3, sticky=W+E)
    
    # Re-create unit_cell_entry in the new frame
    unit_cell_entry = Entry(unit_cell_frame, validate="key", validatecommand=(validation, '%S'), textvariable=unit_cell_var, width=10)
    unit_cell_entry.pack(side="left")
    
    # Re-create wrap_cb in the new frame
    wrap_cb = Checkbutton(unit_cell_frame, text="Wrap atoms into unit cell", variable=wrap_atoms_var)
    wrap_cb.pack(side="left", padx=10)

    status_label = Label(window, text="Ready!", width=24, height=2, anchor='e')
    progress_bar = Progressbar(window, orient='horizontal', mode='determinate')
    status_label.grid(column=1, row=5)
    progress_bar.grid(column=2, row=5, sticky=W+E, padx=10)

    run_button.grid(column=3, row=5, padx=10, sticky=W)

    separator = Separator(window, orient='horizontal')
    separator.grid(column=1, columnspan=3, row=6, sticky=W+E, padx=10, pady=6)

    model_label.grid(column=1, row=7)
    model_menu.grid(column=2, row=7, sticky=W+E)
    delete_button.grid(column=3, row=7, padx=10, sticky=W)

    # Unit Cell Info Label
    unit_cell_info_label = Label(window, textvariable=unit_cell_info_var, justify="center", fg="blue")
    unit_cell_info_label.grid(column=2, row=8, sticky=W+E, pady=2)

    view_label.grid(column=1, row=9)
    view_ab.grid(column=2, row=9, sticky=W)
    view_ac.grid(column=2, row=9)
    view_bc.grid(column=2, row=9, sticky=E)
    export_button.grid(column=3, row=9, padx=10, sticky=W)

    # Frame for Color and Grid Checkbuttons + Parameters
    options_frame = Frame(window)
    options_frame.grid(column=2, row=10, sticky=W+E)
    
    check_frame = Frame(options_frame)
    check_frame.pack(anchor='center')
    
    # Re-create checkbuttons in the new frame
    show_model_cb = Checkbutton(check_frame, text="Show model", variable=show_model_var, command=apply_display_options)
    color_cb = Checkbutton(check_frame, text="Color by symmetry operator", variable=color_asu_var, command=apply_display_options)
    grid_cb = Checkbutton(check_frame, text="Show grid points", variable=show_grid_var, command=apply_display_options)
    
    show_model_cb.pack(side="left", padx=10)
    color_cb.pack(side="left", padx=10)
    grid_cb.pack(side="left", padx=10)

    scheme_frame = Frame(options_frame)
    scheme_frame.pack(anchor='center', pady=5)
    Label(scheme_frame, text="Color scheme:").pack(side="left", padx=5)
    scheme_menu = OptionMenu(scheme_frame, color_scheme_var, *color_scheme_map.keys())
    scheme_menu.configure(state=NORMAL if color_asu_var.get() else DISABLED)
    scheme_menu.pack(side="left", padx=5)

    # Inner frame for parameters
    params_frame = Frame(options_frame)
    params_frame.pack(anchor='center', pady=5)
    
    Label(params_frame, text="Grid point size:").pack(side="left", padx=5)
    grid_size_entry = Entry(params_frame, textvariable=grid_point_size_var, width=5)
    grid_size_entry.pack(side="left", padx=5)
    grid_size_entry.bind('<Return>', apply_display_options)
    grid_size_entry.bind('<FocusOut>', apply_display_options)

    Label(params_frame, text="Cell thickness:").pack(side="left", padx=5)
    thickness_entry = Entry(params_frame, textvariable=cell_thickness_var, width=5)
    thickness_entry.pack(side="left", padx=5)
    thickness_entry.bind('<Return>', apply_display_options)
    thickness_entry.bind('<FocusOut>', apply_display_options)

    color_scheme_var.trace_add('write', lambda *_args: apply_display_options())

    input_entry.bind('<FocusIn>', lambda _evt: input_mode_var.set('file') or update_input_mode())
    pdb_entry.bind('<FocusIn>', lambda _evt: input_mode_var.set('pdb') or update_input_mode())
    update_input_mode()

    print("[pymol_lattice_viewer] entering Tk mainloop...")
    window.mainloop()
    print("[pymol_lattice_viewer] Tk mainloop exited.")


def launch_lattice_viewer_safe():
    try:
        print("[pymol_lattice_viewer] launching GUI...")
        launch_lattice_viewer()
    except Exception:
        import traceback
        print("[pymol_lattice_viewer] GUI launch failed:")
        print(traceback.format_exc())


cmd.extend("pymol_lattice_viewer", launch_lattice_viewer_safe)
print("[pymol_lattice_viewer] loaded.")
launch_lattice_viewer_safe()
