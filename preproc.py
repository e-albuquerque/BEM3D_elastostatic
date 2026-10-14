import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.patches as mpatches
import meshio
import shutil

def input_data():
    """
    Defines input parameters for the 3D Boundary Element Method (BEM) simulation
    (QUADESTH - 3D elastostatics for a 6m x 6m x 6m cube).

    Returns:
        dict: A dictionary containing boundary conditions, material properties, and mesh file name.
    """
    # Boundary conditions for each surface:
    # sup1 (z = 6): Upper face - Uniform normal traction p = 100 N/m^2; zero shear tractions
    # sup3 (z = 0): Lower face - Normal displacement constrained (uz = 0); zero shear tractions
    # sup2 (x = 0) and sup4 (x = 6): Lateral faces - Normal displacement constrained (ux = 0); zero shear tractions
    # sup5 (y = 0) and sup6 (y = 6): Lateral faces - Normal displacement constrained (uy = 0); zero shear tractions
    bound_cond = {
        'sup1': {
            'type_x': 1, 'value_x': 0.0,
            'type_y': 1, 'value_y': 0.0,
            'type_z': 1, 'value_z': 100.0
        },
        'sup3': {
            'type_x': 1, 'value_x': 0.0,
            'type_y': 1, 'value_y': 0.0,
            'type_z': 0, 'value_z': 0.0
        },
        'sup2': {
            'type_x': 0, 'value_x': 0.0,
            'type_y': 1, 'value_y': 0.0,
            'type_z': 1, 'value_z': 0.0
        },
        'sup4': {
            'type_x': 0, 'value_x': 0.0,
            'type_y': 1, 'value_y': 0.0,
            'type_z': 1, 'value_z': 0.0
        },
        'sup5': {
            'type_x': 1, 'value_x': 0.0,
            'type_y': 0, 'value_y': 0.0,
            'type_z': 1, 'value_z': 0.0
        },
        'sup6': {
            'type_x': 1, 'value_x': 0.0,
            'type_y': 0, 'value_y': 0.0,
            'type_z': 1, 'value_z': 0.0
        }
    }

    # Material properties:
    # Shear modulus mu = 1.0e6 N/m^2, Poisson's ratio nu = 0.25 => E = 2 * mu * (1 + nu) = 2.5e6 N/m^2
    E = 2.5e6
    nu = 0.25
    file_name = 'placa'

    return {
        'bound_cond': bound_cond,
        'E': E,
        'nu': nu,
        'file_name': file_name
    }

def compute_inodes(file_name, bound_cond):
    """
    Reads mesh file and extracts nodes, elements, and 3D boundary conditions.
    """
    mesh = meshio.read(file_name + '.msh')
    coordinates = mesh.points
    elem = mesh.cells_dict['quad9']
    segments = mesh.cell_data_dict['gmsh:physical']['quad9']
    edges = mesh.field_data

    bc_info = {}
    for key, value in bound_cond.items():
        segment = mesh.field_data[key][0]
        bc_info[key] = {
            'segment': segment,
            'type_x': value['type_x'],
            'value_x': value['value_x'],
            'type_y': value['type_y'],
            'value_y': value['value_y'],
            'type_z': value['type_z'],
            'value_z': value['value_z'],
            'edge': edges[key]
        }

    boundary_nodes = np.unique(elem)

    return {
        'nodes': boundary_nodes,
        'coordinates': coordinates,
        'elements': elem,
        'segments': segments,
        'bc_info': bc_info
    }

def find_node_by_coord(computed_data, target_coord=(3.0, 3.0, 6.0)):
    """
    Finds the boundary node closest to the specified target coordinates.
    Returns:
        tuple: (1-based node ID, actual coordinates array)
    """
    coords = computed_data['coordinates']
    b_nodes = computed_data['nodes']
    pts = coords[b_nodes]
    dists = np.linalg.norm(pts - np.array(target_coord), axis=1)
    min_idx = np.argmin(dists)
    closest_node_0based = b_nodes[min_idx]
    node_id_1based = int(closest_node_0based + 1)
    return node_id_1based, coords[closest_node_0based]

def mount_bcs(segments, bc_info):
    """
    Constructs boundary condition matrices for 3D elements (directions X, Y, Z).
    """
    num_elements = len(segments)

    type_bcs = np.zeros((num_elements, 3), dtype=np.int64)
    val_bcs = np.zeros((num_elements, 3), dtype=np.float64)

    for i, segment_index in enumerate(segments):
        for bc_name, bc_data in bc_info.items():
            if bc_data['segment'] == segment_index:
                type_bcs[i, 0] = bc_data['type_x']
                val_bcs[i, 0] = bc_data['value_x'] if np.isscalar(bc_data['value_x']) else bc_data['value_x'][0]

                type_bcs[i, 1] = bc_data['type_y']
                val_bcs[i, 1] = bc_data['value_y'] if np.isscalar(bc_data['value_y']) else bc_data['value_y'][0]

                type_bcs[i, 2] = bc_data['type_z']
                val_bcs[i, 2] = bc_data['value_z'] if np.isscalar(bc_data['value_z']) else bc_data['value_z'][0]
                break

    return type_bcs, val_bcs

def format_coord(val):
    val_rounded = round(val, 4)
    if val_rounded.is_integer():
        return f"{int(val_rounded)}."
    return f"{val_rounded}"

def create_input_file(filename, computed_data, int_points, type_bcs, val_bcs, input_params, tracked_node=770):
    """
    Creates an input file for the Fortran BEM program based on computed mesh data,
    boundary conditions, and tracked node on the upper face.
    """
    with open(filename, 'w') as f:
        # 1. Title Line
        f.write("CUBIC DOMAIN UNDER UNIFORM TRACTION - 3D ELASTOSTATICS (QUADESTH)\n")

        # Extract data
        nodes = computed_data['coordinates']
        boundary_nodes = computed_data['nodes']
        num_nodes = len(boundary_nodes)
        elements = computed_data['elements']
        num_elements = elements.shape[0]
        num_internal_points = len(int_points) if int_points is not None else 0

        shear_modulus_real = input_params['E'] / (2 * (1 + input_params['nu']))
        poisson_ratio = input_params['nu']

        # 2. Basic Parameter Line (Elastostatics: num_nodes, num_elements, num_internal_points, shear_modulus, poisson_ratio)
        f.write(f" {num_nodes},{num_elements},{num_internal_points},{shear_modulus_real:.1f},{poisson_ratio:.2f}\n")

        # 3. Boundary Nodes Coordinates Lines
        for node_idx in boundary_nodes:
            x, y, z = nodes[node_idx]
            f.write(f" {node_idx + 1} {format_coord(x)} {format_coord(y)} {format_coord(z)}\n")

        # 4. Element Connectivities Lines (9-node quadratic elements)
        # Convert from Gmsh [c0, c1, c2, c3, m01, m12, m23, m30, center]
        # to Dominguez [c0, m01, c1, m12, c2, m23, c3, m30, center]
        gmsh_to_dominguez = [0, 4, 1, 5, 2, 6, 3, 7, 8]
        for i, elem in enumerate(elements):
            elem_ordered = elem[gmsh_to_dominguez]
            nodes_str = " ".join(str(idx + 1) for idx in elem_ordered)
            f.write(f" {i + 1} 9 {nodes_str}\n")

        # 5. Boundary Conditions Lines (3 lines per element: X, Y, Z)
        for i in range(num_elements):
            # Direction X
            bc_x_str = " ".join(f"{val_bcs[i, 0]:.4f}" for _ in range(9))
            f.write(f" {i + 1} {type_bcs[i, 0]} {bc_x_str}\n")

            # Direction Y
            bc_y_str = " ".join(f"{val_bcs[i, 1]:.4f}" for _ in range(9))
            f.write(f" {i + 1} {type_bcs[i, 1]} {bc_y_str}\n")

            # Direction Z
            bc_z_str = " ".join(f"{val_bcs[i, 2]:.4f}" for _ in range(9))
            f.write(f" {i + 1} {type_bcs[i, 2]} {bc_z_str}\n")

        # 6. Internal Points Coordinates Lines
        if num_internal_points > 0:
            for p in int_points:
                f.write(f" {format_coord(p[0])} {format_coord(p[1])} {format_coord(p[2])}\n")

        # 7. Tracked Node printed on the last line
        f.write(f" {tracked_node}\n")


def generate_preprocessing(tracked_node=None, target_coord=(3.0, 3.0, 6.0), output_filename='cube8.dat'):
    """
    Main preprocessing execution function.
    """
    inp_data = input_data()
    computed_data = compute_inodes(inp_data['file_name'], inp_data['bound_cond'])
    type_bcs, val_bcs = mount_bcs(computed_data['segments'], computed_data['bc_info'])

    if tracked_node is None:
        tracked_node, found_coord = find_node_by_coord(computed_data, target_coord)
        print(f"Detected middle node on upper face: Node {tracked_node} at {found_coord}")
    else:
        coords = computed_data['coordinates'][tracked_node - 1]
        print(f"Using specified tracked node: Node {tracked_node} at {coords}")

    # Internal point at domain center (3m, 3m, 3m)
    int_points = np.array([[3.0, 3.0, 3.0]])

    # Generate primary input file
    create_input_file(output_filename, computed_data, int_points,
                      type_bcs, val_bcs, inp_data, tracked_node=tracked_node)
    print(f"Generated {output_filename} successfully.")

    # Also generate quad8.dat as requested
    if output_filename != 'quad8.dat':
        shutil.copyfile(output_filename, 'quad8.dat')
        print("Created copy as quad8.dat successfully.")

    return computed_data, tracked_node


def plot_mesh_physical_entities(msh_filename='placa.msh', save_filename='mesh_faces.png', show=True):
    """
    Plota em 3D todas as faces da malha do Gmsh, destacando cada face com uma cor distinta
    de acordo com os nomes das Physical Entities definidas no Gmsh.
    """
    if not msh_filename.endswith('.msh'):
        msh_filename += '.msh'

    mesh = meshio.read(msh_filename)
    coords = mesh.points
    field_data = mesh.field_data

    tag_to_name = {val[0]: name for name, val in field_data.items()}

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    cmap = plt.get_cmap('tab10')
    handles = []

    cell_types_2d = ['quad9', 'quad8', 'quad', 'triangle6', 'triangle']

    all_physical_tags = set()
    for ctype in cell_types_2d:
        if ctype in mesh.cell_data_dict.get('gmsh:physical', {}):
            all_physical_tags.update(mesh.cell_data_dict['gmsh:physical'][ctype])

    sorted_tags = sorted(list(all_physical_tags))
    tag_counter = 0

    for tag in sorted_tags:
        name = tag_to_name.get(tag, f'Physical Tag {tag}')
        color = cmap(tag_counter % 10)
        tag_counter += 1

        total_elements_in_tag = 0
        tag_verts = []

        for ctype in cell_types_2d:
            if ctype not in mesh.cells_dict:
                continue
            if ctype not in mesh.cell_data_dict.get('gmsh:physical', {}):
                continue

            cells = mesh.cells_dict[ctype]
            phys_tags = mesh.cell_data_dict['gmsh:physical'][ctype]
            matching_idx = np.where(phys_tags == tag)[0]

            if len(matching_idx) == 0:
                continue

            matching_cells = cells[matching_idx]
            total_elements_in_tag += len(matching_cells)

            if ctype in ('quad9', 'quad8', 'quad'):
                verts = [coords[cell[:4]] for cell in matching_cells]
            elif ctype in ('triangle6', 'triangle'):
                verts = [coords[cell[:3]] for cell in matching_cells]
            else:
                verts = [coords[cell] for cell in matching_cells]

            tag_verts.extend(verts)

        if tag_verts:
            poly = Poly3DCollection(
                tag_verts,
                facecolors=color,
                edgecolors='black',
                linewidths=0.6,
                alpha=0.8
            )
            ax.add_collection3d(poly)

            patch = mpatches.Patch(
                facecolor=color,
                edgecolor='black',
                label=f"{name} ({total_elements_in_tag} elementos)"
            )
            handles.append(patch)

    x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
    y_min, y_max = coords[:, 1].min(), coords[:, 1].max()
    z_min, z_max = coords[:, 2].min(), coords[:, 2].max()

    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])
    ax.set_zlim([z_min, z_max])
    ax.set_box_aspect([1, 1, 1])

    ax.set_xlabel('X (m)', fontsize=11, labelpad=8)
    ax.set_ylabel('Y (m)', fontsize=11, labelpad=8)
    ax.set_zlabel('Z (m)', fontsize=11, labelpad=8)
    ax.set_title(f'Malha Gmsh ({msh_filename}) - Faces por Physical Entity', fontsize=13, pad=15)

    if handles:
        ax.legend(
            handles=handles,
            loc='upper left',
            bbox_to_anchor=(1.05, 1.0),
            title='Physical Entities (Gmsh)',
            fontsize=10,
            title_fontsize=11
        )

    plt.tight_layout()

    if save_filename:
        plt.savefig(save_filename, dpi=300, bbox_inches='tight')
        print(f"Visualização das faces salva em '{save_filename}'.")

    if show:
        plt.show()

    return fig, ax


if __name__ == '__main__':
    # 1. Executa pré-processamento para gerar o arquivo cube8.dat (nó 770 rastreado)
    computed_data, tracked_node = generate_preprocessing(tracked_node=770)

    # 2. Visualiza as faces da malha em 3D coloridas por Physical Entity
    plot_mesh_physical_entities('placa.msh', save_filename='mesh_faces.png', show=False)
