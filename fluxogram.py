# Fluxograma da arquitetura e pipeline completo do código QUADESTH (BEM 3D Elastostático)
import os
import graphviz

# Garante que o executável 'dot' do MSYS2 seja encontrado no Windows, se necessário
msys_bin = r"C:\msys64\ucrt64\bin"
if os.path.isdir(msys_bin) and msys_bin not in os.environ.get("PATH", ""):
    os.environ["PATH"] = msys_bin + os.pathsep + os.environ.get("PATH", "")

# Cria um grafo direcionado (Digraph) com layout de Cima para Baixo (TB)
dot = graphviz.Digraph(
    name='QUADESTH_Flowchart',
    comment='Flowchart BEM 3D Elastostatics (QUADESTH)',
    format='png'
)

dot.attr(
    rankdir='TB',
    compound='true',
    splines='spline',
    nodesep='0.38',
    ranksep='0.46',
    pad='0.25',
    dpi='300',
    fontname='Helvetica',
    bgcolor='white'
)
dot.attr('edge', color='#495057', arrowsize='0.75', penwidth='1.2', fontname='Helvetica', fontsize='9')

# Estilos visuais padronizados (evoluídos a partir do diagrama original)
oval_attr = {
    'shape': 'ellipse',
    'style': 'filled',
    'fillcolor': '#e8e8ff',
    'color': '#7986cb',
    'penwidth': '1.4',
    'fontname': 'Helvetica',
    'fontsize': '10.5',
    'margin': '0.18,0.08'
}

circle_attr = {
    'shape': 'circle',
    'width': '0.14',
    'height': '0.14',
    'fixedsize': 'true',
    'style': 'filled',
    'fillcolor': '#e8e8ff',
    'color': '#7986cb',
    'penwidth': '1.2',
    'label': ''
}

py_attr = {
    'shape': 'box',
    'style': 'filled,rounded',
    'fillcolor': '#d1ecf1',
    'color': '#3a9fbf',
    'penwidth': '1.4',
    'fontname': 'Helvetica',
    'fontsize': '10',
    'margin': '0.18,0.10'
}

node_attr = {
    'shape': 'box',
    'style': 'filled,rounded',
    'fillcolor': '#d4edda',
    'color': '#5cb874',
    'penwidth': '1.4',
    'fontname': 'Helvetica',
    'fontsize': '10',
    'margin': '0.18,0.10'
}

omp_attr = {
    'shape': 'box',
    'style': 'filled,rounded',
    'fillcolor': '#cfe2ff',
    'color': '#4a86e8',
    'penwidth': '1.5',
    'fontname': 'Helvetica',
    'fontsize': '10',
    'margin': '0.18,0.10'
}

aux_attr = {
    'shape': 'box',
    'style': 'filled,rounded,dashed',
    'fillcolor': '#f1f3f5',
    'color': '#868e96',
    'penwidth': '1.1',
    'fontname': 'Helvetica',
    'fontsize': '9',
    'margin': '0.14,0.08'
}

file_attr = {
    'shape': 'parallelogram',
    'style': 'filled',
    'fillcolor': '#fff3cd',
    'color': '#d6a228',
    'penwidth': '1.3',
    'skew': '0.18',
    'fontname': 'Helvetica',
    'fontsize': '9.5',
    'margin': '0.12,0.07'
}

# Nós de Início e Fim
dot.node('START', '<<B>BEGIN</B>>', group='spine', **oval_attr)

# ==============================================================================
# 1. SUBGRAFO: PRÉ-PROCESSAMENTO (Python + Gmsh)
# ==============================================================================
with dot.subgraph(name='cluster_preproc') as pre:
    pre.attr(
        label='<<B>1. PRE-PROCESSING (Python / Gmsh)</B>>',
        labeljust='l',
        fontname='Helvetica',
        fontsize='10.5',
        style='dashed,rounded',
        color='#3a9fbf',
        bgcolor='#f7fcfd',
        margin='14'
    )
    pre.node(
        'MSH',
        '<<B>placa.geo / placa.msh</B><BR/><FONT POINT-SIZE="8.5" COLOR="#444444">3D Volume &amp; quad9 Mesh</FONT>>',
        **file_attr
    )
    pre.node(
        'PREPROC',
        '<<B>preproc.py</B><BR/><FONT POINT-SIZE="8.5" COLOR="#333333">Extracts quad9 boundary, internal<BR/>nodes, node permutation &amp; BCs</FONT>>',
        group='spine',
        **py_attr
    )
    pre.node(
        'MESH_IMG',
        '<<B>mesh_faces.png</B><BR/><FONT POINT-SIZE="8.5" COLOR="#444444">3D Mesh Plot</FONT>>',
        **file_attr
    )

    with pre.subgraph() as s_pre:
        s_pre.attr(rank='same')
        s_pre.node('MSH')
        s_pre.node('PREPROC')
        s_pre.node('MESH_IMG')

    pre.edge('MSH', 'PREPROC')
    pre.edge('PREPROC', 'MESH_IMG')

# Arquivo de entrada gerado para o solver Fortran
dot.node(
    'INP_FILE',
    '<<B>cube8.dat</B><BR/><FONT POINT-SIZE="8.5" COLOR="#444444">Nodes, Connectivities, BCs<BR/>&amp; 1965 Internal Points</FONT>>',
    group='spine',
    **file_attr
)

# ==============================================================================
# 2. SUBGRAFO: PROGRAMA PRINCIPAL FORTRAN 90 (QUADESTH)
# ==============================================================================
with dot.subgraph(name='cluster_main') as c:
    c.attr(
        label='<<B>2. MAIN PROGRAM: quadesth_modern (main.f90)</B>>',
        labeljust='r',
        fontname='Helvetica',
        fontsize='10.5',
        style='dashed,rounded',
        color='#6c757d',
        bgcolor='#fafafa',
        margin='16'
    )

    # Espinha dorsal de execução (elastostática pura, sem loop de frequência)
    c.node(
        'INPUT',
        '<<B>input_data</B> <FONT POINT-SIZE="8.5" COLOR="#555555">[io_mod]</FONT><BR/><FONT POINT-SIZE="8.5" COLOR="#333333">Reads model data, elastic constants<BR/>&amp; allocates real64 dynamic arrays</FONT>>',
        group='spine',
        **node_attr
    )
    c.node(
        'PRECOMP',
        '<<B>precompute_element_geometry</B> <FONT POINT-SIZE="8.5" COLOR="#555555">[assembly_mod]</FONT><BR/><FONT POINT-SIZE="8.5" COLOR="#333333">Precomputes shape functions, Jacobians<BR/>&amp; unit normals at 36 Gauss points</FONT>>',
        group='spine',
        **node_attr
    )
    c.node(
        'GHMAT',
        '<<B>ghmat8</B> <FONT POINT-SIZE="8.5" COLOR="#1d4ed8">[assembly_mod — OpenMP]</FONT><BR/><FONT POINT-SIZE="8.5" COLOR="#333333">Parallel collocation assembly of [G] &amp; [H]<BR/>and rigid-body translation diagonal</FONT>>',
        group='spine',
        **omp_attr
    )
    c.node(
        'SOLVER',
        '<<B>solve_linear_system</B> <FONT POINT-SIZE="8.5" COLOR="#1d4ed8">[solver_mod — LAPACK]</FONT><BR/><FONT POINT-SIZE="8.5" COLOR="#333333">Solves real dense system [A]{x} = {FI}<BR/>via double-precision DGESV</FONT>>',
        group='spine',
        **omp_attr
    )
    c.node(
        'REORG',
        '<<B>reorganize_solution</B> <FONT POINT-SIZE="8.5" COLOR="#555555">[io_mod]</FONT><BR/><FONT POINT-SIZE="8.5" COLOR="#333333">Reconstructs boundary displacements (FI)<BR/>and element surface tractions (DFI)</FONT>>',
        group='spine',
        **node_attr
    )
    c.node(
        'BSTRESS',
        '<<B>compute_boundary_stresses</B> <FONT POINT-SIZE="8.5" COLOR="#1d4ed8">[stress_mod — OpenMP]</FONT><BR/><FONT POINT-SIZE="8.5" COLOR="#333333">Traction-Recovery Method on 9-node<BR/>Lagrange elements &amp; nodal averaging</FONT>>',
        group='spine',
        **omp_attr
    )
    c.node(
        'INTPTS',
        '<<B>compute_internal_points</B> <FONT POINT-SIZE="8.5" COLOR="#1d4ed8">[stress_mod — OpenMP]</FONT><BR/><FONT POINT-SIZE="8.5" COLOR="#333333">Somigliana displacement &amp; stress identities<BR/>for all internal volume points (if l &gt; 0)</FONT>>',
        group='spine',
        **omp_attr
    )
    c.node(
        'OUTPUT',
        '<<B>output_data / output_boundary_stresses /<BR/>output_internal_points / save_disp_trac</B><BR/><FONT POINT-SIZE="8.5" COLOR="#333333">Exports simulation results &amp; deallocates memory</FONT>>',
        group='spine',
        **node_attr
    )

    # Sub-rotinas auxiliares / Kernels conectadas à direita
    c.node(
        'AUX_DATA',
        '<<B>data_mod</B><BR/><FONT POINT-SIZE="8" COLOR="#444444">allocate_bem_data &amp;<BR/>update_material_constants</FONT>>',
        **aux_attr
    )
    c.node(
        'AUX_ETAXJA',
        '<<B>etaxja</B><BR/><FONT POINT-SIZE="8" COLOR="#444444">Intrinsic Lagrange<BR/>geometry &amp; normals</FONT>>',
        **aux_attr
    )
    c.node(
        'AUX_KELVIN',
        '<<B>extin8 / locin8 / funkelvin</B><BR/><FONT POINT-SIZE="8" COLOR="#444444">Regular &amp; polar singular<BR/>3D Kelvin quadrature (U*, T*)</FONT>>',
        **aux_attr
    )
    c.node(
        'AUX_LAPACK',
        '<<B>LAPACK DGESV</B><BR/><FONT POINT-SIZE="8" COLOR="#444444">LU factorization with<BR/>partial pivoting (BLAS)</FONT>>',
        **aux_attr
    )
    c.node(
        'AUX_BSTR',
        '<<B>boundary_stresses_elem</B><BR/><FONT POINT-SIZE="8" COLOR="#444444">shape_deriv_lagrange &amp;<BR/>calc_principal_and_vm</FONT>>',
        **aux_attr
    )
    c.node(
        'AUX_SOMIG',
        '<<B>Kelvin Kernels (U*, T*, D*, S*)</B><BR/><FONT POINT-SIZE="8" COLOR="#444444">Pre-weighted Gauss integration<BR/>&amp; calc_principal_and_vm</FONT>>',
        **aux_attr
    )

    # Alinhamento lateral à direita da espinha principal
    for main_n, aux_n in [
        ('INPUT', 'AUX_DATA'),
        ('PRECOMP', 'AUX_ETAXJA'),
        ('GHMAT', 'AUX_KELVIN'),
        ('SOLVER', 'AUX_LAPACK'),
        ('BSTRESS', 'AUX_BSTR'),
        ('INTPTS', 'AUX_SOMIG'),
    ]:
        with c.subgraph() as s_pair:
            s_pair.attr(rank='same')
            s_pair.node(main_n)
            s_pair.node(aux_n)
        c.edge(main_n, aux_n, dir='back', style='dashed', color='#868e96')

    # Conexões principais dentro do MAIN PROGRAM
    c.edge('INPUT', 'PRECOMP')
    c.edge('PRECOMP', 'GHMAT')
    c.edge('GHMAT', 'SOLVER')
    c.edge('SOLVER', 'REORG')
    c.edge('REORG', 'BSTRESS')
    c.edge('BSTRESS', 'INTPTS')
    c.edge('INTPTS', 'OUTPUT')

# ==============================================================================
# ARQUIVOS DE SAÍDA DO SOLVER FORTRAN
# ==============================================================================
dot.node(
    'OUT_MAIN',
    '<<B>cube8_out.dat</B><BR/><FONT POINT-SIZE="8" COLOR="#444444">Full BEM Report</FONT>>',
    **file_attr
)
dot.node(
    'OUT_STR',
    '<<B>stresses.dat</B><BR/><FONT POINT-SIZE="8" COLOR="#444444">Boundary Stresses</FONT>>',
    **file_attr
)
dot.node(
    'OUT_INT',
    '<<B>internal_points.dat</B><BR/><FONT POINT-SIZE="8" COLOR="#444444">Internal u &amp; σ</FONT>>',
    **file_attr
)
dot.node(
    'OUT_TRK',
    '<<B>disp_trac.dat</B><BR/><FONT POINT-SIZE="8" COLOR="#444444">Tracked Node</FONT>>',
    **file_attr
)

with dot.subgraph() as s_out:
    s_out.attr(rank='same')
    s_out.node('OUT_MAIN')
    s_out.node('OUT_STR')
    s_out.node('OUT_INT')
    s_out.node('OUT_TRK')
    s_out.edge('OUT_MAIN', 'OUT_STR', style='invis')
    s_out.edge('OUT_STR', 'OUT_INT', style='invis')
    s_out.edge('OUT_INT', 'OUT_TRK', style='invis')

dot.node('JOIN_OUT', group='spine', **circle_attr)

# ==============================================================================
# 3. SUBGRAFO: PÓS-PROCESSAMENTO E VALIDAÇÃO (Python)
# ==============================================================================
with dot.subgraph(name='cluster_posproc') as pos:
    pos.attr(
        label='<<B>3. POST-PROCESSING &amp; VALIDATION (Python)</B>>',
        labeljust='r',
        fontname='Helvetica',
        fontsize='10.5',
        style='dashed,rounded',
        color='#3a9fbf',
        bgcolor='#f7fcfd',
        margin='14'
    )
    pos.node(
        'POSPROC',
        '<<B>posproc.py</B><BR/><FONT POINT-SIZE="8.5" COLOR="#333333">Compares numerical u &amp; σ against<BR/>exact 3D analytical solutions</FONT>>',
        group='spine',
        **py_attr
    )
    pos.node(
        'PLOT_IMG',
        '<<B>result_plot.png</B><BR/><FONT POINT-SIZE="8.5" COLOR="#444444">Validation Plots</FONT>>',
        **file_attr
    )

    with pos.subgraph() as s_pos:
        s_pos.attr(rank='same')
        s_pos.node('POSPROC')
        s_pos.node('PLOT_IMG')

    pos.edge('POSPROC', 'PLOT_IMG')

dot.node('END', '<<B>END</B>>', group='spine', **oval_attr)

# ==============================================================================
# CONEXÕES ENTRE OS ESTÁGIOS DO PIPELINE
# ==============================================================================
dot.edge('START', 'PREPROC')
dot.edge('PREPROC', 'INP_FILE')
dot.edge('INP_FILE', 'INPUT')

dot.edge('OUTPUT', 'OUT_MAIN')
dot.edge('OUTPUT', 'OUT_STR')
dot.edge('OUTPUT', 'OUT_INT')
dot.edge('OUTPUT', 'OUT_TRK')

dot.edge('OUT_MAIN', 'JOIN_OUT', arrowhead='none')
dot.edge('OUT_STR', 'JOIN_OUT', arrowhead='none')
dot.edge('OUT_INT', 'JOIN_OUT', arrowhead='none')
dot.edge('OUT_TRK', 'JOIN_OUT', arrowhead='none')
dot.edge('JOIN_OUT', 'POSPROC')

dot.edge('POSPROC', 'END')

# ==============================================================================
# LEGENDA COMPACTA HORIZONTAL NA BASE (2 LINHAS x 3 COLUNAS)
# ==============================================================================
with dot.subgraph(name='cluster_legend') as leg:
    leg.attr(
        label='<<B>Legend</B>>',
        labeljust='l',
        fontname='Helvetica',
        fontsize='10',
        color='#adb5bd',
        style='rounded',
        bgcolor='#ffffff',
        margin='10'
    )
    leg.node('leg_se', 'Start / End', **oval_attr)
    leg.node('leg_py', 'Python Script', **py_attr)
    leg.node('leg_file', 'Input / Output File', **file_attr)
    leg.node('leg_sub', 'Fortran Subroutine', **node_attr)
    leg.node('leg_omp', 'OpenMP / LAPACK Routine', **omp_attr)
    leg.node('leg_aux', 'Auxiliary Kernel / Module', **aux_attr)

    with leg.subgraph() as leg_r1:
        leg_r1.attr(rank='same')
        leg_r1.node('leg_se')
        leg_r1.node('leg_py')
        leg_r1.node('leg_file')
        leg_r1.edge('leg_se', 'leg_py', style='invis')
        leg_r1.edge('leg_py', 'leg_file', style='invis')

    with leg.subgraph() as leg_r2:
        leg_r2.attr(rank='same')
        leg_r2.node('leg_sub')
        leg_r2.node('leg_omp')
        leg_r2.node('leg_aux')
        leg_r2.edge('leg_sub', 'leg_omp', style='invis')
        leg_r2.edge('leg_omp', 'leg_aux', style='invis')

    leg.edge('leg_se', 'leg_sub', style='invis')
    leg.edge('leg_py', 'leg_omp', style='invis')
    leg.edge('leg_file', 'leg_aux', style='invis')

dot.edge('END', 'leg_py', style='invis')

output_name = 'flowchart'
dot.render(output_name, format='png', cleanup=True)
print(f"Fluxograma exportado com sucesso como '{output_name}.png'!")
