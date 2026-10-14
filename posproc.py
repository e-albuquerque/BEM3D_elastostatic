import numpy as np
import matplotlib.pyplot as plt
import re

def compute_analytical_static(z, L=6.0, E=2.5e6, nu=0.25, p=100.0):
    """
    Computes the analytical 1D vertical displacement for a prism/cube of height L
    with constrained lateral roller supports (uniaxial strain condition)
    under uniform surface traction p.

    The constrained modulus M (oedometric modulus) is:
        M = E * (1 - nu) / ((1 + nu) * (1 - 2*nu))
        u_z(z) = (p / M) * z
    """
    M = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu))
    uz = (p / M) * z
    return uz, M

def parse_cube8_out(filename='cube8_out.dat'):
    """
    Parses node coordinates and computed displacements from cube8_out.dat.
    """
    nodes_coords = {}
    nodes_disp = {}
    internal_pts_disp = {}

    section = None
    with open(filename, 'r') as f:
        for line in f:
            line_str = line.strip()
            if not line_str:
                continue

            if 'BOUNDARY NODES COORDINATES' in line_str:
                section = 'COORDS'
                continue
            elif 'ELEMENT CONNECTIVITIES' in line_str:
                section = 'CONNECT'
                continue
            elif 'BOUNDARY NODES DISPLACEMENTS' in line_str:
                section = 'DISP'
                continue
            elif 'BOUNDARY ELEMENTS TRACTIONS' in line_str:
                section = 'TRACTIONS'
                continue
            elif 'INTERNAL POINTS DISPLACEMENTS' in line_str:
                section = 'INTERNAL'
                continue
            elif '*******************************************************************************' in line_str:
                if section == 'INTERNAL':
                    section = None
                continue

            parts = line_str.split()

            if section == 'COORDS':
                if len(parts) == 4 and parts[0].isdigit():
                    nid = int(parts[0])
                    coords = [float(parts[1]), float(parts[2]), float(parts[3])]
                    nodes_coords[nid] = coords

            elif section == 'DISP':
                if len(parts) == 4 and parts[0].isdigit():
                    nid = int(parts[0])
                    disp = [float(parts[1]), float(parts[2]), float(parts[3])]
                    nodes_disp[nid] = disp

            elif section == 'INTERNAL':
                if len(parts) == 4 and parts[0].isdigit():
                    pid = int(parts[0])
                    disp = [float(parts[1]), float(parts[2]), float(parts[3])]
                    internal_pts_disp[pid] = disp

    return nodes_coords, nodes_disp, internal_pts_disp

def main():
    L = 6.0          # Comprimento do domínio cúbico (m)
    nu = 0.25        # Coeficiente de Poisson
    mu = 1.0e6       # Módulo de cisalhamento (N/m^2)
    p = 100.0        # Tração normal uniforme aplicada (N/m^2)
    E = 2.0 * mu * (1.0 + nu) # E = 2.5e6 N/m^2

    uz_top_analyt, M = compute_analytical_static(L, L=L, E=E, nu=nu, p=p)
    uz_mid_analyt, _ = compute_analytical_static(3.0, L=L, E=E, nu=nu, p=p)

    print("=" * 65)
    print(" QUADESTH: Pós-processamento e Validação Elastoestática 3D")
    print("=" * 65)
    print(f"Propriedades do material:")
    print(f"  Módulo de Young (E)      = {E:.2e} N/m^2")
    print(f"  Módulo de cisalhamento (mu) = {mu:.2e} N/m^2")
    print(f"  Coeficiente de Poisson (nu) = {nu:.2f}")
    print(f"  Módulo confinado (M)     = {M:.2e} N/m^2")
    print(f"  Tração aplicada no topo  = {p:.1f} N/m^2")
    print(f"  Solução analítica topo (z=6m): u_z = {uz_top_analyt:.8e} m")
    print(f"  Solução analítica meio (z=3m): u_z = {uz_mid_analyt:.8e} m")
    print("-" * 65)

    # 1. Carrega dados do nó monitorado (disp_trac.dat)
    try:
        data_track = np.loadtxt('disp_trac.dat')
        ux_num, uy_num, uz_num = data_track[0], data_track[1], data_track[2]
        tx_num, ty_num, tz_num = data_track[3], data_track[4], data_track[5]
    except Exception as e:
        print(f"Erro ao ler disp_trac.dat: {e}")
        return

    err_uz_top = abs(uz_num - uz_top_analyt) / uz_top_analyt * 100.0

    print("Resultados no Nó Monitorado (Nó 770, centro da face superior z=6m):")
    print(f"  Deslocamento u_x = {ux_num:14.6e} m")
    print(f"  Deslocamento u_y = {uy_num:14.6e} m")
    print(f"  Deslocamento u_z = {uz_num:14.6e} m  (Analítico: {uz_top_analyt:14.6e} m)")
    print(f"  Tração t_x       = {tx_num:14.6e} N/m^2")
    print(f"  Tração t_y       = {ty_num:14.6e} N/m^2")
    print(f"  Tração t_z       = {tz_num:14.6e} N/m^2")
    print(f"  -> Erro relativo em u_z (topo): {err_uz_top:.4f} %")
    print("-" * 65)

    # 2. Carrega perfil completo de nós a partir de cube8_out.dat
    nodes_coords, nodes_disp, int_pts = parse_cube8_out('cube8_out.dat')

    if int_pts:
        uz_mid_num = int_pts[1][2]
        err_uz_mid = abs(uz_mid_num - uz_mid_analyt) / uz_mid_analyt * 100.0
        print(f"Resultado no Ponto Interno 1 (z = 3.0 m, centro do cubo):")
        print(f"  Deslocamento u_z = {uz_mid_num:14.6e} m  (Analítico: {uz_mid_analyt:14.6e} m)")
        print(f"  -> Erro relativo em u_z (meio): {err_uz_mid:.4f} %")
        print("-" * 65)

    # 3. Gráficos de Validação
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Painel 1: Perfil de deslocamento u_z ao longo da altura z
    z_dense = np.linspace(0, L, 200)
    uz_dense, _ = compute_analytical_static(z_dense, L=L, E=E, nu=nu, p=p)

    ax1.plot(uz_dense * 1e4, z_dense, 'r-', linewidth=2.5, label='Solução Analítica: $u_z(z) = \\frac{p \\cdot z}{M}$')

    # Plota nós de contorno
    if nodes_coords and nodes_disp:
        z_vals = [nodes_coords[nid][2] for nid in nodes_coords if nid in nodes_disp]
        uz_vals = [nodes_disp[nid][2] for nid in nodes_coords if nid in nodes_disp]
        ax1.scatter(np.array(uz_vals) * 1e4, z_vals, color='royalblue', s=18, alpha=0.5, label='Nós BEM 3D (QUADESTH)')

    # Destaca nó monitorado e ponto interno
    ax1.plot(uz_num * 1e4, L, 's', color='darkblue', markersize=9, label=f'Nó Topo 770 ($z=6$m, erro={err_uz_top:.3f}%)')
    if int_pts:
        ax1.plot(uz_mid_num * 1e4, 3.0, '^', color='darkgreen', markersize=9, label=f'Pto. Interno ($z=3$m, erro={err_uz_mid:.3f}%)')

    ax1.set_xlabel('Deslocamento Vertical $u_z \\times 10^4$ (m)', fontsize=11)
    ax1.set_ylabel('Cota $Z$ (m)', fontsize=11)
    ax1.set_title('Perfil de Deslocamento Vertical $u_z(z)$', fontsize=12, pad=10)
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.legend(fontsize=9, loc='upper left')

    # Painel 2: Deslocamento no topo - Comparação BEM vs Analítico
    labels = ['Analítico', 'BEM 3D (Topo)', 'BEM 3D (Pto. Int.)']
    vals = [uz_top_analyt * 1e4, uz_num * 1e4, (uz_mid_num if int_pts else 0.0) * 1e4]
    colors = ['#d9534f', '#337ab7', '#5cb85c']

    bars = ax2.bar(labels, vals, color=colors, width=0.45, edgecolor='black', linewidth=1)
    for bar in bars:
        h = bar.get_height()
        ax2.annotate(f'{h:.5f}',
                     xy=(bar.get_x() + bar.get_width() / 2, h),
                     xytext=(0, 4), textcoords="offset points",
                     ha='center', va='bottom', fontsize=10, fontweight='bold')

    ax2.set_ylabel('Deslocamento Vertical $u_z \\times 10^4$ (m)', fontsize=11)
    ax2.set_title('Concordância Numérica vs. Analítica', fontsize=12, pad=10)
    ax2.set_ylim(0, max(vals) * 1.25)
    ax2.grid(True, axis='y', linestyle='--', alpha=0.6)

    plt.suptitle('Validação da Formulação Elastoestática 3D (QUADESTH - Solução de Kelvin)', fontsize=14, y=0.98)
    plt.tight_layout()

    output_plot = 'result_plot.png'
    plt.savefig(output_plot, dpi=300)
    print(f"Gráfico de validação salvo com sucesso em '{output_plot}'.")
    print("=" * 65)

if __name__ == '__main__':
    main()
