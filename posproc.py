import numpy as np
import matplotlib.pyplot as plt
import os

def compute_analytical_static(z, L=6.0, E=2.5e6, nu=0.25, p=100.0):
    """
    Computes the analytical 1D vertical displacement and stresses for a prism/cube of height L
    with constrained lateral roller supports (uniaxial strain condition)
    under uniform surface traction p.

    The constrained modulus M (oedometric modulus) is:
        M = E * (1 - nu) / ((1 + nu) * (1 - 2*nu))
        u_z(z) = (p / M) * z
        sigma_zz = p
        sigma_xx = sigma_yy = (nu / (1 - nu)) * p
        sigma_vm = sigma_zz - sigma_xx
    """
    M = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu))
    uz = (p / M) * z
    sig_zz = p
    sig_xx = (nu / (1.0 - nu)) * p
    sig_yy = sig_xx
    sig_vm = abs(sig_zz - sig_xx)
    return uz, M, sig_zz, sig_xx, sig_yy, sig_vm

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

def parse_stresses(filename='stresses.dat'):
    """
    Parses boundary nodal stresses from stresses.dat.
    """
    if not os.path.exists(filename):
        return None
    data = np.loadtxt(filename)
    return data

def parse_internal_points(filename='internal_points.dat'):
    """
    Parses internal points displacements and stresses from internal_points.dat.
    Col 0: Point ID
    Col 1-3: X, Y, Z coordinates
    Col 4-6: u_x, u_y, u_z displacements
    Col 7-12: Sigma_xx, Sigma_yy, Sigma_zz, Sigma_xy, Sigma_yz, Sigma_zx
    Col 13: Von Mises Equivalent Stress
    Col 14-16: Principal Stresses (Sigma_1 >= Sigma_2 >= Sigma_3)
    """
    if not os.path.exists(filename):
        return None
    data = np.loadtxt(filename)
    return data

def main():
    L = 6.0          # Comprimento do domínio cúbico (m)
    nu = 0.25        # Coeficiente de Poisson
    mu = 1.0e6       # Módulo de cisalhamento (N/m^2)
    p = 100.0        # Tração normal uniforme aplicada (N/m^2)
    E = 2.0 * mu * (1.0 + nu) # E = 2.5e6 N/m^2

    uz_top_analyt, M, sig_zz_exact, sig_xx_exact, sig_yy_exact, sig_vm_exact = compute_analytical_static(L, L=L, E=E, nu=nu, p=p)
    uz_mid_analyt, _, _, _, _, _ = compute_analytical_static(3.0, L=L, E=E, nu=nu, p=p)

    print("=" * 70)
    print(" QUADESTH: Pós-processamento e Validação Elastoestática 3D")
    print(" Deslocamentos e Tensões de Contorno (Traction-Recovery)")
    print("=" * 70)
    print(f"Propriedades do material:")
    print(f"  Módulo de Young (E)         = {E:.2e} N/m^2")
    print(f"  Módulo de cisalhamento (mu) = {mu:.2e} N/m^2")
    print(f"  Coeficiente de Poisson (nu) = {nu:.2f}")
    print(f"  Módulo confinado (M)        = {M:.2e} N/m^2")
    print(f"  Tração aplicada no topo     = {p:.1f} N/m^2")
    print(f"  Solução analítica topo (z=6m): u_z = {uz_top_analyt:.8e} m")
    print(f"  Tensões analíticas: sigma_zz = {sig_zz_exact:.2f} Pa, sigma_xx = sigma_yy = {sig_xx_exact:.2f} Pa, Von Mises = {sig_vm_exact:.2f} Pa")
    print("-" * 70)

    # 1. Carrega dados do nó monitorado (disp_trac.dat)
    try:
        data_track = np.loadtxt('disp_trac.dat')
        ux_num, uy_num, uz_num = data_track[0], data_track[1], data_track[2]
        tx_num, ty_num, tz_num = data_track[3], data_track[4], data_track[5]
    except Exception as e:
        print(f"Erro ao ler disp_trac.dat: {e}")
        return

    err_uz_top = abs(uz_num - uz_top_analyt) / uz_top_analyt * 100.0

    print("Resultados de Deslocamento no Nó 770 (centro da face superior z=6m):")
    print(f"  Deslocamento u_x = {ux_num:14.6e} m")
    print(f"  Deslocamento u_y = {uy_num:14.6e} m")
    print(f"  Deslocamento u_z = {uz_num:14.6e} m  (Analítico: {uz_top_analyt:14.6e} m)")
    print(f"  Tração t_x       = {tx_num:14.6e} N/m^2")
    print(f"  Tração t_y       = {ty_num:14.6e} N/m^2")
    print(f"  Tração t_z       = {tz_num:14.6e} N/m^2")
    print(f"  -> Erro relativo em u_z (topo): {err_uz_top:.4f} %")
    print("-" * 70)

    # 2. Carrega perfil completo de nós a partir de cube8_out.dat
    nodes_coords, nodes_disp, int_pts = parse_cube8_out('cube8_out.dat')

    # 3. Carrega e valida pontos internos (internal_points.dat)
    int_pts_data = parse_internal_points('internal_points.dat')
    has_int_pts = int_pts_data is not None and len(int_pts_data) > 0

    if has_int_pts:
        pids_int = int_pts_data[:, 0].astype(int)
        x_int = int_pts_data[:, 1]
        y_int = int_pts_data[:, 2]
        z_int = int_pts_data[:, 3]
        ux_int = int_pts_data[:, 4]
        uy_int = int_pts_data[:, 5]
        uz_int = int_pts_data[:, 6]
        sxx_int = int_pts_data[:, 7]
        syy_int = int_pts_data[:, 8]
        szz_int = int_pts_data[:, 9]
        sxy_int = int_pts_data[:, 10]
        syz_int = int_pts_data[:, 11]
        szx_int = int_pts_data[:, 12]
        svm_int = int_pts_data[:, 13]

        # Ponto interno mais próximo do centro geométrico (3.0, 3.0, 3.0)
        dist_center = np.sqrt((x_int - 3.0)**2 + (y_int - 3.0)**2 + (z_int - 3.0)**2)
        idx_c = np.argmin(dist_center)
        zc = z_int[idx_c]
        uz_c_exact = (p / M) * zc
        err_uz_c = abs(uz_int[idx_c] - uz_c_exact) / uz_c_exact * 100.0
        err_szz_c = abs(szz_int[idx_c] - sig_zz_exact) / sig_zz_exact * 100.0
        err_sxx_c = abs(sxx_int[idx_c] - sig_xx_exact) / sig_xx_exact * 100.0
        err_syy_c = abs(syy_int[idx_c] - sig_yy_exact) / sig_yy_exact * 100.0
        err_svm_c = abs(svm_int[idx_c] - sig_vm_exact) / sig_vm_exact * 100.0

        print(f"Resultados no Ponto Interno Central #{pids_int[idx_c]} (Coord: x={x_int[idx_c]:.3f}, y={y_int[idx_c]:.3f}, z={zc:.3f} m):")
        print(f"  Deslocamento u_z = {uz_int[idx_c]:14.6e} m  (Analítico: {uz_c_exact:14.6e} m, Erro: {err_uz_c:.4f}%)")
        print(f"  sigma_zz         = {szz_int[idx_c]:10.4f} Pa  (Analítico: {sig_zz_exact:10.4f} Pa, Erro: {err_szz_c:.4f}%)")
        print(f"  sigma_xx         = {sxx_int[idx_c]:10.4f} Pa  (Analítico: {sig_xx_exact:10.4f} Pa, Erro: {err_sxx_c:.4f}%)")
        print(f"  sigma_yy         = {syy_int[idx_c]:10.4f} Pa  (Analítico: {sig_yy_exact:10.4f} Pa, Erro: {err_syy_c:.4f}%)")
        print(f"  sigma_xy         = {sxy_int[idx_c]:10.4e} Pa, sigma_yz = {syz_int[idx_c]:10.4e} Pa, sigma_zx = {szx_int[idx_c]:10.4e} Pa")
        print(f"  Von Mises        = {svm_int[idx_c]:10.4f} Pa  (Analítico: {sig_vm_exact:10.4f} Pa, Erro: {err_svm_c:.4f}%)")
        print("-" * 70)

        # Erros relativos de deslocamento para z > 0.05m (todos os pontos internos)
        mask_z = z_int > 0.05
        uz_exact_all = (p / M) * z_int
        err_uz_all = np.abs(uz_int[mask_z] - uz_exact_all[mask_z]) / uz_exact_all[mask_z] * 100.0

        print(f"Estatísticas Globais de TODOS os Pontos Internos ({len(z_int)} pontos, incluindo próximos ao contorno):")
        print(f"  u_z:       Erro Médio = {np.mean(err_uz_all):.5f}%, Erro Máx = {np.max(err_uz_all):.5f}%")
        print(f"  sigma_zz:  Média = {np.mean(szz_int):.4f} +/- {np.std(szz_int):.4f} Pa (Analítico: {sig_zz_exact:.2f} Pa)")
        print(f"  sigma_xx:  Média = {np.mean(sxx_int):.4f} +/- {np.std(sxx_int):.4f} Pa (Analítico: {sig_xx_exact:.2f} Pa)")
        print(f"  sigma_yy:  Média = {np.mean(syy_int):.4f} +/- {np.std(syy_int):.4f} Pa (Analítico: {sig_yy_exact:.2f} Pa)")
        print(f"  Von Mises: Média = {np.mean(svm_int):.4f} +/- {np.std(svm_int):.4f} Pa (Analítico: {sig_vm_exact:.2f} Pa)")
        print("-" * 70)

    # 4. Carrega e valida tensões de contorno (stresses.dat)
    stresses = parse_stresses('stresses.dat')
    if stresses is not None:
        row770 = stresses[stresses[:, 0] == 770]
        if len(row770) > 0:
            s770 = row770[0]
            sxx_770, syy_770, szz_770 = s770[4], s770[5], s770[6]
            sxy_770, syz_770, szx_770 = s770[7], s770[8], s770[9]
            svm_770 = s770[10]

            err_szz = abs(szz_770 - sig_zz_exact) / sig_zz_exact * 100.0
            err_sxx = abs(sxx_770 - sig_xx_exact) / sig_xx_exact * 100.0
            err_svm = abs(svm_770 - sig_vm_exact) / sig_vm_exact * 100.0

            print("Tensões de Contorno no Nó 770 (Recuperação no Contorno):")
            print(f"  sigma_zz  = {szz_770:10.4f} Pa  (Analítico: {sig_zz_exact:10.4f} Pa, Erro: {err_szz:.4f}%)")
            print(f"  sigma_xx  = {sxx_770:10.4f} Pa  (Analítico: {sig_xx_exact:10.4f} Pa, Erro: {err_sxx:.4f}%)")
            print(f"  sigma_yy  = {syy_770:10.4f} Pa  (Analítico: {sig_yy_exact:10.4f} Pa, Erro: {abs(syy_770-sig_yy_exact)/sig_yy_exact*100.0:.4f}%)")
            print(f"  sigma_xy  = {sxy_770:10.4e} Pa, sigma_yz = {syz_770:10.4e} Pa, sigma_zx = {szx_770:10.4e} Pa")
            print(f"  Von Mises = {svm_770:10.4f} Pa  (Analítico: {sig_vm_exact:10.4f} Pa, Erro: {err_svm:.4f}%)")
            print("-" * 70)

        szz_all = stresses[:, 6]
        sxx_all = stresses[:, 4]
        syy_all = stresses[:, 5]
        svm_all = stresses[:, 10]
        print(f"Estatísticas Globais das Tensões de Contorno ({len(stresses)} nós):")
        print(f"  sigma_zz:  Média = {np.mean(szz_all):.4f} Pa, DesvPad = {np.std(szz_all):.4f} Pa (Analítico: {sig_zz_exact:.2f} Pa)")
        print(f"  sigma_xx:  Média = {np.mean(sxx_all):.4f} Pa, DesvPad = {np.std(sxx_all):.4f} Pa (Analítico: {sig_xx_exact:.2f} Pa)")
        print(f"  sigma_yy:  Média = {np.mean(syy_all):.4f} Pa, DesvPad = {np.std(syy_all):.4f} Pa (Analítico: {sig_yy_exact:.2f} Pa)")
        print(f"  Von Mises: Média = {np.mean(svm_all):.4f} Pa, DesvPad = {np.std(svm_all):.4f} Pa (Analítico: {sig_vm_exact:.2f} Pa)")
        print("-" * 70)

    # 5. Gráficos de Validação (2x2)
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    ((ax1, ax2), (ax3, ax4)) = axes

    # Painel 1: Perfil de deslocamento u_z ao longo da altura z
    z_dense = np.linspace(0, L, 200)
    uz_dense, _, _, _, _, _ = compute_analytical_static(z_dense, L=L, E=E, nu=nu, p=p)

    ax1.plot(uz_dense * 1e4, z_dense, 'r-', linewidth=2.5, label='Solução Analítica: $u_z(z) = \\frac{p \\cdot z}{M}$')

    if nodes_coords and nodes_disp:
        z_vals = [nodes_coords[nid][2] for nid in nodes_coords if nid in nodes_disp]
        uz_vals = [nodes_disp[nid][2] for nid in nodes_coords if nid in nodes_disp]
        ax1.scatter(np.array(uz_vals) * 1e4, z_vals, color='royalblue', s=16, alpha=0.45, label='Nós de Contorno (QUADESTH)')

    if has_int_pts:
        ax1.scatter(uz_int * 1e4, z_int, color='mediumseagreen', s=12, alpha=0.35, label=f'Pontos Internos ({len(z_int)} pts)')
        ax1.plot(uz_int[idx_c] * 1e4, zc, '^', color='forestgreen', markersize=9, label=f'Centro Interno ({zc:.2f}m, erro={err_uz_c:.4f}%)')

    ax1.plot(uz_num * 1e4, L, 's', color='darkblue', markersize=9, label=f'Nó Topo 770 ($z=6$m, erro={err_uz_top:.4f}%)')

    ax1.set_xlabel('Deslocamento Vertical $u_z \\times 10^4$ (m)', fontsize=11)
    ax1.set_ylabel('Cota $Z$ (m)', fontsize=11)
    ax1.set_title('(a) Perfil de Deslocamento Vertical $u_z(z)$', fontsize=12, pad=10)
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.legend(fontsize=9, loc='upper left')

    # Painel 2: Tensão Vertical sigma_zz
    ax2.axhline(sig_zz_exact, color='black', linestyle='--', linewidth=2.0, label=rf'Analítico $\sigma_{{zz}} = {sig_zz_exact:.1f}$ Pa')
    if stresses is not None:
        z_nodes = stresses[:, 3]
        szz_nodes = stresses[:, 6]
        ax2.scatter(z_nodes, szz_nodes, color='crimson', s=16, alpha=0.4, label=r'Nós de Contorno ($\sigma_{zz}$)')
        if len(row770) > 0:
            ax2.plot(6.0, szz_770, 's', color='darkred', markersize=9, label=rf'Nó 770 ({szz_770:.2f} Pa, erro={err_szz:.4f}%)')

    if has_int_pts:
        ax2.scatter(z_int, szz_int, color='mediumpurple', s=12, alpha=0.45, label=rf'Pontos Internos ({len(z_int)} pts, $\mu$={np.mean(szz_int):.2f} Pa)')
        ax2.plot(zc, szz_int[idx_c], '^', color='darkmagenta', markersize=9, label=rf'Centro ({szz_int[idx_c]:.2f} Pa, erro={err_szz_c:.4f}%)')

    ax2.set_xlabel('Cota $Z$ (m)', fontsize=11)
    ax2.set_ylabel(r'Tensão $\sigma_{zz}$ (Pa)', fontsize=11)
    ax2.set_title(r'(b) Tensão Normal Vertical $\sigma_{zz}$', fontsize=12, pad=10)
    ax2.set_ylim(85.0, 115.0)
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.legend(fontsize=9, loc='lower right')

    # Painel 3: Tensões Laterais sigma_xx e sigma_yy
    ax3.axhline(sig_xx_exact, color='black', linestyle='--', linewidth=2.0, label=rf'Analítico $\sigma_{{xx}}=\sigma_{{yy}} = {sig_xx_exact:.2f}$ Pa')
    if stresses is not None:
        sxx_nodes = stresses[:, 4]
        syy_nodes = stresses[:, 5]
        ax3.scatter(z_nodes, sxx_nodes, color='dodgerblue', s=14, alpha=0.35, label=r'$\sigma_{xx}$ (Contorno)')
        ax3.scatter(z_nodes, syy_nodes, color='mediumseagreen', s=14, alpha=0.35, label=r'$\sigma_{yy}$ (Contorno)')
        if len(row770) > 0:
            ax3.plot(6.0, sxx_770, 'D', color='blue', markersize=8, label=rf'Nó 770 $\sigma_{{xx}}$ ({sxx_770:.2f} Pa)')

    if has_int_pts:
        ax3.scatter(z_int, sxx_int, color='steelblue', s=10, alpha=0.3, label=rf'$\sigma_{{xx}}$ Internos ({len(z_int)} pts)')
        ax3.scatter(z_int, syy_int, color='seagreen', s=10, alpha=0.3, label=rf'$\sigma_{{yy}}$ Internos ({len(z_int)} pts)')
        ax3.plot(zc, sxx_int[idx_c], '^', color='navy', markersize=8, label=rf'Centro $\sigma_{{xx}}$ ({sxx_int[idx_c]:.2f} Pa)')

    ax3.set_xlabel('Cota $Z$ (m)', fontsize=11)
    ax3.set_ylabel(r'Tensões Laterais $\sigma_{xx}, \sigma_{yy}$ (Pa)', fontsize=11)
    ax3.set_title(r'(c) Tensões de Confinamento Lateral $\sigma_{xx}$ e $\sigma_{yy}$', fontsize=12, pad=10)
    ax3.set_ylim(20.0, 45.0)
    ax3.grid(True, linestyle='--', alpha=0.6)
    ax3.legend(fontsize=9, loc='lower right')

    # Painel 4: Tensão Equivalente de Von Mises
    ax4.axhline(sig_vm_exact, color='black', linestyle='--', linewidth=2.0, label=rf'Analítico $\sigma_{{vm}} = {sig_vm_exact:.2f}$ Pa')
    if stresses is not None:
        svm_nodes = stresses[:, 10]
        ax4.scatter(z_nodes, svm_nodes, color='purple', s=16, alpha=0.4, label=r'$\sigma_{vm}$ (Contorno)')
        if len(row770) > 0:
            ax4.plot(6.0, svm_770, '*', color='indigo', markersize=11, label=rf'Nó 770 $\sigma_{{vm}}$ ({svm_770:.2f} Pa, erro={err_svm:.4f}%)')

    if has_int_pts:
        ax4.scatter(z_int, svm_int, color='darkorange', s=12, alpha=0.4, label=rf'$\sigma_{{vm}}$ Internos ({len(z_int)} pts, $\mu$={np.mean(svm_int):.2f} Pa)')
        ax4.plot(zc, svm_int[idx_c], '^', color='chocolate', markersize=9, label=rf'Centro ({svm_int[idx_c]:.2f} Pa, erro={err_svm_c:.4f}%)')

    ax4.set_xlabel('Cota $Z$ (m)', fontsize=11)
    ax4.set_ylabel(r'Tensão de Von Mises $\sigma_{vm}$ (Pa)', fontsize=11)
    ax4.set_title(r'(d) Tensão Equivalente de Von Mises $\sigma_{vm}$', fontsize=12, pad=10)
    ax4.set_ylim(50.0, 85.0)
    ax4.grid(True, linestyle='--', alpha=0.6)
    ax4.legend(fontsize=9, loc='lower right')

    plt.suptitle('Validação da Formulação Elastoestática 3D (QUADESTH - Solução de Kelvin)\nDeslocamentos e Tensões no Contorno e nos Pontos Internos (Somigliana)', fontsize=14, y=0.99)
    plt.tight_layout()

    output_plot = 'result_plot.png'
    plt.savefig(output_plot, dpi=300)
    print(f"Gráfico de validação salvo com sucesso em '{output_plot}'.")
    print("=" * 70)

if __name__ == '__main__':
    main()
