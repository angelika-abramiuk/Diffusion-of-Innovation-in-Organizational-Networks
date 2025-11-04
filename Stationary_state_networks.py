#!/usr/bin/env python3
"""Plot stationary states from simulations on real networks (S1..L2).

Place this script in your repo and run it from the project root (so relative
paths to Results/... work). It writes a PDF into ./figures/.

Requires:
    numpy, matplotlib
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import os
import warnings

# --- Configuration: update paths or p_eng here if needed ---
BASE_DIR = Path(__file__).resolve().parent
NETWORKS = [
    {'name': 'S1', 'folder': BASE_DIR / '2024-11-20 14-34 S1', 'N': 320,  'k': 14.81, 'q_eff': 3.809375},
    {'name': 'M1', 'folder': BASE_DIR / '2024-11-20 17-29 M1', 'N': 1429, 'k': 8.8,   'q_eff': 3.966410},
    {'name': 'L1', 'folder': BASE_DIR / '2025-03-17 15-39 L1 T4000', 'N': 5793, 'k': 27.09, 'q_eff': 2.611255},
    {'name': 'S2', 'folder': BASE_DIR / '2024-11-20 16-30 S2', 'N': 165,  'k': 45.22, 'q_eff': 3.339394},
    {'name': 'M2', 'folder': BASE_DIR / '2024-11-21 19-22 M2', 'N': 3862, 'k': 10.62, 'q_eff': 4.0},
    {'name': 'L2', 'folder': BASE_DIR / '2024-11-23 21-39 L2', 'N': 5524, 'k': 34.11, 'q_eff': 3.999638},
]

# default parameters
Q = 4
P_ENG = 0.75          # change if you want to plot other p_eng
P_IND = np.arange(0, 1.02, 0.02)

# plotting style
try:
    plt.rcParams['text.usetex'] = True
except Exception:
    # systems without LaTeX will fall back silently
    plt.rcParams['text.usetex'] = False

plt.rcParams['font.family'] = 'serif'
plt.rcParams.update({
    'font.size': 24,
    'axes.titlesize': 24,
    'axes.labelsize': 24,
    'legend.fontsize': 24,
})

BLUE = np.array([87, 117, 144]) / 255.0
RED  = np.array([158, 42, 43])   / 255.0

# --- Helper functions ---
def load_stationary_values(input_folder: Path, N: int, q: int, p_eng: float, p_ind_array):
    """Load averaged trajectory last values for each p_ind; return numpy array."""
    values = []
    files_found = 0
    for ip in p_ind_array:
        fname = input_folder / f"Averaged_Trajectory_N_{N}_q_{q}_p_eng_{p_eng:.2f}_p_ind_{ip:.2f}.csv"
        if fname.exists():
            try:
                data = np.loadtxt(fname, delimiter=',')
                if np.size(data) == 0:
                    values.append(np.nan)
                elif data.ndim == 1:
                    values.append(data[-1])
                else:
                    # if 2D, take last element of last row
                    values.append(data[-1].ravel()[-1])
                files_found += 1
            except Exception:
                # try genfromtxt as fallback
                try:
                    data = np.genfromtxt(fname, delimiter=',')
                    if np.size(data) == 0:
                        values.append(np.nan)
                    elif data.ndim == 1:
                        values.append(data[-1])
                    else:
                        values.append(data[-1].ravel()[-1])
                except Exception:
                    values.append(np.nan)
        else:
            values.append(np.nan)
    return np.array(values), files_found

def compute_mfa_pa_curves(q_eff: float, k: float, p_eng: float, npoints=1000):
    """Return (c_mfa, p_mfa, valid_mfa, c_pa, p_pa, valid_pa)."""
    c_mfa = np.linspace(0, 1, npoints)
    with np.errstate(divide='ignore', invalid='ignore'):
        num = c_mfa * (1 - c_mfa)**q_eff - (1 - c_mfa) * c_mfa**q_eff
        den = num - c_mfa + p_eng
        p_mfa = num / den
    valid_mfa = np.isfinite(p_mfa) & (p_mfa >= 0) & (p_mfa <= 1)

    c_pa = np.linspace(0, 1, npoints)
    with np.errstate(divide='ignore', invalid='ignore'):
        term1 = ((1 - c_pa)**q_eff) * p_eng
        term2 = (c_pa**q_eff) * (1 - p_eng)
        numerator = 2 * (c_pa * (1 - c_pa) * (term1 - term2) -
                         (q_eff / k) * (p_eng - c_pa) * (c_pa * (1 - c_pa)**q_eff + (1 - c_pa) * c_pa**q_eff))
        denominator = (term1 - term2 -
                       (q_eff / k) * (p_eng - c_pa) * ((1 - c_pa)**q_eff + c_pa**q_eff))
        b = numerator / denominator
        theta_up = np.where(c_pa > 0, b / (2 * c_pa), np.nan)
        theta_down = np.where((1 - c_pa) > 0, b / (2 * (1 - c_pa)), np.nan)
        with np.errstate(invalid='ignore', divide='ignore'):
            num_pa = c_pa * theta_up**q_eff - (1 - c_pa) * theta_down**q_eff
            den_pa = num_pa - c_pa + p_eng
            p_pa = num_pa / den_pa
    valid_pa = np.isfinite(p_pa) & (p_pa >= 0) & (p_pa <= 1)

    return (c_mfa, p_mfa, valid_mfa, c_pa, p_pa, valid_pa)

# --- Plotting routine ---
def plot_all(networks, q, p_eng, p_ind):
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()

    missing_report = []

    for idx, info in enumerate(networks):
        ax = axes[idx]
        name = info['name']
        N = info['N']
        k = info['k']
        qeff = info['q_eff']
        main_folder = Path(info['folder'])

        input_folder = main_folder / f"Results_Organizations_N_{N}_q_{q}_p_eng_{p_eng:.2f}"
        stationary_values, found = load_stationary_values(input_folder, N, q, p_eng, p_ind)
        missing = len(p_ind) - found
        missing_report.append((name, found, missing, input_folder))

        # compute theory curves
        c_mfa, p_mfa, valid_mfa, c_pa, p_pa, valid_pa = compute_mfa_pa_curves(qeff, k, p_eng)

        # Plot PA/MFA according to p_eng cases
        if abs(p_eng - 0.75) < 1e-9:
            # determine switching points safely
            valid_lower_mfa = valid_mfa & (c_mfa < 0.5)
            idxs_mfa = np.where(valid_lower_mfa)[0]
            if idxs_mfa.size > 0:
                p_vals = p_mfa[idxs_mfa]
                i_rel = int(np.nanargmax(p_vals))
                p_star_mfa = p_vals[i_rel]
                c_star_mfa = c_mfa[idxs_mfa][i_rel]
            else:
                p_star_mfa = np.nan; c_star_mfa = np.nan

            valid_lower_pa = valid_pa & (c_pa < 0.5)
            idxs_pa = np.where(valid_lower_pa)[0]
            if idxs_pa.size > 0:
                p_vals_pa = p_pa[idxs_pa]
                i_rel_pa = int(np.nanargmax(p_vals_pa))
                p_star_pa = p_vals_pa[i_rel_pa]
                c_star_pa = c_pa[idxs_pa][i_rel_pa]
            else:
                p_star_pa = np.nan; c_star_pa = np.nan

            # PA segments
            if not np.isnan(p_star_pa):
                seg1 = valid_pa & (p_pa <= p_star_pa) & (c_pa <= c_star_pa)
                seg2 = valid_pa & (p_pa <= p_star_pa) & (c_pa >= c_star_pa) & (c_pa <= 0.5)
                seg_up = valid_pa & (c_pa >= 0.5)
                ax.plot(p_pa[seg2], c_pa[seg2], 'k--', linewidth=1)
                ax.plot(p_pa[seg1], c_pa[seg1], 'k-', linewidth=1)
                ax.plot(p_pa[seg_up], c_pa[seg_up], 'k-', linewidth=1, label='PA')
            else:
                # fallback: plot valid PA
                ax.plot(p_pa[valid_pa], c_pa[valid_pa], 'k-', linewidth=1, label='PA')

            # MFA segments
            if not np.isnan(p_star_mfa):
                seg1m = valid_mfa & (p_mfa <= p_star_mfa) & (c_mfa <= c_star_mfa)
                seg2m = valid_mfa & (p_mfa <= p_star_mfa) & (c_mfa >= c_star_mfa) & (c_mfa <= 0.5)
                seg_upm = valid_mfa & (c_mfa >= 0.5)
                ax.plot(p_mfa[seg2m], c_mfa[seg2m], 'k--', linewidth=2)
                ax.plot(p_mfa[seg1m], c_mfa[seg1m], 'k-', linewidth=2)
                ax.plot(p_mfa[seg_upm], c_mfa[seg_upm], 'k-', linewidth=2, label='MFA')
            else:
                ax.plot(p_mfa[valid_mfa], c_mfa[valid_mfa], 'k-', linewidth=2, label='MFA')

        elif abs(p_eng - 0.5) < 1e-9:
            # symmetric case: horizontal lines at c_star
            if np.any(valid_mfa):
                p_star_mfa = np.nanmax(p_mfa[valid_mfa])
                c_star_mfa = c_mfa[valid_mfa][np.nanargmax(p_mfa[valid_mfa])]
                ax.plot(p_mfa[valid_mfa], c_mfa[valid_mfa], 'k-', linewidth=2, label='MFA')
                ax.hlines(y=c_star_mfa, xmin=0, xmax=p_star_mfa, colors='k', linestyles='dashed', linewidth=2)
                ax.hlines(y=c_star_mfa, xmin=p_star_mfa, xmax=0.6, colors='k', linewidth=2)
            if np.any(valid_pa):
                p_star_pa = np.nanmax(p_pa[valid_pa])
                c_star_pa = c_pa[valid_pa][np.nanargmax(p_pa[valid_pa])]
                ax.plot(p_pa[valid_pa], c_pa[valid_pa], 'k-', linewidth=1, label='PA')
                ax.hlines(y=c_star_pa, xmin=0, xmax=p_star_pa, colors='k', linestyles='dashed', linewidth=1)
                ax.hlines(y=c_star_pa, xmin=p_star_pa, xmax=0.6, colors='k', linewidth=1)
        else:
            # general fallback: just plot valid curves
            ax.plot(p_pa[valid_pa], c_pa[valid_pa], 'k-', linewidth=1, label='PA')
            ax.plot(p_mfa[valid_mfa], c_mfa[valid_mfa], 'k-', linewidth=2, label='MFA')

        # simulation points
        ax.plot(P_IND := p_ind, stationary_values, marker='o', linestyle='none', markersize=7,
                markerfacecolor=RED, markeredgecolor='k', markeredgewidth=0.8, label=name)

        # formatting
        ax.set_xlim(0, 0.6)
        ax.set_ylim(0, 1)
        ax.grid(False)

        # formatting
        ax.set_xlim(0, 0.6)
        ax.set_ylim(0, 1)
        ax.grid(False)

        # only left column gets y-labels, only bottom row gets x-labels
        if idx % 3 == 0:  # left column
            ax.set_ylabel(r'$c$')
        else:
            ax.set_ylabel('')
            ax.set_yticklabels([])

        if idx // 3 == 1:  # bottom row
            ax.set_xlabel(r'$p^{\mathrm{ind}}$')
        else:
            ax.set_xlabel('')
            ax.set_xticklabels([])

        ax.legend(loc='lower right')

    plt.tight_layout()
    outdir = Path('figures')
    outdir.mkdir(parents=True, exist_ok=True)
    out_file = outdir / f"Stationary_state_networks_p_eng_{p_eng:.2f}_q_{q}_min_q_k_effective_q.pdf"
    # plt.savefig(out_file, bbox_inches='tight')
    # print(f"Saved figure to: {out_file}")
    plt.show()

    # summary of missing files
    # print("\nFiles summary (network, found, missing, expected folder):")
    # for name, found, missing, folder in missing_report:
    #     print(f"  {name:3s}: found={found:3d}, missing={missing:3d}, folder={folder}")

# --- Main ---
if __name__ == "__main__":
    # change P_ENG here or add argparse if you want CLI
    plot_all(NETWORKS, Q, P_ENG, P_IND)
