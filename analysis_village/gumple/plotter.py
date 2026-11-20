#!/usr/bin/python3
# Imports
from scipy.interpolate import interp1d
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from cycler import cycler
import scipy.stats
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
import argparse
import json
import os

extended_colors = [
    # Original matplotlib defaults (tab10)
    "#1f77b4",  # blue
    "#ff7f0e",  # orange
    "#2ca02c",  # green
    "#d62728",  # red
    "#9467bd",  # purple
    "#8c564b",  # brown
    "#e377c2",  # pink
    "#7f7f7f",  # gray
    "#bcbd22",  # olive
    "#17becf",  # cyan
    # Extensions
    "#aec7e8",  # light blue
    "#ffbb78",  # light orange
    "#98df8a",  # light green
    "#ff9896",  # light red
    "#c5b0d5",  # light purple
    "#c49c94",  # light brown
    "#f7b6d2",  # light pink
    "#c7c7c7",  # light gray
    "#dbdb8d",  # light olive
    "#9edae5",  # light cyan
    "#393b79",  # navy
    "#843c39",  # maroon
    "#637939",  # moss
    "#e7ba52",  # gold
]

#def envelope(x_per_file, y_per_file, y_min=-1, y_max=2, n_points=500):
#    y_vals = np.logspace(y_min, y_max, n_points)
#    nfiles = len(x_per_file)
#    for y_val in y_vals:
#        for i in range(nfiles): 

def line_from_csv(file_path, x_column, y_column):
    try:
        # Load the data
        df = pd.read_csv(file_path)
        # Check if columns exist
        if x_column not in df.columns or y_column not in df.columns:
            print(f"Error: Columns '{x_column}' or '{y_column}' not found in CSV.")
            return

    except FileNotFoundError:
        print("Error: The file was not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

    x_unsort = df[x_column].to_numpy()
    y_unsort = df[y_column].to_numpy()

    sort_ind = np.argsort(y_unsort)

    return x_unsort[sort_ind], y_unsort[sort_ind]

def plot_line_from_csv(file_path, label, x_column, y_column, color):
    try:
        # Load the data
        df = pd.read_csv(file_path)
        # Check if columns exist
        if x_column not in df.columns or y_column not in df.columns:
            print(f"Error: Columns '{x_column}' or '{y_column}' not found in CSV.")
            return

        plt.plot(df[x_column], df[y_column], label=label, lw=3, linestyle='--', color=color)
        
    except FileNotFoundError:
        print("Error: The file was not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
    

def plotProposal(file_path, c='black'):
    """
    Reads x and y coordinates from a JSON file, sorts them by y-value, 
    and creates a log-log line plot.

    Args:
        file_path (str): The path to the JSON file.
    """
    try:
        # Open the JSON file and load its contents
        with open(file_path, 'r') as f:
            data = json.load(f)

        # Extract 'x' and 'y' lists from the data
        x_coords = data.get('x')
        y_coords = data.get('y')

        # Check if the data is valid
        if x_coords is None or y_coords is None:
            print("Error: JSON file must contain 'x' and 'y' keys.")
            return
        if len(x_coords) != len(y_coords):
            print("Error: 'x' and 'y' lists must have the same length.")
            return

        # Combine x and y coordinates into a list of pairs
        points = zip(x_coords, y_coords)

        # Sort the points based on the y-value (the second element of each pair)
        sorted_points = sorted(points, key=lambda point: point[1])

        # Unzip the sorted points back into separate lists
        x_sorted, y_sorted = zip(*sorted_points)

        # Create a new figure and axis for the plot
        # plt.figure(figsize=(10, 6))
        if "nova" in file_path or "minos" in file_path:
            x_sorted = np.sin(2*np.arcsin(np.sqrt(x_sorted)))**2

        # Create the line plot with the sorted data
        plt.plot(x_sorted, y_sorted, linewidth=3, color=c, ls='--')

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from the file '{file_path}'.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


def plotFile(FileStr, linestyleStr, labelStr, colorStr, level):
  with open(FileStr) as f:
      skip = 0
      adding_params = False
      params = []
      for line in f:
          if line[:4] == 'xval':
              skip += 1
              if not 'p0' in line:
                  params = []
              break
          elif 'Parameters:' in line:
              skip += 1
              adding_params = True 
          elif adding_params == True and line != '\n':
              skip += 1
              params.append(line)
          elif adding_params == True and line == '\n':
              skip += 1
              adding_params = False
          else:
              skip += 1

  vals = np.genfromtxt(FileStr, skip_header=skip)
  params = ['chi'] + params
  params = ['yval'] + params
  params = ['xval'] + params
  df = pd.DataFrame(vals, columns=params)
  x = np.power(10, vals[:,0])
  y = np.power(10, vals[:,1])
  z = vals[:,2]

  # Create a regular grid
  xi = np.unique(x)
  yi = np.unique(y)
  xg, yg = np.meshgrid(xi, yi)
  z = df.pivot(index='xval', columns='yval', values='chi').to_numpy()
  interp_func = RegularGridInterpolator((xi, yi), z, method='pchip')
  zi = interp_func((xg, yg))
  nlevels = [level]

  plt.ylim(max(min(yi), 0.1), min(max(yi), 100))
  plt.xlim(max(min(xi), 0.01), min(max(xi), 1))

  #plt.contourf(xi, yi, zi, levels=[nlevels[0], np.inf], colors=[colorStr], alpha=0.3)
  plt.contour(xi, yi, zi, nlevels, linestyles=linestyleStr, colors = [colorStr], linewidths=3)

DEFAULT_BASE_DIR = "/home/nathanielerowe/SBN/PROfit"

def main(files, labels, linestyles, colors, output, proposal_path, sigma_bool, icarus_dir, excluded, title, level, points):
  if labels:
    labels = [label.encode('utf-8').decode('unicode_escape') for label in labels]
  fig=plt.figure(figsize=(8,7), dpi=150)#plt.figure(figsize=(8,7), dpi=150, layout='constrained')
  ax=plt.subplot(111)
  for axis in ['top','bottom','left','right']:
      ax.spines[axis].set_linewidth(2)
  if sigma_bool and title == "$90\\%$ CL PROfit for SBN":
      plt.title("$5\\sigma$ PROfit for SBN", fontsize=30, pad=40)
  else:
      plt.title(title, fontsize=20)
  plt.tick_params(size=8, which='minor')
  plt.tick_params(size=15, labelsize=15, which='major')
  ax.set_xlabel(r'$\sin^2 2 \theta_{\mu \mu}$',fontsize=20)
  ax.set_ylabel(r'$\Delta m_{41}^2$  [eV$^2$]',fontsize=20)
  plt.loglog()

  if not colors:
      def_colors = extended_colors#list(mcolors.TABLEAU_COLORS.keys())*2
  else:
      def_colors = colors

  if not linestyles:
      def_linestyles = ['solid']*len(labels)
  else:
      def_linestyles = linestyles 

  patches = []
  i = 0
  for f in files:
      plotFile(f, def_linestyles[i], "Stat", def_colors[i], level)
      patch = mpatches.Patch(color=def_colors[i], label=labels[i], fill=False, linewidth=3, ls=def_linestyles[i])
      patches.append(patch)
      i += 1
  if proposal_path is not None:
      # Resolve default path based on sigma mode if no explicit path was given
      if proposal_path == "":
          proposal_path = os.path.join(DEFAULT_BASE_DIR, "tn.json" if sigma_bool else "data.json")
      plotProposal(proposal_path)
      patch = mpatches.Patch(color='black', label='Proposal', fill=False, linewidth=3)
      patches.append(patch)

  if icarus_dir is not None:
      tools = icarus_dir if icarus_dir != "" else os.path.join(DEFAULT_BASE_DIR, "tools")

      plot_line_from_csv(os.path.join(tools, "spine_icarus_run2.csv"), "ICARUS\nSPINE", 'x', 'y', 'black')
      patch = mpatches.Patch(color='black', label='ICARUS\nSPINE', fill=False, linewidth=3, ls='--')
      patches.append(patch)

      plotProposal(os.path.join(tools, "minos.json"), c='purple')
      patch = mpatches.Patch(color='purple', label='MINOS', fill=False, linewidth=3, ls='--')
      patches.append(patch)

      plot_line_from_csv(os.path.join(tools, "nova_2024_mumu.csv"), "NOVA\n2024", 'x', 'y', 'cyan')
      patch = mpatches.Patch(color='cyan', label='NOVA\n2024', fill=False, linewidth=3, ls='--')
      patches.append(patch)

      plot_line_from_csv(os.path.join(tools, "icecube_2024.csv"), "IceCube\n2024", 'x', 'y', 'brown')
      patch = mpatches.Patch(color='brown', label='IceCube\n2024', fill=False, linewidth=3, ls='--')
      patches.append(patch)

      plot_line_from_csv(os.path.join(tools, "numu_MiniBoone_Sciboone_2012.csv"), "MiniBooNE", 'x', 'y', 'gold')
      patch = mpatches.Patch(color='gold', label='MiniBooNE', fill=False, linewidth=3, ls='--')
      patches.append(patch)

      #plot_line_from_csv(os.path.join(tools, "numu_CCFR_2012.csv"), "CCFR", 'x', 'y', 'red')
      #patch = mpatches.Patch(color='red', label='CCFR', fill=False, linewidth=3, ls='--')
      #patches.append(patch)

  if excluded:
      patch = mpatches.Patch(color='red', label='Excluded', fill=False, linewidth=3, linestyle='-', lw=2)
      patches.append(patch)
      tools = icarus_dir if (icarus_dir is not None and icarus_dir != "") else os.path.join(DEFAULT_BASE_DIR, "tools")

      x2cat = []
      y2cat = []

      for f in ["spine_icarus_run2.csv", "nova_2024_mumu.csv", "numu_MiniBoone_Sciboone_2012.csv", "numu_CCFR_2012.csv"]:
          x, y = line_from_csv(os.path.join(tools, f), 'x', 'y')

          # add these for closure
          x2cat.append([1.0])
          y2cat.append([y[0]])

          x2cat.append(x)
          y2cat.append(y)

          # add these for closure
          x2cat.append([1.0])
          y2cat.append([y2cat[-1][-1]])

      file_path = os.path.join(tools, "minos.json")
      with open(file_path, 'r') as f:
          data = json.load(f)

      # Extract 'x' and 'y' lists from the data
      x = data.get('x')
      x = np.sin(2*np.arcsin(np.sqrt(x)))**2
      y = data.get('y')

      # add these for closure
      x2cat.append([1.0])
      y2cat.append([y[0]])

      x2cat.append(x)
      y2cat.append(y)

      # add these for closure
      x2cat.append([1.0])
      y2cat.append([y2cat[-1][-1]])

      all_x = np.concatenate(x2cat)
      all_y = np.concatenate(y2cat)
      poly = plt.fill_betweenx(all_y, all_x, 1.0, alpha=0.2, color='red', edgecolor='none')

      for i in range(len(x2cat)):
          inds = np.argsort(y2cat[i])
          x2cat[i] = np.array(x2cat[i])[inds]
          y2cat[i] = np.array(y2cat[i])[inds]

      x_verts = []
      y_verts = []
      glob_y_smooth = np.linspace(-1.0, 2.0, 500)
      for x, y in zip(x2cat, y2cat):
          unique_y, indices = np.unique(y, return_index=True)
          unique_x = x[indices] 
          
          # 2. Convert to Log-Space
          # We interpolate Log(X) vs Log(Y) because it's a log-log plot
          log_x = np.log10(unique_x)
          log_y = np.log10(unique_y)
          
          # 3. Create the Interpolation Function
          # 'cubic' gives a smooth curve; 'linear' gives a clean connect-the-dots
          f_interp = interp1d(log_y, log_x, kind='linear', fill_value="extrapolate")
         
          y_smooth_log = glob_y_smooth[(glob_y_smooth < log_y.max()) & (glob_y_smooth > log_y.min())]
          # 4. Generate smooth points for plotting
          x_smooth_log = f_interp(y_smooth_log)
          
          # 5. Convert back to Linear Space
          x_smooth = 10**x_smooth_log
          y_smooth = 10**y_smooth_log

          x_verts.append(x_smooth)
          y_verts.append(y_smooth)

          # 6. Plot the interpolated line
          #plt.plot(x_smooth, y_smooth, color='red', lw=2, label='Interpolated Edge')

      all_x_smooth = np.concatenate(x_verts)
      all_y_smooth = np.concatenate(y_verts)
      #path = poly.get_paths()[0]
      #verts = path.vertices
      # 1. Capture the collection from your plot call
      # poly = ax.fill_betweenx(y_data, x1_data, x2_data, ...)
      
      # 2. Extract vertices
      # This gets the actual [x, y] coordinates used to draw the shape
      path = poly.get_paths()[0]
      vertices = path.vertices  # Array of shape (N, 2)
      
      # 3. Find the minimum X for every unique Y
      # We use a dictionary or pandas to group by Y and find the min X
      unique_ys = np.unique(all_y_smooth)
      left_edge_x = []
      left_edge_y = []
      
      for y in unique_ys:
          # Find all x-coordinates associated with this specific y-level
          mask = (all_y_smooth == y)
          min_x = np.min(all_x_smooth[mask])
          
          left_edge_x.append(min_x)
          left_edge_y.append(y)
      
      # 4. Plot the clean boundary
      plt.plot(left_edge_x, left_edge_y, color='red', lw=2, linestyle='-', label='Excluded') # uncomment for envelope

      #envelope_x, envelope_y = envelope(x2cat, y2cat)
      #plt.plot(envelope_x, envelope_y, color='red', lw=2, alpha=1.0) 

  point_patches = []
  if points:
      for x, y, label, color in points:
          plt.scatter([x], [y], color=color, zorder=5)
          point_patches.append(mlines.Line2D([], [], color=color, marker='o', linestyle='None',
                                             markersize=8, label=label))

  # Interleave contour patches with point patches so each row of the legend
  # groups related entries: [contour_a, contour_b, point] per row (ncol=3).
  # This only applies when the number of contour patches is an exact multiple
  # of the number of points (e.g. 6 contours + 3 points → 3 rows of 3).
  # Otherwise fall back to contours first, then points.
  n_pts = len(point_patches)
  n_con = len(patches)
  if n_pts > 0 and n_con % n_pts == 0:
      stride = n_con // n_pts          # contour entries per point (e.g. 2)
      legend_handles = []
      for i, pp in enumerate(point_patches):
          legend_handles.extend(patches[i * stride:(i + 1) * stride])
          legend_handles.append(pp)
      ncol = stride + 1                # e.g. 3 columns
  else:
      legend_handles = patches + point_patches
      ncol = min(len(legend_handles), 3)

  plt.legend(handles=legend_handles, fontsize=16, loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=ncol, framealpha=0)
  legend = ax.get_legend()
  legend_height = legend.get_window_extent().height / fig.dpi
  fig_height = 7 + legend_height
  fig.set_size_inches(8, fig_height)
  plt.savefig(output, bbox_inches='tight', pad_inches=0.3)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='SensPlotter',
        description='Plotter for txt file outputs from PROfit surfaces.',
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument('-f', '--file', nargs='+', type=str,
                        help='Input PROfit surface file(s).')
    parser.add_argument('-l', '--label', nargs='+', type=str,
                        help='Legend labels, one per --file.')
    parser.add_argument('-ls', '--linestyle', nargs='+', type=str,
                        help='Matplotlib linestyles, one per --file.')
    parser.add_argument('-c', '--color', nargs='+', type=str,
                        help='Colors, one per --file.')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Output image path.')
    parser.add_argument('-t', '--title', type=str, default="$90\\%$ CL PROfit for SBN",
                        help='Plot title. (default: 90%% CL PROfit for SBN)')
    parser.add_argument('-s', '--sigma', action='store_true',
                        help='Use 5-sigma level (23.66) and title. Overridden by --level.')
    parser.add_argument('--level', type=float, default=None,
                        help='Contour delta-chi2 level (1-sided, 1 DOF). Overrides --sigma.\n'
                             'Common values: 1.64 (90%% CL), 0.23 (1σ), 2.86 (2σ),\n'
                             '               7.74 (3σ), 14.70 (4σ), 23.40 (5σ)\n'
                             '(default: 23.40 with --sigma, else 1.64)')
    parser.add_argument('-p', '--proposal', nargs='?', const='', default=None,
                        metavar='PATH',
                        help='Overlay proposal contour. Optionally provide a JSON path;\n'
                             'defaults to <base>/tn.json (--sigma) or <base>/data.json.')
    parser.add_argument('-i', '--icarus', nargs='?', const='', default=None,
                        metavar='DIR',
                        help='Overlay ICARUS/MINOS/NOvA/MiniBooNE/CCFR limits.\n'
                             'Optionally provide a tools directory path;\n'
                             'defaults to <base>/tools/.')
    parser.add_argument('-e', '--excluded', action='store_true',
                        help='Shade the excluded region from external limits.')
    parser.add_argument('--point', nargs=4, action='append',
                        metavar=('X', 'Y', 'LABEL', 'COLOR'),
                        help='Plot a point at (X, Y) with the given label and color.\n'
                             'Repeatable: --point 0.16 3.5 "Best fit" red --point 0.116 2.07 "Lo dm2" blue')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # Resolve contour level: explicit --level > --sigma > default
    if args.level is not None:
        level = args.level
    elif args.sigma:
        level = 23.40
    else:
        level = 1.64

    # Parse x/y strings to floats for each point, preserve label and color as-is
    points = [(float(x), float(y), label, color)
              for x, y, label, color in args.point] if args.point else []

    main(args.file, args.label, args.linestyle, args.color, args.output,
         args.proposal, args.sigma, args.icarus, args.excluded, args.title, level, points)
