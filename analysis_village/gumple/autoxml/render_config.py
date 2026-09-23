import os
import sys
import glob
import argparse
from jinja2 import Environment, FileSystemLoader

# Resolve and add ../../gump to sys.path
script_dir = os.path.dirname(os.path.abspath(__file__))
gump_dir = os.path.abspath(os.path.join(script_dir, "../../gump"))

if gump_dir not in sys.path:
    sys.path.insert(0, gump_dir)

# Import grab_pot function from pot.py
from pot import grab_pot

def format_sci(value, precision=3):
    """Formats floats in scientific notation without the '+' in the exponent (e.g. 5.272e20)."""
    return f"{value:.{precision}e}".replace("e+", "e")

def get_sample_pot(file_pattern, use_pot=True):
    """
    Finds all files matching pattern and passes them to grab_pot
    to return total aggregated POT.
    """
    files = glob.glob(file_pattern)
    if not files:
        print(f"Warning: No files found matching pattern: {file_pattern}")
        return 0.0
    
    pot_bools = [use_pot] * len(files)
    total_pot = grab_pot(files, onbeam_bools=pot_bools, sep_bool=False)
    return total_pot

def get_scale(sample_pot, target_pot):
    """Calculates scaling factor: scale = target_pot / sample_pot"""
    if sample_pot == 0:
        return 1.0
    return target_pot / sample_pot

def render_template(template_path, output_path, base_dir, hdf_dir, target_sbnd_pot=1e20, target_icarus_r2_pot=2e20, target_icarus_r4_pot=3e20):
    print("--- Calculating POT and Scale Factors ---")

    # 1. SBND MC (Wildcard)
    sbnd_mc_pattern = os.path.join(hdf_dir, "SBNDMCCV_*.df")
    sbnd_mc_pot = get_sample_pot(sbnd_mc_pattern, use_pot=True)

    # 6. SBND OffBeam (Data)
    sbnd_offbeam_pattern = os.path.join(hdf_dir, "SBND_SpringBNBOffData.df")
    sbnd_offbeam_pot = get_sample_pot(sbnd_offbeam_pattern, use_pot=False)

    # 9. SBND Dirt (MC)
    sbnd_dirt_pattern = os.path.join(hdf_dir, "SBND_SpringLowEMC.df")
    sbnd_dirt_pot = get_sample_pot(sbnd_dirt_pattern, use_pot=True)

    # 3. ICARUS Run 2 MC (Wildcard)
    icarus_r2_mc_pattern = os.path.join(hdf_dir, "ICARUSRun2_SpringMCOverlay_rewgt_*.df")
    icarus_r2_mc_pot = get_sample_pot(icarus_r2_mc_pattern, use_pot=True)
    icarus_r2_mc_scale = get_scale(icarus_r2_mc_pot, target_icarus_r2_pot)

    # 4. ICARUS Run 2 OffBeam (Data)
    icarus_r2_offbeam_pattern = os.path.join(hdf_dir, "ICARUS_SpringRun2BNBOff_unblind.df")
    icarus_r2_offbeam_pot = get_sample_pot(icarus_r2_offbeam_pattern, use_pot=False)
    icarus_r2_offbeam_scale = get_scale(icarus_r2_offbeam_pot, target_icarus_r2_pot)

    # 7. ICARUS Run 2 Dirt (MC)
    icarus_r2_dirt_pattern = os.path.join(hdf_dir, "ICARUSRun2_Spring_Overlay_Dirt.df")
    icarus_r2_dirt_pot = get_sample_pot(icarus_r2_dirt_pattern, use_pot=True)
    icarus_r2_dirt_scale = get_scale(icarus_r2_dirt_pot, target_icarus_r2_pot)

    # 2. ICARUS Run 4 MC (Wildcard)
    icarus_r4_mc_pattern = os.path.join(hdf_dir, "ICARUSRun4_SpringMCOverlay_rewgt_*.df")
    icarus_r4_mc_pot = get_sample_pot(icarus_r4_mc_pattern, use_pot=True)
    icarus_r4_mc_scale = get_scale(icarus_r4_mc_pot, target_icarus_r4_pot)

    # 5. ICARUS Run 4 OffBeam (Data)
    icarus_r4_offbeam_pattern = os.path.join(hdf_dir, "ICARUS_SpringRun4BNBOff_unblind.df")
    icarus_r4_offbeam_pot = get_sample_pot(icarus_r4_offbeam_pattern, use_pot=False)
    icarus_r4_offbeam_scale = get_scale(icarus_r4_offbeam_pot, target_icarus_r4_pot)

    # 8. ICARUS Run 4 Dirt (MC)
    icarus_r4_dirt_pattern = os.path.join(hdf_dir, "ICARUSRun4_Spring_Overlay_Dirt.df")
    icarus_r4_dirt_pot = get_sample_pot(icarus_r4_dirt_pattern, use_pot=True)
    icarus_r4_dirt_scale = get_scale(icarus_r4_dirt_pot, target_icarus_r4_pot)

    print("\n--- Rendering Jinja2 Configuration ---")

    env = Environment(
        loader=FileSystemLoader('.'),
        comment_start_string='{%#',
        comment_end_string='#%}'
    )
    template = env.get_template(template_path)

    context = {
        "BASE_DIR": base_dir,
        "sbnd_target_pot": format_sci(target_sbnd_pot, 1),
        "icarus_target_pot": format_sci(target_icarus_r2_pot+target_icarus_r4_pot, 1),
        
        # SBND POT values
        "sbnd_mc_pot": format_sci(sbnd_mc_pot, 3),
        "sbnd_offbeam_pot": format_sci(sbnd_offbeam_pot, 3),
        "sbnd_dirt_pot": format_sci(sbnd_dirt_pot, 3),

        # ICARUS Run 2 scales & raw POTs
        "icarus_r2_mc_pot": format_sci(icarus_r2_mc_pot, 3),
        "icarus_r2_mc_scale": f"{icarus_r2_mc_scale:.3f}",
        "icarus_r2_offbeam_scale": f"{icarus_r2_offbeam_scale:.3f}",
        "icarus_r2_dirt_scale": f"{icarus_r2_dirt_scale:.3f}",

        # ICARUS Run 4 scales & raw POTs
        "icarus_r4_mc_pot": format_sci(icarus_r4_mc_pot, 3),
        "icarus_r4_mc_scale": f"{icarus_r4_mc_scale:.3f}",
        "icarus_r4_offbeam_scale": f"{icarus_r4_offbeam_scale:.3f}",
        "icarus_r4_dirt_scale": f"{icarus_r4_dirt_scale:.3f}",
    }

    rendered_content = template.render(**context)

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(rendered_content)
    
    print(f"Successfully written configuration to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Render XML config with dataset POT calculations.")
    parser.add_argument("-t", "--template", default="GumpleTemplate.xml.j2", help="Input Jinja2 template")
    parser.add_argument("-o", "--output", default="run_config.xml", help="Output XML filename")
    parser.add_argument("-d", "--dir", required=True, help="Absolute path to storage ROOT directory")
    parser.add_argument("--hdf-dir", default="/exp/sbnd/data/users/gputnam/GUMPLE/sbn-rewgted-22/", help="Directory containing HDF files")

    args = parser.parse_args()
    render_template(args.template, args.output, args.dir, args.hdf_dir)
