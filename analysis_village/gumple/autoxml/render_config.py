import os
import numpy as np
import sys
import glob
import argparse
from jinja2 import Environment, FileSystemLoader

# Resolve and add ../../gump to sys.path
script_dir = os.path.dirname(os.path.abspath(__file__))
gump_dir = os.path.abspath(os.path.join(script_dir, "../../gump"))

if gump_dir not in sys.path:
    sys.path.insert(0, gump_dir)

import loaddf as loaddf

def format_sci(value, precision=3):
    """Formats floats in scientific notation without the '+' in the exponent (e.g. 5.272e20)."""
    return f"{value:.{precision}e}".replace("e+", "e")

def get_sample_pot(file_pattern, detector, offbeampot=False, data_quality=False, beam_quality=False):
    """
    Finds all files matching pattern and returns the total POT via loaddf.loadl,
    including dedup correction. For offbeam samples pass offbeampot=True.
    For data samples that need good-run / beam-quality filtering pass the
    corresponding flags.
    """
    files = glob.glob(file_pattern)
    if not files:
        print(f"Warning: No files found matching pattern: {file_pattern}")
        return 0.0

    _, _, total_pot = loaddf.loadl(
        files,
        detector=detector,
        load_truth=False,
        load_evtrec=False,
        include_syst=False,
        reweight_aFF=False,
        match_Enu=True,
        offbeampot=offbeampot,
        data_quality=data_quality,
        beam_quality=beam_quality,
    )
    print(f"  {file_pattern}  ->  POT = {total_pot:.4e}")
    return total_pot

def get_scale(sample_pot, target_pot):
    """Calculates scaling factor: scale = target_pot / sample_pot"""
    if sample_pot == 0:
        return 1.0
    return target_pot / sample_pot

def build_context(base_dir, hdf_dir, target_sbnd_pot=1e20, target_icarus_r2_pot=2e20, target_icarus_r4_pot=3e20):
    """
    Reads every sample's POT once and returns the Jinja2 context. The result
    depends only on base_dir/hdf_dir, so it can be reused for any number of
    templates without re-reading the HDF files.
    """
    print("--- Calculating POT and Scale Factors ---")

    # 1. SBND MC (Wildcard)
    sbnd_mc_pattern = os.path.join(hdf_dir, "SBNDMCCV_*.df")
    sbnd_mc_pot = get_sample_pot(sbnd_mc_pattern, detector="SBND")

    # 6. SBND OffBeam (Data)
    sbnd_offbeam_pattern = os.path.join(hdf_dir, "SBND_SpringBNBOffData.df")
    sbnd_offbeam_pot = get_sample_pot(sbnd_offbeam_pattern, detector="SBND",
                                      offbeampot=True, data_quality=True)

    # 9. SBND Dirt (MC)
    sbnd_dirt_pattern = os.path.join(hdf_dir, "SBND_SpringLowEMC.df")
    sbnd_dirt_pot = get_sample_pot(sbnd_dirt_pattern, detector="SBND")

    # SBND OnBeam (Data)
    sbnd_onbeam_pattern = os.path.join(hdf_dir, "SBND_SpringBNBData_FixedDev.df")
    sbnd_onbeam_pot = get_sample_pot(sbnd_onbeam_pattern, detector="SBND",
                                     data_quality=True, beam_quality=True)

    # 3. ICARUS Run 2 MC (Wildcard)
    icarus_r2_mc_pattern = os.path.join(hdf_dir, "ICARUSRun2_SpringMCOverlay_rewgt_*.df")
    icarus_r2_mc_pot = get_sample_pot(icarus_r2_mc_pattern, detector="ICARUS Run2")

    # 4. ICARUS Run 2 OffBeam (Data)
    icarus_r2_offbeam_pattern = os.path.join(hdf_dir, "ICARUS_SpringRun2BNBOff_unblind.df")
    icarus_r2_offbeam_pot = get_sample_pot(icarus_r2_offbeam_pattern, detector="ICARUS Run2",
                                           offbeampot=True, data_quality=True)

    # 7. ICARUS Run 2 Dirt (MC)
    icarus_r2_dirt_pattern = os.path.join(hdf_dir, "ICARUSRun2_Spring_Overlay_Dirt.df")
    icarus_r2_dirt_pot = get_sample_pot(icarus_r2_dirt_pattern, detector="ICARUS Run2")

    # ICARUS Run 2 OnBeam (Data)
    icarus_r2_onbeam_pattern = os.path.join(hdf_dir, "ICARUS_SpringRun2BNB_unblind.df")
    icarus_r2_onbeam_pot = get_sample_pot(icarus_r2_onbeam_pattern, detector="ICARUS Run2",
                                          data_quality=True, beam_quality=True)

    # 2. ICARUS Run 4 MC (Wildcard)
    icarus_r4_mc_pattern = os.path.join(hdf_dir, "ICARUSRun4_SpringMCOverlay_rewgt_*.df")
    icarus_r4_mc_pot = get_sample_pot(icarus_r4_mc_pattern, detector="ICARUS Run4")

    # 5. ICARUS Run 4 OffBeam (Data)
    icarus_r4_offbeam_pattern = os.path.join(hdf_dir, "ICARUS_SpringRun4BNBOff_ReCAF2026.df")
    icarus_r4_offbeam_pot = get_sample_pot(icarus_r4_offbeam_pattern, detector="ICARUS Run4",
                                           offbeampot=True, data_quality=True)

    # 8. ICARUS Run 4 Dirt (MC)
    icarus_r4_dirt_pattern = os.path.join(hdf_dir, "ICARUSRun4_Spring_Overlay_Dirt.df")
    icarus_r4_dirt_pot = get_sample_pot(icarus_r4_dirt_pattern, detector="ICARUS Run4")

    # ICARUS Run 4 OnBeam (Data)
    icarus_r4_onbeam_pattern = os.path.join(hdf_dir, "ICARUS_SpringRun4BNB_unblind.df")
    icarus_r4_onbeam_pot = get_sample_pot(icarus_r4_onbeam_pattern, detector="ICARUS Run4",
                                          data_quality=True, beam_quality=True)

    context = {
        "BASE_DIR": base_dir,
        "sbnd_target_pot": format_sci(target_sbnd_pot, 1),
        "icarus_target_pot": format_sci(target_icarus_r2_pot+target_icarus_r4_pot, 1),
        
        # SBND POT values
        "sbnd_mc_pot": format_sci(sbnd_mc_pot, 3),
        "sbnd_offbeam_pot": format_sci(sbnd_offbeam_pot, 3),
        "sbnd_dirt_pot": format_sci(sbnd_dirt_pot, 3),

        # On-beam data POT (the target POT each detector is normalised to)
        "sbnd_data_pot": format_sci(sbnd_onbeam_pot, 3),
        "icarus_r2_data_pot": format_sci(icarus_r2_onbeam_pot, 3),
        "icarus_r4_data_pot": format_sci(icarus_r4_onbeam_pot, 3),

        # ICARUS Run 2 sample POTs
        "icarus_r2_mc_pot": format_sci(icarus_r2_mc_pot, 3),
        "icarus_r2_offbeam_pot": format_sci(icarus_r2_offbeam_pot, 3),
        "icarus_r2_dirt_pot": format_sci(icarus_r2_dirt_pot, 3),

        # ICARUS Run 4 sample POTs
        "icarus_r4_mc_pot": format_sci(icarus_r4_mc_pot, 3),
        "icarus_r4_offbeam_pot": format_sci(icarus_r4_offbeam_pot, 3),
        "icarus_r4_dirt_pot": format_sci(icarus_r4_dirt_pot, 3),
    }

    return context

def render_templates(template_paths, output_paths, context, base_dirs=None):
    """
    Renders each template with an already-built context. base_dirs, if given,
    overrides BASE_DIR per template, so templates pointing at different storage
    roots can share a single POT read.
    """
    print("\n--- Rendering Jinja2 Configuration ---")

    for i, (template_path, output_path) in enumerate(zip(template_paths, output_paths)):
        template_context = dict(context)
        if base_dirs:
            template_context["BASE_DIR"] = base_dirs[i]

        env = Environment(
            loader=FileSystemLoader(os.path.dirname(template_path) or '.'),
            comment_start_string='{%#',
            comment_end_string='#%}'
        )
        template = env.get_template(os.path.basename(template_path))

        rendered_content = template.render(**template_context)

        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(rendered_content)

        print(f"Successfully written configuration to {output_path}")

def render_template(template_path, output_path, base_dir, hdf_dir, **kwargs):
    """Single-template entry point, kept for existing callers."""
    context = build_context(base_dir, hdf_dir, **kwargs)
    render_templates([template_path], [output_path], context)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Render XML config with dataset POT calculations.")
    parser.add_argument("-t", "--template", nargs="+", default=["GumpleTemplate.xml.j2"],
                        help="Input Jinja2 template(s). POT is read once and reused for all of them.")
    parser.add_argument("-o", "--output", nargs="+",
                        help="Output XML filename(s), one per template. Defaults to each template "
                             "with its .j2 suffix stripped, written into --outdir.")
    parser.add_argument("--outdir", default=".", help="Directory for derived output names when -o is omitted")
    parser.add_argument("-d", "--dir", nargs="+", required=True,
                        help="Absolute path to storage ROOT directory. Give one value for all "
                             "templates, or one per template to render several roots in one pass.")
    parser.add_argument("--hdf-dir", default="/exp/sbnd/data/users/gputnam/GUMPLE/sbn-rewgted-24/", help="Directory containing HDF files")

    args = parser.parse_args()

    if args.output:
        if len(args.output) != len(args.template):
            parser.error(f"got {len(args.template)} template(s) but {len(args.output)} output name(s); "
                         "pass one -o per -t, or omit -o to derive the names")
        outputs = args.output
    else:
        outputs = []
        for t in args.template:
            name = os.path.basename(t)
            if name.endswith(".j2"):
                name = name[:-len(".j2")]
            outputs.append(os.path.join(args.outdir, name))

    if len(args.dir) == 1:
        base_dirs = args.dir * len(args.template)
    elif len(args.dir) == len(args.template):
        base_dirs = args.dir
    else:
        parser.error(f"got {len(args.template)} template(s) but {len(args.dir)} -d value(s); "
                     "pass a single -d for all of them, or one per -t")

    context = build_context(base_dirs[0], args.hdf_dir)
    render_templates(args.template, outputs, context, base_dirs=base_dirs)
