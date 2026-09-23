import uproot
import awkward as ak
import numpy as np
import os
import pandas as pd

def export_dataframe_to_uproot(df, output_filename, tree_name="SelectedEvents"):
    """Dynamically steps through a Pandas DataFrame, cleans up union types,

    and writes the columns to a flat/jagged TTree.
    """
    uproot_payload = {}

    print(f"Processing {len(df.columns)} columns for TTree conversion...")

    for col in df.columns:
        if df[col].empty:
            continue

        sample_element = df[col].iloc[0]
        # Check for nested structures (lists, arrays, vectors)
        if isinstance(sample_element, (list, tuple, np.ndarray)):
            print(f"  -> Column '{col}': Detected nested list structure.")
            ak_arr = ak.Array(df[col].tolist())

            if "union" in str(ak.type(ak_arr)):
                print(f"     [!] Warning: Union type detected in '{col}'. Enforcing float64 conversion.")
                ak_arr = ak.values_astype(ak_arr, "float64")

            uproot_payload[col] = ak_arr

        else:
            print(f"  -> Column '{col}': Detected flat scalar structure.")
            if df[col].dtype == object:
                print(
                    f"     [!] Warning: Column '{col}' has object dtype. Forcing numeric conversion."
                )
                uproot_payload[col] = pd.to_numeric(df[col], errors="coerce").fillna(-999).to_numpy()
            elif df[col].dtypes != 'category' and df[col].dtypes != 'Int8':
                uproot_payload[col] = df[col].fillna(-999).to_numpy()
            elif df[col].dtypes == 'Int8':
                uproot_payload[col] = df[col].fillna(-128).to_numpy()

    # Write out the sanitized payload
    print(f"Writing staged columns to {output_filename} [{tree_name}]...")
    with uproot.recreate(output_filename) as f:
        f[tree_name] = uproot_payload

    print("Stage 1 write complete.")

def run_makesbruce_macro(input_file, output_file, macro_dir="/exp/sbnd/app/users/nrowe/cafpyana/analysis_village/PROfit/", root_path="/cvmfs/larsoft.opensciencegrid.org/spack-fnal-v1.0.0/opt/spack/linux-x86_64_v2/root-6.28.12-vgs6mjswsg36hl3oarsrsyc2dcua6khe/bin/thisroot.sh"):
    """
    Invokes the MakesBruce.C ROOT macro from Python using the os package.
    """

    # Path to the actual macro file
    macro_path = os.path.join(macro_dir, "MakesBruceNew.C")

    # Construct the exact command string.
    # We escape the internal quotes so ROOT parses the file strings properly.
    command = f'source {root_path} && root -q \'{macro_path}("{input_file}", "{output_file}")\''

    print(f"Executing: {command}")

    # Run the command via the shell
    exit_status = os.system(command)

    # os.system returns the shell exit status (0 means success)
    if exit_status != 0:
        print(
            f"Error: ROOT macro failed with shell exit code {exit_status}."
        )
    else:
        print("ROOT macro executed successfully.")

    return exit_status
