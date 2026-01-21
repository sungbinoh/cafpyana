import subprocess

# Interface for a class to pre-process CAF inputs before running dataframe makers
class PreProcessor(object):
    def __init__(self):
        pass

    def run(self, input_file, output_file):
        pass


# pre-process CAF input by running a bash script
class Script(PreProcessor):
    def __init__(self, script):
        self.script = script

    def run(self, input_file, output_file):
        result = subprocess.run([self.script, input_file, output_file], 
                              capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(
                f"Preprocess script failed with return code {result.returncode}.\n"
                f"Script: {self.script}\n"
                f"Input: {input_file}\n"
                f"Output: {output_file}\n"
                f"stderr: {result.stderr}\n"
                f"stdout: {result.stdout}"
            )