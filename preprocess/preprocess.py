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
        subprocess.run([self.script, input_file, output_file])
