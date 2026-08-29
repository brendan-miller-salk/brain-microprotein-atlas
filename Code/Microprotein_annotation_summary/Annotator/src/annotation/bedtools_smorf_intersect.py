import subprocess

class BedtoolsRunner:
    def __init__(self, args):
        self.smorf_gtf = args.smorf_gtf
        self.ensembl_gtf = args.ensembl_gtf
        self.intersect_file = args.intersect_output
        self.non_intersect_file = args.non_intersect_output

    def run_intersect(self):
        # Define the command for intersect
        intersect_command = [
            'bedtools', 'intersect', '-wo', '-s', '-a', self.smorf_gtf, '-b', self.ensembl_gtf
        ]

        # Execute the intersect command and redirect the output to a file
        with open(self.intersect_file, 'w') as outfile:
            subprocess.run(intersect_command, stdout=outfile, check=True)
            
        print(f"Intersect output file '{self.intersect_file}' created successfully.")

    def run_non_intersect(self):
        # Define the command for non-intersect
        non_intersect_command = [
            'bedtools', 'intersect', '-v', '-s', '-a', self.smorf_gtf, '-b', self.ensembl_gtf
        ]

        # Execute the non-intersect command and redirect the output to a file
        with open(self.non_intersect_file, 'w') as outfile:
            subprocess.run(non_intersect_command, stdout=outfile, check=True)

        print(f"Non-intersect output file '{self.non_intersect_file}' created successfully.")
