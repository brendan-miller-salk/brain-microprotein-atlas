import os
import shutil
from ..annotation import BedtoolsRunner, smORFAnnotator

class Pipeline:
    def __init__(self, args):
        self.args = args
        self.outdir = args.outdir
    
    def annotate(self):
        print("▶️ You have successfully initiated smORF annotation...")
        run = BedtoolsRunner(args=self.args)
        run.run_intersect()
        run.run_non_intersect()

        annotate = smORFAnnotator(args=self.args)
        annotate.process_gtf_files()

    def __cleanup_output_directory(self, outdir):
        """Removes all directories and files in the specified output directory and recreates the directory."""
        try:
            # Remove the output directory and all its contents
            shutil.rmtree(outdir)
            print(f"🧹 Output directory '{outdir}' cleaned up successfully.")
        except Exception as e:
            print(f"Error cleaning up output directory '{outdir}': {e}")

        