#!/usr/bin/python3

import os
import sys
import argparse
from src.pipeline import Pipeline


class Annotator:
    def __init__(self):
        self.args = self.__get_args()

    def __get_args(self):
        self.parser = argparse.ArgumentParser(
            description="Annotator - Platform for annotating smORF types",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        # Mode selection
        self.mode_parser = self.parser.add_argument_group("Mode input options")
        self.mode_parser.add_argument(
            "mode",
            metavar="Mode",
            help="Mode to run the pipeline for.\nList of Modes: smorf_types"
        )

        # General parameters
        self.general_args = self.parser.add_argument_group("General Parameters")
        self.general_args.add_argument("--outdir", "-o", help="Inform the output directory", default="Annotator_output")
        self.general_args.add_argument("--threads", "-p", help="Number of threads to be used.", type=int, default=1)

        # Add smorf_types mode arguments (since it's the only supported mode)
        self.__set_annotator_mode()

        # Parse all arguments at once
        args = self.parser.parse_args()
        
        # Validate mode after parsing
        self.mode = args.mode
        supported_modes = ["smorf_types"]
        if self.mode not in supported_modes:
            self.parser.error(f"Unsupported mode '{self.mode}'. Supported modes: {', '.join(supported_modes)}")

        # Ensure output directory exists
        os.makedirs(args.outdir, exist_ok=True)

        return args

        # Ensure output directory exists
        os.makedirs(args.outdir, exist_ok=True)

        return args

    def __set_annotator_mode(self):
        self.modeArguments = self.parser.add_argument_group("smORF annotation mode options")
        self.modeArguments.add_argument("--smorf_gtf", help="Provide the smORF GTF file", required=True)
        self.modeArguments.add_argument("--ensembl_gtf", help="Provide the ENSEMBL GTF file", required=True)
        self.modeArguments.add_argument("--intersect_output", help="Provide the intersect output file", default="Annotator_output/intersect.gtf")
        self.modeArguments.add_argument("--non_intersect_output", help="Provide the non-intersect output file", default="Annotator_output/nonintersect.gtf")
        self.modeArguments.add_argument("--output_file", help="Provide the output file", default="Annotator_output/smORF_annotation.txt")


    def execute(self):
        if self.mode == 'smorf_types':
            pipeline = Pipeline(args=self.args)
            pipeline.annotate()

if __name__ == '__main__':
    print("""
 ▗▄▖ ▗▖  ▗▖▗▖  ▗▖ ▗▄▖▗▄▄▄▖▗▄▖▗▄▄▄▖▗▄▖ ▗▄▄▖ 
▐▌ ▐▌▐▛▚▖▐▌▐▛▚▖▐▌▐▌ ▐▌ █ ▐▌ ▐▌ █ ▐▌ ▐▌▐▌ ▐▌
▐▛▀▜▌▐▌ ▝▜▌▐▌ ▝▜▌▐▌ ▐▌ █ ▐▛▀▜▌ █ ▐▌ ▐▌▐▛▀▚▖
▐▌ ▐▌▐▌  ▐▌▐▌  ▐▌▝▚▄▞▘ █ ▐▌ ▐▌ █ ▝▚▄▞▘▐▌ ▐▌
                                           
                                           
                                           """)

    print("""Annotator is a platform that annotates smORF types.""")
    Annotator = Annotator()
    Annotator.execute()
