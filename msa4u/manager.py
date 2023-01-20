"""
This module provides managing classes and methods for the tool.
"""
import subprocess
import traceback
import argparse
import tempfile
import configs
import typing
import time
import sys
import os

import Bio.Seq
import Bio.SeqIO
import Bio.Align
import Bio.AlignIO
import Bio.Align.AlignInfo
import Bio.Data.CodonTable

from reportlab.lib.units import cm, mm
import reportlab.pdfbase.pdfmetrics
import reportlab.pdfbase.ttfonts
import reportlab.pdfgen.canvas
import reportlab.rl_config
import reportlab.pdfgen

import msa4u.methods
import msa4u.drawing

reportlab.rl_config.warnOnMissingFontGlyphs = 0


class MSA4uError(Exception):
    """A helper for exceptions parsing inherited from the Exception class.

    """
    pass


class Parameters:
    """A Parameters object holds and parse cmd's and config's arguments for the tool.

    Note:
        A Parameters object have to be created in each script since it's used by each
            class of the tool as a mandatory argument.

    """

    def __init__(self, api=True, config="standard"):
        self.arguments = dict(debug=False, verbose=False)
        self.cmd_arguments = {"output_filename_aln": "auto", "output_filename": "auto", "verbose": True,
                              "debug": False, "sequence_type": "auto", "alignments": None, "fasta": None}
        if api:
            self.load_config(config)

    def parse_cmd_arguments(self) -> None:
        """Parse command-line arguments

        Returns:
            None

        """
        parser = argparse.ArgumentParser(prog="msa4u", add_help=False,
                                         usage="msa4u [-fa path | -aln path] [optional arguments]")
        mutually_exclusive_group = parser.add_mutually_exclusive_group()
        mutually_exclusive_group.add_argument("-fa", dest="fasta", type=str, default=None)
        mutually_exclusive_group.add_argument("-aln", dest="alignments", type=str, default=None)
        parser.add_argument("-data", "--data", dest="msa4u_data", action="store_true")
        parser.add_argument("-linux", "--linux", dest="linux", action="store_true", default=None)
        parser.add_argument("-label", dest="label", choices=["id", "description", "all"], default="all")
        parser.add_argument("-o-aln", dest="output_filename_aln", type=str, default="auto")
        parser.add_argument("-o", dest="output_filename", type=str, default="auto")
        parser.add_argument("-st", dest="sequence_type", choices=["nt", "aa", "auto"], type=str, default="auto")
        parser.add_argument("-c", dest="config_file", type=str, default="standard")
        parser.add_argument("-v", "--version", action='version', version='%(prog)s 0.4.0')
        parser.add_argument("-q", "--quiet", dest="verbose", default=True, action="store_false")
        parser.add_argument("--debug", "-debug", dest="debug", action="store_true")
        parser.add_argument("-h", "--help", dest="help", action="store_true")
        args = parser.parse_args()
        args = vars(args)

        if len(sys.argv[1:]) == 0:
            args["help"] = True

        if args["msa4u_data"]:
            msa4u.methods.copy_package_data()
            sys.exit()

        if args["linux"]:
            msa4u.methods.adjust_paths_for_linux()
            sys.exit()

        if args["help"]:
            help_message_path = os.path.join(os.path.dirname(__file__), 'msa4u_data', "help.txt")
            with open(help_message_path, "r") as help_message:
                print(help_message.read(), file=sys.stdout)
                sys.exit()

        filtered_args = {k: v for k, v in args.items() if (v is not None or (k == "alignments" or k == "fasta"))}

        self.cmd_arguments = filtered_args
        return None

    def load_config(self, path: str = "standard") -> None:
        """Load configuration file parameters.

        Arguments
            path (str): path to a config file or name (only standard available at this moment)

        Returns:
            None

        """
        try:
            if path == "standard":
                path = os.path.join(os.path.dirname(__file__), "msa4u_data", f"standard.cfg")
            config = configs.load(path)
            config = config.get_config()
            internal_dir = os.path.dirname(__file__)
            for key in config["root"].keys():
                if type(config["root"][key]) is str and "{internal}" in config["root"][key]:
                    config["root"][key] = config["root"][key].replace("{internal}",
                                                                      os.path.join(internal_dir, "msa4u_data"))
            self.arguments.update(config["root"])
            self.arguments.update(self.cmd_arguments)
            self.load_palette()
            self.load_color_config()
            return None
        except Exception as error:
            raise MSA4uError(
                "Unable to parse the specified config file. Please check your config file or written name.") from error

    def load_palette(self) -> None:
        """Load palette file.

        Returns:
            None

        """
        palette_path = self.arguments[f"palette"]
        self.arguments[f"palette"] = configs.load(palette_path).get_config()["root"]

    def load_color_config(self) -> None:
        """Load color config files.

        Returns:
            None

        """
        for seq_type in ["nt", "aa"]:
            path = self.arguments[f"colors_{seq_type}"]
            colors_pre_dict = configs.load(path).get_config()["root"]
            colors_dict = dict()
            for elements, color in colors_pre_dict.items():
                for element in elements:
                    colors_dict[element] = color
            self.arguments[f"colors_{seq_type}"] = colors_dict

    def update(self, parameters):
        self.arguments.update(parameters)


class Fasta:
    """Fasta class objects holds and manages a fasta file with unaligned sequences.

    Attributes:

    """

    def __init__(self, fasta, parameters):
        """Create a Fasta class object.

        Arguments:
            fasta (str): path to a fasta file.
            parameters (Parameters): Parameters' class object.

        """
        self.fasta = fasta
        self.parameters = parameters
        if not self.parameters.arguments["fasta"]:
            self.parameters.arguments["fasta"] = self.fasta

    def run_mafft(self) -> Bio.Align.MultipleSeqAlignment:
        """Run MAFFT tool to get MSA.

        Returns:
            Bio.Align.MultipleSeqAlignment: MSA
        """
        try:
            if self.parameters.arguments["save_msa_results"]:
                output_name = self.parameters.arguments["output_filename_aln"]
                if output_name == "auto":
                    output_name = msa4u.methods.update_path_extension(self.fasta, "aln.fa")
                mafft_output = open(output_name, "w")
                if self.parameters.arguments["verbose"]:
                    print(f"💌 Alignments file was saved as: {os.path.basename(output_name)}", file=sys.stdout)
            else:
                mafft_output = tempfile.NamedTemporaryFile()
            mafft = self.parameters.arguments["mafft_binary"]
            if self.parameters.arguments["mafft_reorder_parameter"]:
                subprocess.run([mafft, "--auto", "--reorder", self.fasta], stdout=mafft_output,
                               stderr=subprocess.DEVNULL)
            else:
                subprocess.run([mafft, "--auto", self.fasta], stdout=mafft_output, stderr=subprocess.DEVNULL)
            msa = Bio.AlignIO.read(mafft_output.name, "fasta")
            mafft_output.close()
            for record in msa:
                record.description = " ".join(record.description.split(" ")[1:])
            return msa
        except Exception as error:
            raise MSA4uError("Unable to run MAFFT.") from error


class MSA:
    """MSA object holds and manages multiple sequence alignments.

    Attributes:
        parameters (Parameters): Parameters' class object.
        msa (Bio.Align.MultipleSeqAlignment): Biopython' MSA object.

    """

    def __init__(self, msa: typing.Union[str, Bio.Align.MultipleSeqAlignment], parameters: Parameters) -> None:
        """Create a MSA object.

        Arguments:
            msa (str or Bio.Align.MultipleSeqAlignment): Path to a fasta file with MSA or
                Bio.Align.MultipleSeqAlignment object.
            parameters (Parameters): Parameters' class object.

        """
        try:
            self.parameters = parameters
            if isinstance(msa, str):
                self.msa = Bio.AlignIO.read(msa, "fasta")
                for record in self.msa:
                    record.description = " ".join(record.description.split(" ")[1:])
            if isinstance(msa, Bio.Align.MultipleSeqAlignment):
                self.msa = msa
            for record in self.msa:
                if self.parameters.arguments["label"] == "id":
                    record.annotations["label"] = record.id
                elif self.parameters.arguments["label"] == "description":
                    record.annotations["label"] = record.description
                elif self.parameters.arguments["label"] == "all":
                    record.annotations["label"] = f"{record.id} {record.description}"
                else:
                    raise MSA4uError("Incorrect label parameter (can be id, desription or all)\n"
                                     f"Check your config file.")
            if parameters.arguments["sequence_type"] == "auto":
                used_alphabet = set("".join([str(record.seq) for record in self.msa]))
                ambiguous_codon_table = Bio.Data.CodonTable.ambiguous_dna_by_name["Standard"]
                ambiguous_alphabet = dict(nt=set(ambiguous_codon_table.nucleotide_alphabet),
                                          aa=set(ambiguous_codon_table.protein_alphabet))
                if len((used_alphabet - ambiguous_alphabet["nt"]) & ambiguous_alphabet["aa"]):
                    parameters.arguments["sequence_type"] = "aa"
                else:
                    parameters.arguments["sequence_type"] = "nt"
        except Exception as error:
            raise MSA4uError("Unable to create an MSA object.") from error

    def plot(self) -> None:
        """Run MSA visualisation.

        Returns: None

        """
        try:
            msa_plot_manager = MSAPlotManager(self.msa, self.parameters)
            msa_plot_manager.define_x_axis_coordinate_system()
            msa_plot_manager.create_tracks()
            if self.parameters.arguments["output_filename"] == "auto":
                if self.parameters.arguments["alignments"]:
                    self.parameters.arguments["output_filename"] = msa4u.methods.update_path_extension(
                        self.parameters.arguments["alignments"], "pdf")
                elif self.parameters.arguments["fasta"]:
                    self.parameters.arguments["output_filename"] = msa4u.methods.update_path_extension(
                        self.parameters.arguments["fasta"], "pdf")
                else:
                    self.parameters.arguments["output_filename"] = f"msa4u_{time.strftime('%Y_%m_%d-%H_%M')}.pdf"
            msa_plot_manager.plot(self.parameters.arguments["output_filename"])
            return None
        except Exception as error:
            raise MSA4uError("Unable to produce the MSA plot.") from error


class MSAPlotManager:
    """
    AnnotationPlotManager object holds needed information for MSA visualisation and controls it.

    Attributes:
        msa (FILL IN): Multiple sequence alignment.
        parameters (Parameters): Parameters' class object.
        coordinate_system (dict): coordinate system of figure.
        additional_data (dict): dict with data for visualisation tracks.

    """

    def __init__(self, msa, parameters: Parameters):
        """Create a AnnotationPlotManager object.

        Arguments:
            msa (FILL IN): Multiple sequence alignment.
            parameters (Parameters): Parameters' class object.

        """
        self.msa = msa
        self.parameters = parameters
        self.coordinate_system = dict()
        self.additional_data = dict()

    def define_x_axis_coordinate_system(self) -> None:
        """Define coordinate system.

        Returns:
            None

        """
        label_height = self.parameters.arguments["label_size"] * self.parameters.arguments["tile_size"] * cm
        label_font_size = msa4u.methods.string_height_to_font_size(label_height, "regular", self.parameters.arguments)
        self.additional_data["label_font_size"] = label_font_size
        msa_length = self.msa.get_alignment_length()
        max_label_width = max(
            [reportlab.pdfbase.pdfmetrics.stringWidth(i.annotations["label"], "regular",  # to regulate
                                                      label_font_size) for i in self.msa])

        char_height = self.parameters.arguments["char_size"] * self.parameters.arguments["tile_size"] * cm
        char_font_size = msa4u.methods.string_height_to_font_size(char_height, "mono", self.parameters.arguments)
        self.additional_data["label_size"] = self.parameters.arguments["label_size"]
        self.additional_data["char_font_size"] = char_font_size
        self.additional_data["char_size"] = self.parameters.arguments["char_size"]
        self.additional_data["tile_size"] = self.parameters.arguments["tile_size"] * cm
        self.additional_data["number_of_sequences"] = len(self.msa)
        self.coordinate_system["x_labels_start"] = self.parameters.arguments["margin"] * cm
        self.coordinate_system["x_labels_stop"] = self.coordinate_system["x_labels_start"] + max_label_width
        self.coordinate_system["x_msa_start"] = self.coordinate_system["x_labels_stop"] + \
                                                self.parameters.arguments["label_gap"] * cm
        msa_width = self.parameters.arguments["tile_size"] * msa_length * cm
        self.coordinate_system["x_msa_stop"] = self.coordinate_system["x_msa_start"] + msa_width
        self.coordinate_system["figure_width"] = 2 * self.parameters.arguments["margin"] * cm + msa_width + \
                                                 max_label_width + self.parameters.arguments["label_gap"] * cm
        self.coordinate_system["figure_height"] = self.parameters.arguments["margin"] * cm
        self.additional_data["palette"] = self.parameters.arguments[
            f"colors_{self.parameters.arguments['sequence_type']}"]
        self.additional_data["palette"] = {k: msa4u.methods.color_name_to_hex(v, self.parameters.arguments) for k, v in
                                           self.additional_data["palette"].items()}

        return None

    def create_tracks(self) -> None:
        """Create visualisation tracks.

        Returns:
            None

        """
        self.tracks = []
        """
        title_loader = TitleLoader(self.parameters)
        title_loader.prepare_data(self.coordinate_system, self.additional_data)
        title_track = title_loader.create_track()
        self.tracks.append(title_track)
        self.coordinate_system["figure_height"] += title_track.needed_y_space()
        """
        for record in self.msa:
            sequence_loader = SequencesLoader(self.parameters)
            sequence_loader.prepare_data(record, self.coordinate_system, self.additional_data)
            track = sequence_loader.create_track()
            self.tracks.append(track)
            self.coordinate_system["figure_height"] += track.needed_y_space()
            # self.coordinate_system["figure_height"] += self.parameters.arguments["gap"] * cm
            # if index < self.additional_data["number_of_sequences"] - 1:
        self.coordinate_system["figure_height"] += self.parameters.arguments["margin"] * cm

    def plot(self, filename):
        image = Image(filename, self.coordinate_system["figure_width"], self.coordinate_system["figure_height"])
        current_y_top = self.coordinate_system["figure_height"] - self.parameters.arguments["margin"] * cm
        for track in self.tracks:
            track.visualisation_data["y_top"] = current_y_top
            track.draw(image.canvas)
            current_y_top -= (track.needed_space)
        image.save()
        if self.parameters.arguments["verbose"]:
            print(f"🎨 MSA plot was saved as: {os.path.basename(filename)}", file=sys.stdout)
        return None


class Loader:
    """Parent class for tracks loaders.

    Attributes:
        parameters (Parameters): Parameters' class object.
        prepared_data (dict): dict with data needed for visualisation tracks.

    """

    def __init__(self, parameters: Parameters):
        """Parent's constructor for creating a Loader class object.

        Arguments:
            parameters (Parameters): Parameters' class object.
        """
        self.parameters = parameters
        self.prepared_data = None

    def prepare_data(self) -> None:
        """Empty parent's method for data preparation.

        Returns:
            None

        """
        pass

    def create_track(self) -> None:
        """Empty parent's method for initialisation of a track.

        Returns:
            None

        """
        pass


class SequencesLoader(Loader):
    """A SequencesLoader object prepares data for a Sequence track object.


    Attributes:
        parameters (Parameters): Parameters' class object.
        prepared_data (dict): dict with data needed for a visualisation track.

    """

    def __init__(self, parameters):
        """Create a SequenceLoader object.

        Arguments:
            parameters Parameters): Parameters' class object.

        """
        super().__init__(parameters)

    def prepare_data(self, record, coordinate_system: dict, additional_data: dict) -> dict:
        """Prepare data for a Title visualisation track.

        Attributes:
            record (FILL in): record of blablabla
            coordinate_system (dict): coordinate system of a figure page.
            additional_data (dict): data needed for a track initialisation.

        Returns:
            dict: dictionary with prepared data for visualisation.

        """
        prepared_data = dict()
        prepared_data["coordinate_system"] = coordinate_system
        prepared_data["label_size"] = additional_data["label_size"]
        prepared_data["char_size"] = additional_data["char_size"]
        prepared_data["tile_size"] = additional_data["tile_size"]
        prepared_data["label_font_size"] = additional_data["label_font_size"]
        prepared_data["char_font_size"] = additional_data["char_font_size"]
        prepared_data["sequence"] = record.seq
        prepared_data["label"] = record.annotations["label"]
        prepared_data["palette"] = additional_data["palette"]
        self.prepared_data = prepared_data

        return prepared_data

    def create_track(self) -> msa4u.drawing.SequenceVis:
        """Initialise a Sequence track object.

        Returns:
            SequenceVis: visualisation track.

        """
        return msa4u.drawing.SequenceVis(self.prepared_data, self.parameters)


class Image:
    """An Image object holds pdf.

    Attributes:
        canvas (reportlab.pdfgen.canvas.Canvas): pdf object of the reportlab library.

    """

    def __init__(self, filename: str, width: float, height: float):
        """Create an Image object.

        Arguments:
            filename (str): path and name of a pdf.
            width (float): width of a pdf.
            height (float): height of a pdf.

        """
        self.canvas = reportlab.pdfgen.canvas.Canvas(filename, pagesize=(width, height))

    def save(self) -> None:
        """Save a pdf file.

        Returns:
            None

        """
        self.canvas.save()
        return None
