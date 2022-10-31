"""
This module provides visualisation of loci annotation.
"""
from reportlab.lib.units import cm, mm
import reportlab.pdfbase.pdfmetrics
import reportlab.pdfbase.ttfonts
import reportlab.pdfgen.canvas
import reportlab.rl_config
import reportlab.pdfgen

reportlab.rl_config.warnOnMissingFontGlyphs = 0

import msa4u.methods
import msa4u.manager


class Track:
    """Parent clas for visualisation Tracks.

    Attributes:
        visualisation_data (dict): a dictionary with data needed for visualisation.
        parameters (uorf4u.manager.Parameters): Parameters' class object.

    """

    def __init__(self, visualisation_data: dict, parameters):
        """Parent's constructor for creating a Track object.

        Arguments:
            visualisation_data (dict): a dictionary with data needed for visualisation.

        """
        self.visualisation_data = visualisation_data
        self.parameters = parameters

    def needed_y_space(self) -> None:
        """Empy parent's method for calculation needed vertical space for a track.

        Returns:
            None

        """
        pass

    def draw(self, canvas: reportlab.pdfgen.canvas.Canvas) -> None:
        """Empy parent's method for track visualisation.

        Arguments:
            canvas (reportlab.pdfgen.canvas.Canvas): a pdf object.

        Returns:
            None
        """
        pass


class SequenceVis(Track):
    """SequenceVis track draws sequences and annotation.

    Attributes:
        visualisation_data (dict): a dictionary with data needed for visualisation.
        parameters (uorf4u.manager.Parameters): Parameters' class object.
        needed_space (float): needed vertical space for a track.

    """

    def __init__(self, visualisation_data: dict, parameters):
        """Create a SequenceVis object.

        Arguments:
            visualisation_data (dict): a dictionary with data needed for visualisation.

        """
        super().__init__(visualisation_data, parameters)
        self.needed_space = None

    def needed_y_space(self) -> float:
        """Calculate needed vertical space for a SequenceVis track.

        Returns:
            float: needed vertical space.

        """
        self.needed_space = self.visualisation_data["tile_size"]
        return self.needed_space

    def draw(self, canvas: reportlab.pdfgen.canvas.Canvas) -> None:
        """Draw a Sequence track.

        Arguments:
            canvas (reportlab.pdfgen.canvas.Canvas): a pdf object.

        Returns:
            None

        """
        tile_size = self.visualisation_data["tile_size"]
        y_c = self.visualisation_data["y_top"] - (tile_size * 0.5)
        y_l = self.visualisation_data["y_top"] - tile_size

        y_gap_label = tile_size * (1 - self.visualisation_data["label_size"]) * 0.5
        y_gap_char = tile_size * (1 - self.visualisation_data["char_size"]) * 0.5
        # Labels
        canvas.setFillColorRGB(*msa4u.methods.get_color("label_color", self.parameters.arguments))
        canvas.setFont("regular", self.visualisation_data["label_font_size"])
        canvas.drawRightString(self.visualisation_data["coordinate_system"]["x_labels_stop"], y_l + y_gap_label,
                               self.visualisation_data["label"])

        canvas.setLineWidth(0.05 * tile_size)
        canvas.setStrokeColorRGB(1, 1, 1)
        canvas.setFont("mono", self.visualisation_data["char_font_size"])
        x = self.visualisation_data["coordinate_system"]["x_msa_start"]
        for symbol in self.visualisation_data["sequence"]:
            x_c = x + tile_size * 0.5
            symbol = symbol.upper()
            try:
                color = self.visualisation_data["palette"][symbol]
            except:
                color = "#FFFFFF"
            canvas.setFillColorRGB(*msa4u.methods.hex_to_rgb(color), self.parameters.arguments["tile_alpha"])
            canvas.rect(x, y_l, tile_size, tile_size, fill=1, stroke=self.parameters.arguments["stroke"])
            canvas.setFillColorRGB(0, 0, 0, 0.8)  # to change
            canvas.drawCentredString(x_c, y_l + y_gap_char, symbol)
            x += tile_size
        return None
