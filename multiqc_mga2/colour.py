#!/usr/bin/env python

""" MultiQC Multi Genome Alignment 2 module """

class Colour(object):
    '''
    A simple class to represent RGB colour values, 0-255 on each primary colour.
    '''

    @staticmethod
    def fromBytes(r: int, g: int, b: int) -> 'Colour':
        return Colour(r, g, b)


    def __init__(self, r: int, g: int, b: int):
        self.red = r
        self.green = g
        self.blue = b


    def applyAlpha(self, alpha: float) -> 'Colour':
        '''
        Apply an alpha change to the current colour based on a white background.

        :param float alpha: The alpha value to apply, 0 <= alpha <= 1.

        :return A new Colour object for the altered colour.
        :rtype Colour
        '''
        return Colour(self._alphaValue(self.red, alpha),
                      self._alphaValue(self.green, alpha),
                      self._alphaValue(self.blue, alpha))


    def toHtml(self) -> str:
        '''
        Get an HTML representation of this colour.

        :return An HTML string for this colour.
        :rtype str
        '''
        return "#{:02x}{:02x}{:02x}".format(self.red, self.green, self.blue)


    def _alphaValue(self, c: int, a: float) -> int:
        '''
        Apply an alpha to a single primary colour, based on a white background.

        :param int c: The original primary colour.
        :param float a: The alpha to apply.

        :return The colour with the alpha applied.
        :rtype int
        '''
        # See https://stackoverflow.com/a/746937
        # 255 because we're always considering a white background.
        bgc = 255
        return int(a * c + (1.0 - a) * bgc)
