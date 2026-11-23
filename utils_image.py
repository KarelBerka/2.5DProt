#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import math
import svgwrite
import utils_sses_types
import utils_vector

class Image(object):
    """Just a wrapper around svgwrite"""

    ZOOM = 20
    DEFAULT_TEXT_SIZE = 1
    CIRCLE_MULTIPLIER = 0.1 # have to be 1/2 from the rel number #.20#  05
    HELIC_TEXT_MULTIPLIER = 0.04

    def __init__(self, filename, dimension, fill="white", center=None):
        if fill == True:
            fill = "white"
        self.dimension = dimension
        if center is None:
            center = (utils_vector.frac(dimension, 2))
        self.center = center

        # Create SVG drawing instance
        self.dwg = svgwrite.Drawing(filename,
            size=(self.dimension[0] * self.ZOOM, self.dimension[1] * self.ZOOM))

        # leave the background transparent
        # Create white background
        if fill:
            background = self.dwg.rect(insert=(0, 0), size=('100%', '100%'), rx=None, ry=None, fill=fill)
            self.dwg.add(background)

    def get_upper_right(self):
        """Return coords of upper right corner of a canvas"""
        return (self.dimension[0] - self.center[0], self.dimension[1] - self.center[1])

    def get_upper_left(self):
        """Return coords of upper left corner of a canvas"""
        return (-1 * self.center[0], self.dimension[1] - self.center[1])

    def get_lower_left(self):
        """Return coords of upper left corner of a canvas"""
        return (-1 * self.center[0], -1 *  self.center[1])

    def add_group(self, group_name, color="red", opacity=1, classname="", ggroup = False, without_id=False):
        if ggroup == False:
            if without_id:
                gr = self.dwg.add(self.dwg.g(class_=classname, style='stroke: %s; fill: %s; opacity:%s' \
                    % (color, color, opacity)))
            else:
                gr = self.dwg.add(self.dwg.g(id=group_name, class_=classname, style='stroke: %s; fill: %s; opacity:%s' \
                    % (color, color, opacity)))
        else:
            if without_id:
                gr = ggroup.add(self.dwg.g(class_=classname, style='stroke: %s; fill: %s; opacity:%s' \
                    % (color, color, opacity)))
            else:
                gr = ggroup.add(self.dwg.g(id=group_name, class_=classname, style='stroke:%s; fill: %s; opacity:%s' \
                    % (color, color, opacity)))
        return gr

    def save(self):
        self.dwg.save()
        return self.dwg.filename

    # Pure geometric shape methods:
    def _transform(self, coord):
        coord = utils_vector.add(coord, self.center)   # move coordinates center to center of a drawing area
        coord = (coord[0], self.dimension[1] - coord[1])   # our coordinates y goes from bottom to up, but in svg it is from up to bottom
        coord = utils_vector.mult(coord, self.ZOOM)   # apply zoom
        coord = [int(coord[0]),int(coord[1])]
        return coord

    def _line(self, start, end, **args):
        return self.dwg.line(self._transform(start), self._transform(end), **args)

    def line(self, start, end, group=False,  **args):
        line = self._line(start, end, **args)
        if group == False:
            return self.dwg.add(line)
        else:
            return group.add(line)

    def polygon(self, coord, group=False, **args):
        ###print ">>> line() from %s to %s with %s" % (start, end, args)
        pol = self.dwg.polygon([self._transform(coord[0]), self._transform(coord[1]), self._transform(coord[2]), self._transform(coord[3])], **args)
        if group:
            return group.add(pol)
        else:
            return self.dwg.add(pol)

    def rect(self, center, size, rx=None, ry=None, group=False, **args):
        start = (center[0] - size[0] / 2, center[1] + size[1] / 2)
        start = self._transform(start)
        size = utils_vector.mult(size, self.ZOOM)
        rect = self.dwg.rect(insert=start, size=size, rx=rx, ry=ry, **args)
        if group:
            return group.add(rect)
        else:
            return self.dwg.add(rect)

    def circle(self, coord, r, group=False, **args):
        ###print ">>> circle() on %s with r %s with %s" % (coord, r, args)
        r = r * self.ZOOM
        circle = self.dwg.circle(center=self._transform(coord), r=r, **args)
        if group:
            return group.add(circle)
        else:
            return self.dwg.add(circle)

    def text(self, text, coord, group=False, size=None, **args):
        if size is None:
            size = self.DEFAULT_TEXT_SIZE
        ###print ">>> text() of %s to %s with %s" % (text, coord, args)
        size = size * self.ZOOM
        txt = self.dwg.text(text, insert=self._transform(coord), font_size=size, **args)
        if group:
            return group.add(txt)
        else:
            return self.dwg.add(txt)

    def arrow(self, start, end, sse_group = False, **args):
        angle = utils_vector.vectors_angle(utils_vector.rm(end, start)) - 180
        a1 = utils_vector.angle_target(end, angle-45, 1)
        a2 = utils_vector.angle_target(end, angle+45, 1)
        c1 = utils_vector.angle_target(end, angle-27, 0.8)
        c2 = utils_vector.angle_target(end, angle+27, 0.8)
        start1 = utils_vector.angle_target(start, angle-90, 0.4)
        start2 = utils_vector.angle_target(start, angle+90, 0.4)
        pom = (self._transform(c1), self._transform(start1), self._transform(start2), \
             self._transform(c2), self._transform(a2), \
             self._transform(end), self._transform(a1))
        arrow = self.dwg.polygon(pom, **args)
        if sse_group:
            return sse_group.add(arrow)
        else:
            return self.dwg.add(arrow)

    def box(self, start, end, sse_group=False, **args):
        # helix box drawing
        angle = utils_vector.vectors_angle(utils_vector.rm(end, start)) -180
        start2 = utils_vector.angle_target(start, angle+90, 0.5)
        rcenter=self._transform((start2[0],start2[1]))
        rtransform="rotate(%s %s %s)" % (int(90+angle),rcenter[0],rcenter[1])
        box = self.dwg.rect(insert=rcenter, \
            size=(int(math.sqrt((self._transform(start)[0]-self._transform(end)[0])*(self._transform(start)[0]-self._transform(end)[0])+
            (self._transform(start)[1]-self._transform(end)[1])*(self._transform(start)[1]-self._transform(end)[1]))), 10), \
            transform=rtransform,  rx=7, ry=45, **args)
        if sse_group:
            return sse_group.add(box)
        else:
            return self.dwg.add(box)

    def arrow_frac(self, start, end, order, from_number, **args):
        # per res part of arrow (input start, end, number of res from_number is the length of sses in res)
        # the forst sign is arrow the rest are bars
        if order == 0:
            angle = utils_vector.vectors_angle(utils_vector.rm(end, start)) - 180
            a1 = utils_vector.angle_target(end, angle-45, 1)
            a2 = utils_vector.angle_target(end, angle+45, 1)
            pom=(self._transform(a2), self._transform(end), self._transform(a1))
            arrow=self.dwg.polygon(pom, **args)
            return self.dwg.add(arrow)
        else:
            size = (math.sqrt((self._transform(start)[0]-self._transform(end)[0])*(self._transform(start)[0]-self._transform(end)[0])+
                   (self._transform(start)[1]-self._transform(end)[1])*(self._transform(start)[1]-self._transform(end)[1])), 20)
            angle = utils_vector.vectors_angle(utils_vector.rm(end, start)) - 180
            start2 = utils_vector.angle_target(start, angle+90, 1.5)
            start3 = utils_vector.angle_target(start2, angle, 0.05*(order+1)/from_number*size[0])
            rcenter=self._transform((start3[0], start3[1]))
            rtransform="rotate(%s %s %s)" % (90+angle, rcenter[0], rcenter[1])
            box = self.dwg.rect(insert=[rcenter[0]+size[0], rcenter[1]+size[1]], \
               size=[size[0]/from_number, 20], \
               transform=rtransform,  rx=0, ry=0, **args)
            return self.dwg.add(box)


    def secstr(self, name, coord, size, angle=None, h_type=None, \
          color=None, opacity=1.0, pdb_file="", threeddata = "", sse_group=False):
        # set stroke width, filling and colorh
        if h_type == 1:
            stroke_width = 4
        elif h_type == 2:
            stroke_width = 4
        else:
            stroke_width = 1
        fill="fill:white" if opacity == 1 else "fill-opacity:0"

        angle_start = utils_vector.angle_target(
            coord, 180+angle, size*self.CIRCLE_MULTIPLIER)
        angle_target = utils_vector.angle_target(
            coord, angle, size*self.CIRCLE_MULTIPLIER)

        # draw sses based on sse type
        if utils_sses_types.is_sheet(name):
            self.arrow(angle_start, angle_target, sse_group=sse_group, id=name, \
                style='stroke: %s; stroke-width:%s; stroke-linecap:round; %s; opacity:%s' \
                % (color, stroke_width, fill, opacity))
        if utils_sses_types.is_helic(name):
            self.box(angle_start, angle_target, sse_group=sse_group, id=name, \
                style='stroke: %s; stroke-width:%s; stroke-linecap:round; %s; opacity:%s' \
                 % (color, stroke_width, fill, opacity))


    def box_frac(self, start, end, order, from_number, **args):
        # helix box drawing
        size = (math.sqrt((self._transform(start)[0]-self._transform(end)[0])*(self._transform(start)[0]-self._transform(end)[0])+
           (self._transform(start)[1]-self._transform(end)[1])*(self._transform(start)[1]-self._transform(end)[1])), 20)
        angle = utils_vector.vectors_angle(utils_vector.rm(end, start)) - 180
        start2 = utils_vector.angle_target(start, angle+90, 1.5)
        start3 = utils_vector.angle_target(start2, angle, 0.05*(order+1)/from_number*size[0])
        rcenter=self._transform((start3[0], start3[1]))
        rtransform="rotate(%s %s %s)" % (90+angle, rcenter[0], rcenter[1])
        if order == 0:
            box = self.dwg.rect(insert=[rcenter[0]+size[0], rcenter[1]+size[1]], \
               size=[size[0]/from_number/2, 20], transform=rtransform,  rx=0, ry=0, **args)
            self.dwg.add(box)
            box = self.dwg.rect(insert=[rcenter[0]+size[0], rcenter[1]+size[1]], \
               size=[size[0]/from_number, 20], transform=rtransform,  rx=7, ry=45, **args)
        elif order + 1 == from_number:
            box = self.dwg.rect(insert=[rcenter[0]+size[0], rcenter[1]+size[1]], \
               size=[size[0]/from_number, 20], transform=rtransform,  rx=7, ry=45, **args)
            self.dwg.add(box)
            start5 = utils_vector.angle_target(start2, angle, 0.05*(order+0.5)/from_number*size[0])
            rcenter2=self._transform((start5[0], start5[1]))
            rtransform2="rotate(%s %s %s)" % (90+angle, rcenter2[0], rcenter2[1])
            box = self.dwg.rect(insert=[rcenter2[0]+size[0], rcenter2[1]+size[1]], \
               size=[size[0]/from_number/2, 20], transform=rtransform2,  rx=0, ry=0, **args)
        else:
            box = self.dwg.rect(insert=[rcenter[0]+size[0], rcenter[1]+size[1]], \
               size=[size[0]/from_number, 20], \
               transform=rtransform,  rx=0, ry=0, **args)
        return self.dwg.add(box)


    def secstrperres(self, name, coord, size, angle=None, h_type=None, \
          color=None, opacity=1.0, pdb_file="", threeddata = ""):

        # set stroke width, filling and color
        if h_type == 1:
            stroke_width = 0
        elif h_type == 2:
            stroke_width = 0
        else:
            stroke_width = 0
        fill="fill:white" if opacity == 1 else "fill-opacity:0"

        angle_start = utils_vector.angle_target(
            coord, 180+angle, size*self.CIRCLE_MULTIPLIER)
        angle_target = utils_vector.angle_target(
            coord, angle, size*self.CIRCLE_MULTIPLIER)

        # draw sses based on sse type
        if utils_sses_types.is_sheet(name):
            for c in color:
                fill=f"fill:{color[c]}" if opacity == 1 else f"fill:{color[c]}" # fill-opacity:0"
                self.arrow_frac(angle_start, angle_target, c, len(color), \
                    style='stroke: %s; stroke-width:%s; stroke-linecap:round; %s; opacity:%s' \
                    % (color[c], stroke_width, fill, opacity))
#            self.text(name, angle_start)

        if utils_sses_types.is_helic(name):
            for c in color:
                fill=f"fill:{color[c]} " if opacity == 1 else f"fill:{color[c]}" # fill-opacity:0"
                self.box_frac(angle_start, angle_target, c, len(color), \
                    style='stroke: %s; stroke-width:%s; stroke-linecap:round; %s; opacity:%s' \
                     % (color[c], stroke_width, fill, opacity))
#            self.text(name, angle_start)


    def lig_text(self, edsc, ecmp):
        if ecmp != None:
            name = ecmp
        else:
            name = edsc
        if name.startswith("DNA"):
           name = "DNA"
        elif name.startswith("RNA"):
           name = "RNA"
        elif len(name) > 6:
           name = name[:6]
        return(name)

    def ligand_diagram_size(self, size):
        r = math.sqrt(size)/10+0.1
        if size <10:
            r = 0.2+0.1
        return r

    def ligand(self, name, coord, size, e_type=None, ept=None, edsc=None, \
               ecmp=None, add_descriptions="", input_color = None, \
               ocolor = None, ligand_group=False):

        color = input_color
        cc = ocolor

        r = self.ligand_diagram_size(size)
        self.circle(coord, r, group=ligand_group, \
                style='stroke: %s; stroke-width:%s; stroke-linecap:round; fill: %s' \
                 % (cc, 2, color))

        if (add_descriptions != None) and (("ligands" not in add_descriptions) or ("protein" in add_descriptions)):
            name = self.lig_text(edsc=edsc, ecmp=ecmp)
            # if the ligand is big (e.g. HEM in 1tqn) the ligand name should
            # not collide with ligand circle - the position of text have to
            # depend on radius r as well
            self.text(str(name), (coord[0]+r, coord[1]+r+0.3), group=ligand_group, style="text-anchor: middle;")


    def channels_start(self, name, coord, previous_coord, size, text_name, opacity=1.0, \
               input_color=None, ocolor=None, channel_group=False):
        color = input_color
        cc = ocolor

        if size < 1:
            r = 0.14 *3
        else:
            r = (size)/7 *3
        r = r / 2
        self.polygon([coord[0], coord[1], previous_coord[1], previous_coord[0]], group=channel_group, \
                style='stroke: %s; stroke-width:%s; stroke-linecap:round; opacity:%s; fill: %s' \
                 % (cc, 2, opacity, color))


    def channels_middle(self, name, coord, previous_coord, p_p_coord, size, text_name, opacity=1.0, \
               input_color=None, ocolor=None, channel_group=False):
        color = input_color
        cc = ocolor

        if size < 1:
            r = 0.14 *3
        else:
            r =  (size)/7 *3
        r = r / 2
        self.polygon([coord[0], coord[1], previous_coord[1], previous_coord[0]], group = channel_group, \
                style='stroke: %s; stroke-width:%s; stroke-linecap:round; opacity:%s; fill: %s' \
                 % (cc, 2, opacity, color))


    def channels_end(self, name, coord, previous_coord, size, text_name, opacity=1.0, \
               input_color=None, ocolor=None, channel_group=False):
        color = input_color
        cc = ocolor

        if size < 1:
            r = 0.14 *3
        else:
            r =  (size)/7 *3
        r = r / 2
        self.polygon([coord[0], coord[1], previous_coord[1], previous_coord[0]], group = channel_group, \
                style='stroke: %s; stroke-width:%s; stroke-linecap:round; opacity:%s; fill: %s' \
#                 % (cc, 2, opacity, "yellow"))
                 % (cc, 2, opacity, color))


    def channels(self, name, coord, size, text_name, threed_layout,  opacity=1.0, \
               input_color = None, ocolor = None):

        color = input_color
        cc = ocolor

        if size < 1:
            r = 0.14 *3
        else:
            r =  (size)/7 *3
        self.circle(coord, r, \
                style='stroke: %s; stroke-width:%s; stroke-linecap:round; opacity:%s; fill: %s' \
                 % (cc, 2, opacity, color))
        self.circle(coord, r, \
                style='stroke: %s; stroke-width:%s; stroke-linecap:round; fill: %s' \
                 % (cc, 2, color))


    def templsecstr(self, name, coord, size,h_type=None, color=None, angle=None):
        stroke_width = 4

        # Draw sheet arrow
        if h_type == 1:
           if utils_sses_types.is_sheet(name):
               angle_start = utils_vector.angle_target(coord, 180+angle, size*self.CIRCLE_MULTIPLIER)
               angle_target = utils_vector.angle_target(coord, angle, size*self.CIRCLE_MULTIPLIER)
               self.arrow(angle_start, angle_target, \
                  style='stroke: %s; stroke-width:%s; stroke-linecap:round; fill:black; opacity:0.5' % (color, stroke_width))
#               self.circle(coord, size* self.CIRCLE_MULTIPLIER*0.25,  style='stroke: %s; stroke-width:%s; stroke-linecap:round; fill:black; opacity:0.5' % (color, stroke_width))
#               self.text(name, coord)
           if utils_sses_types.is_helic(name):
               angle_start = utils_vector.angle_target(coord,180+angle, size*self.CIRCLE_MULTIPLIER)
               angle_target = utils_vector.angle_target(coord, angle, size*self.CIRCLE_MULTIPLIER)
               self.box(angle_start, angle_target, \
                  style='stroke: %s; stroke-width:%s; stroke-linecap:round; fill:black; opacity:0.5' % (color, stroke_width))
#               self.circle(coord, size* self.CIRCLE_MULTIPLIER*0.25,  style='stroke: %s; stroke-width:%s; stroke-linecap:round; fill:black; opacity:0.5' % (color, stroke_width))
#               self.text(name, coord)


    def helic_connector(self, coord_from, coord_to, size_from, size_to, angle_from, angle_to, color='black', opacity = 1.0, sse_group=False):
        angle_start = utils_vector.angle_target(coord_from,
                                                     180+angle_from,
                                                     size_from*self.CIRCLE_MULTIPLIER)
        angle_target = utils_vector.angle_target(coord_to,
                                                     angle_to,
                                                     size_to*self.CIRCLE_MULTIPLIER)
        self.line(angle_start, angle_target, group=sse_group, style="stroke: %s; opacity: %s" % (color, opacity))


    def channel_connector(self, coord_from, coord_to, color, size1, size2, opacity = 1.0):
        size = min(size1,size2) * 7 *3
        self.line(coord_from, coord_to, style="stroke: %s; stroke-width: %s; " % (color, size))
#        self.line(coord_from, coord_to, style="stroke: %s; opacity: %s; stroke-width: %s; " % (color, opacity, size))


    def ligand_connector(self, coord_from, coord_to, color='black'):
        self.line(coord_from, coord_to, style="stroke: %s" % (color))


if __name__ == "__main__":
    img = Image('aaa.svg', (50, 50))
    img.line((0, 0), (10, 10), style="stroke: red;")
    img.text("hello", (10, 0), 1)
    img.text("hello", (10, -10), 5)
    img.text("hello", (10, -20), 2.5, style="fill: red; font-family: Verdana; font-weight: bold; alignment-baseline: central; text-anchor: middle;")
    img.circle((-10, -10), 5, style="stroke: green;")
    img.arrow((0, 5), (-10, 10), style="stroke: blue;")
    img.arrow((0, 5), (-8, 10), marker='default', style="stroke: green;")
    img.rect((-19, 19), (10, 10), style="fill: white; fill-opacity: 0.9; stroke: black; stroke-width: 3;")
    img.rect((-15, 15), (7.5, 7.5), rx=5, ry=5, style="fill: green; fill-opacity: 0.3; stroke: gray; stroke-width: 1;")
    img.save()
