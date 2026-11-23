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
        coord = [round(coord[0], 2), round(coord[1], 2)]
        return coord

    @staticmethod
    def _parse_rgb(color_str):
        """Parse hex or common CSS color name into (r, g, b) tuple in 0..255."""
        if not color_str or not isinstance(color_str, str):
            return (160, 160, 160)
        c = color_str.strip().lower()
        named = {
            "red": (220, 50, 47),
            "blue": (38, 139, 210),
            "green": (42, 161, 152),
            "grey": (135, 142, 150),
            "gray": (135, 142, 150),
            "black": (40, 44, 52),
            "white": (255, 255, 255),
            "yellow": (215, 165, 30),
            "orange": (203, 75, 22),
            "purple": (108, 113, 196),
            "cyan": (42, 161, 152),
            "magenta": (211, 54, 130),
            "brown": (140, 95, 60),
        }
        if c in named:
            return named[c]
        if c.startswith("#") and len(c) == 7:
            try:
                return (int(c[1:3], 16), int(c[3:5], 16), int(c[5:7], 16))
            except ValueError:
                pass
        return (135, 142, 150)

    @classmethod
    def _tint_color(cls, color_str, blend=0.28):
        """Blend color with white to create a soft, readable fill for 2D SSE elements."""
        r, g, b = cls._parse_rgb(color_str)
        tr = int(round(r * blend + 255 * (1.0 - blend)))
        tg = int(round(g * blend + 255 * (1.0 - blend)))
        tb = int(round(b * blend + 255 * (1.0 - blend)))
        return "#%02x%02x%02x" % (tr, tg, tb)

    @classmethod
    def _darken_color(cls, color_str, factor=0.82):
        """Slightly darken stroke color for crisp SSE borders."""
        r, g, b = cls._parse_rgb(color_str)
        return "#%02x%02x%02x" % (
            max(0, min(255, int(round(r * factor)))),
            max(0, min(255, int(round(g * factor)))),
            max(0, min(255, int(round(b * factor)))),
        )

    def membrane_band(self, y_center, half_height=3.0, label="MEMBRÁNA"):
        """Draw a horizontal semi-transparent membrane bilayer band across the canvas."""
        ll = self.get_lower_left()
        ur = self.get_upper_right()
        margin = 0.6
        x_left = ll[0] + margin
        x_right = ur[0] - margin
        width = max(1.0, x_right - x_left)
        height = max(1.0, 2.0 * half_height)

        mem_group = self.dwg.add(self.dwg.g(id="membrane_layer", class_="membrane"))

        # Main semi-transparent hydrophobic core band
        self.rect(
            center=((x_left + x_right) / 2.0, y_center),
            size=(width, height),
            rx=4,
            ry=4,
            group=mem_group,
            fill_opacity=0.14,
            style="fill:#fdf8e6;fill-opacity:0.14;stroke:none;",
        )

        # Upper and lower lipid headgroup boundary lines
        y_top = y_center + half_height
        y_bot = y_center - half_height
        self.line(
            (x_left, y_top),
            (x_right, y_top),
            group=mem_group,
            stroke_opacity=0.35,
            stroke_dasharray="7,4",
            style="stroke:#dcc47e;stroke-width:1.2;stroke-dasharray:7,4;opacity:0.35;",
        )
        self.line(
            (x_left, y_bot),
            (x_right, y_bot),
            group=mem_group,
            stroke_opacity=0.35,
            stroke_dasharray="7,4",
            style="stroke:#dcc47e;stroke-width:1.2;stroke-dasharray:7,4;opacity:0.35;",
        )
        # Subtle midplane line
        self.line(
            (x_left, y_center),
            (x_right, y_center),
            group=mem_group,
            stroke_opacity=0.22,
            stroke_dasharray="3,5",
            style="stroke:#e8d8a0;stroke-width:0.8;stroke-dasharray:3,5;opacity:0.22;",
        )

        if label:
            self.text(
                label,
                (x_left + 0.5, y_top - 0.65),
                group=mem_group,
                size=0.58,
                fill_opacity=0.45,
                style="fill:#b59b52;stroke:none;font-family:Verdana, sans-serif;font-weight:bold;opacity:0.45;",
            )

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
        angle = utils_vector.vectors_angle(utils_vector.rm(end, start)) - 180
        start2 = utils_vector.angle_target(start, angle+90, 0.5)
        rcenter = self._transform((start2[0], start2[1]))
        rtransform = "rotate(%s %s %s)" % (round(90 + angle, 2), rcenter[0], rcenter[1])
        t_start = self._transform(start)
        t_end = self._transform(end)
        length_px = round(math.hypot(t_start[0] - t_end[0], t_start[1] - t_end[1]), 2)
        box = self.dwg.rect(insert=rcenter, \
            size=(length_px, 11), \
            transform=rtransform, rx=7, ry=45, **args)
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
          color=None, opacity=1.0, pdb_file="", threeddata = "", sse_group=False, is_tm=False):
        # Set stroke width, tinted fill and halo highlighting based on SSE importance & TM status
        if is_tm:
            stroke_width = 4.6
            blend = 0.34
            stroke_color = self._darken_color(color, 0.85) if color == "grey" else color
        elif h_type == 1:
            stroke_width = 4.0
            blend = 0.26
            stroke_color = color
        elif h_type == 2:
            stroke_width = 2.6
            blend = 0.14
            stroke_color = "#7a828a" if color == "grey" else color
        else:
            stroke_width = 1.2
            blend = 0.08
            stroke_color = color

        if opacity == 1:
            tinted_fill = self._tint_color(stroke_color, blend=blend)
            fill = f"fill:{tinted_fill}"
        else:
            fill = "fill-opacity:0"

        angle_start = utils_vector.angle_target(
            coord, 180+angle, size*self.CIRCLE_MULTIPLIER)
        angle_target = utils_vector.angle_target(
            coord, angle, size*self.CIRCLE_MULTIPLIER)

        # Draw white halo underlay on single diagrams so connector lines don't visually clash with SSEs
        if opacity == 1:
            halo_width = stroke_width + 3.5
            halo_style = f"stroke: #ffffff; stroke-width:{halo_width}; stroke-linecap:round; stroke-linejoin:round; fill:#ffffff; opacity:0.92"
            if utils_sses_types.is_sheet(name):
                self.arrow(angle_start, angle_target, sse_group=sse_group, style=halo_style)
            if utils_sses_types.is_helic(name):
                self.box(angle_start, angle_target, sse_group=sse_group, style=halo_style)

        # Draw main highlighted SSE shape
        fill_op = 1.0 if opacity == 1 else 0.0
        if utils_sses_types.is_sheet(name):
            self.arrow(angle_start, angle_target, sse_group=sse_group, id=name,
                stroke_opacity=opacity, fill_opacity=fill_op,
                style='stroke:%s;stroke-width:%s;stroke-linecap:round;stroke-linejoin:round;%s;opacity:%s'
                % (stroke_color, stroke_width, fill, opacity))
        if utils_sses_types.is_helic(name):
            self.box(angle_start, angle_target, sse_group=sse_group, id=name,
                stroke_opacity=opacity, fill_opacity=fill_op,
                style='stroke:%s;stroke-width:%s;stroke-linecap:round;stroke-linejoin:round;%s;opacity:%s'
                 % (stroke_color, stroke_width, fill, opacity))

        # Add clean short label for primary or transmembrane SSEs in single-diagram view
        if opacity == 1 and (h_type == 1 or is_tm) and size >= 7:
            short_label = name.split("_")[0]
            label_color = self._darken_color(stroke_color, 0.38)
            self.text(
                short_label,
                (coord[0], coord[1] - 0.14),
                group=sse_group,
                size=0.44,
                style=f"fill: {label_color}; stroke: none; font-family: Verdana, sans-serif; font-weight: bold; text-anchor: middle;",
            )


    def ghost_secstr(self, name, coord, size, angle=0.0, color="#667788",
                     partner_coord=None, partner_size=None, partner_angle=None,
                     origin_coord=None, show_arc=False):
        """
        Draws a translucent, dashed 'ghost' secondary structure (`duch interagující SSE`)
        to indicate local 3D interactions (such as beta-barrel closure, inward-folded pore
        loops, or split beta-sheets) that could not fit directly together in 2D.
        """
        ghost_group = self.dwg.g(class_="ghost_sse")
        stroke_color = "#6b7682" if (not color or color == "grey") else self._darken_color(color, 0.85)
        tinted_fill = self._tint_color(stroke_color, blend=0.22)

        angle_start = utils_vector.angle_target(
            coord, 180.0 + angle, size * self.CIRCLE_MULTIPLIER
        )
        angle_target = utils_vector.angle_target(
            coord, angle, size * self.CIRCLE_MULTIPLIER
        )

        # 1. Optional curved dashed folding arc from the real SSE to its docked ghost position
        if show_arc and origin_coord is not None:
            p0 = self._transform(origin_coord)
            p2 = self._transform(coord)
            mx = 0.5 * (p0[0] + p2[0])
            my = min(p0[1], p2[1]) - 22.0
            d_arc = f"M {p0[0]:.1f},{p0[1]:.1f} Q {mx:.1f},{my:.1f} {p2[0]:.1f},{p2[1]:.1f}"
            ghost_group.add(
                self.dwg.path(
                    d=d_arc,
                    stroke_dasharray="7,5",
                    stroke_opacity=0.62,
                    fill="none",
                    style=f"fill:none;stroke:{stroke_color};stroke-width:1.8;stroke-dasharray:7,5;opacity:0.62;",
                )
            )

        # 2. Short dashed H-bond / contact ladder rungs between ghost SSE and partner SSE
        if partner_coord is not None and partner_size is not None and partner_angle is not None:
            p_start = utils_vector.angle_target(
                partner_coord, 180.0 + partner_angle, partner_size * self.CIRCLE_MULTIPLIER
            )
            p_target = utils_vector.angle_target(
                partner_coord, partner_angle, partner_size * self.CIRCLE_MULTIPLIER
            )
            # Ensure g_a->g_b and p_a->p_b traverse in the same geometric direction for ladder rungs
            gv = (angle_target[0] - angle_start[0], angle_target[1] - angle_start[1])
            pv = (p_target[0] - p_start[0], p_target[1] - p_start[1])
            if gv[0] * pv[0] + gv[1] * pv[1] >= 0.0:
                p_a, p_b = p_start, p_target
            else:
                p_a, p_b = p_target, p_start

            pdx, pdy = p_b[0] - p_a[0], p_b[1] - p_a[1]
            plen2 = pdx * pdx + pdy * pdy
            fracs = (0.25, 0.50, 0.75)
            t_projs = []
            for frac in fracs:
                gx = angle_start[0] + frac * (angle_target[0] - angle_start[0])
                gy = angle_start[1] + frac * (angle_target[1] - angle_start[1])
                if plen2 > 1e-6:
                    tp = ((gx - p_a[0]) * pdx + (gy - p_a[1]) * pdy) / plen2
                    tp = max(0.12, min(0.88, tp))
                else:
                    tp = frac
                t_projs.append(tp)

            # Avoid V-fan pinching when ghost SSE sits near one end of a partner strand
            if abs(t_projs[-1] - t_projs[0]) < 0.22:
                mean_tp = sum(t_projs) / 3.0
                if mean_tp >= 0.5:
                    t_projs = [0.64, 0.76, 0.88]
                else:
                    t_projs = [0.12, 0.24, 0.36]

            for frac, t_proj in zip(fracs, t_projs):
                gx = angle_start[0] + frac * (angle_target[0] - angle_start[0])
                gy = angle_start[1] + frac * (angle_target[1] - angle_start[1])
                px = p_a[0] + t_proj * pdx
                py = p_a[1] + t_proj * pdy
                # Shorten rung ends slightly so they touch the outer edges of the two SSEs
                rx, ry = px - gx, py - gy
                rlen = math.hypot(rx, ry)
                if 0.80 <= rlen <= 5.5:
                    pad = min(0.42, 0.22 * rlen)
                    ux, uy = rx / rlen, ry / rlen
                    r0 = (gx + ux * pad, gy + uy * pad)
                    r1 = (px - ux * pad, py - uy * pad)
                    ghost_group.add(
                        self._line(
                            r0,
                            r1,
                            stroke_dasharray="4,3.5",
                            stroke_opacity=0.72,
                            style="stroke:#52606e;stroke-width:1.7;stroke-dasharray:4,3.5;opacity:0.72;",
                        )
                    )

        # 3. Soft white backdrop halo so background connector lines don't cut through the ghost
        halo_style = "stroke:#ffffff;stroke-width:5.0;stroke-linecap:round;stroke-linejoin:round;fill:#ffffff;opacity:0.78"
        if utils_sses_types.is_sheet(name):
            self.arrow(angle_start, angle_target, sse_group=ghost_group, fill_opacity=0.78, stroke_opacity=0.78, style=halo_style)
        if utils_sses_types.is_helic(name):
            self.box(angle_start, angle_target, sse_group=ghost_group, fill_opacity=0.78, stroke_opacity=0.78, style=halo_style)

        # 4. Dashed translucent ghost SSE body (pass stroke_dasharray and opacities as SVG attrs too)
        ghost_style = (
            f"stroke:{stroke_color};stroke-width:2.3;stroke-dasharray:6.5,4.5;"
            f"stroke-linecap:butt;stroke-linejoin:round;fill:{tinted_fill};"
            f"fill-opacity:0.28;stroke-opacity:0.80"
        )
        if utils_sses_types.is_sheet(name):
            self.arrow(
                angle_start,
                angle_target,
                sse_group=ghost_group,
                stroke_dasharray="6.5,4.5",
                fill_opacity=0.28,
                stroke_opacity=0.80,
                style=ghost_style,
            )
        if utils_sses_types.is_helic(name):
            self.box(
                angle_start,
                angle_target,
                sse_group=ghost_group,
                stroke_dasharray="6.5,4.5",
                fill_opacity=0.28,
                stroke_opacity=0.80,
                style=ghost_style,
            )

        # 5. Bracketed italic ghost label e.g. [E34]
        short_label = f"[{name.split('_')[0]}]"
        label_color = self._darken_color(stroke_color, 0.42)
        self.text(
            short_label,
            (coord[0], coord[1] - 0.14),
            group=ghost_group,
            size=0.40,
            style=(
                "fill:#ffffff;stroke:#ffffff;stroke-width:2.5;stroke-linejoin:round;"
                "font-family:Verdana, sans-serif;font-style:italic;font-weight:bold;"
                "text-anchor:middle;opacity:0.85;"
            ),
        )
        self.text(
            short_label,
            (coord[0], coord[1] - 0.14),
            group=ghost_group,
            size=0.40,
            style=(
                f"fill:{label_color};stroke:none;font-family:Verdana, sans-serif;"
                "font-style:italic;font-weight:bold;text-anchor:middle;opacity:0.82;"
            ),
        )
        self.dwg.add(ghost_group)



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
        # Subtle white halo around ligand circle so it stands out cleanly against connectors/membrane
        self.circle(coord, r + 0.10, group=ligand_group,
                    fill_opacity=0.85, stroke_opacity=0.85,
                    style="stroke:#ffffff;stroke-width:2.5;fill:#ffffff;opacity:0.85;")
        self.circle(coord, r, group=ligand_group, \
                style='stroke: %s; stroke-width:%s; stroke-linecap:round; fill: %s' \
                 % (cc, 2, color))

        if (add_descriptions != None) and (("ligands" not in add_descriptions) or ("protein" in add_descriptions)):
            name = self.lig_text(edsc=edsc, ecmp=ecmp)
            # Center ligand label directly above the ligand circle with a clean readability halo
            lbl_coord = (coord[0], coord[1] + r + 0.30)
            self.text(str(name), lbl_coord, group=ligand_group, size=0.82,
                      style="text-anchor: middle; stroke: #ffffff; stroke-width: 3px; stroke-linejoin: round; fill: #ffffff; font-family: Verdana, sans-serif; font-weight: bold;")
            self.text(str(name), lbl_coord, group=ligand_group, size=0.82,
                      style="text-anchor: middle; stroke: none; font-family: Verdana, sans-serif; font-weight: bold;")


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


    def tunnel_ribbon(self, points, radii, color="#00838f", opacity=1.0, channel_group=False):
        """
        Render a protein tunnel/channel as a smooth variable-width SVG ribbon with:
          - Catmull-Rom centerline and radius interpolation,
          - rounded semicircular endcaps,
          - soft white outer halo separating the tunnel from background connectors,
          - translucent tinted body fill and crisp contour,
          - subtle dashed inner centerline indicating the tunnel flow axis.
        """
        if not points or len(points) < 2:
            return

        # Filter near-duplicate consecutive points
        clean_pts = [points[0]]
        clean_rads = [float(radii[0]) if radii else 1.6]
        for i in range(1, len(points)):
            p = points[i]
            r = float(radii[i]) if radii and i < len(radii) else clean_rads[-1]
            if math.hypot(p[0] - clean_pts[-1][0], p[1] - clean_pts[-1][1]) > 0.12 or i == len(points) - 1:
                clean_pts.append(p)
                clean_rads.append(r)

        if len(clean_pts) > 2 and math.hypot(clean_pts[-1][0] - clean_pts[-2][0], clean_pts[-1][1] - clean_pts[-2][1]) < 0.08:
            clean_pts.pop(-2)
            clean_rads.pop(-2)

        if len(clean_pts) < 2:
            return

        # Subdivide along a smooth Catmull-Rom spline
        sub_pts = []
        sub_rads = []
        n = len(clean_pts)
        steps_per_seg = 5
        for i in range(n - 1):
            p0 = clean_pts[max(0, i - 1)]
            p1 = clean_pts[i]
            p2 = clean_pts[i + 1]
            p3 = clean_pts[min(n - 1, i + 2)]
            r1 = clean_rads[i]
            r2 = clean_rads[i + 1]
            for s in range(steps_per_seg):
                t = float(s) / steps_per_seg
                t2 = t * t
                t3 = t2 * t
                # Catmull-Rom basis (tension 0.5)
                x = 0.5 * (
                    (2.0 * p1[0])
                    + (-p0[0] + p2[0]) * t
                    + (2.0 * p0[0] - 5.0 * p1[0] + 4.0 * p2[0] - p3[0]) * t2
                    + (-p0[0] + 3.0 * p1[0] - 3.0 * p2[0] + p3[0]) * t3
                )
                y = 0.5 * (
                    (2.0 * p1[1])
                    + (-p0[1] + p2[1]) * t
                    + (2.0 * p0[1] - 5.0 * p1[1] + 4.0 * p2[1] - p3[1]) * t2
                    + (-p0[1] + 3.0 * p1[1] - 3.0 * p2[1] + p3[1]) * t3
                )
                sub_pts.append((x, y))
                sub_rads.append(r1 + (r2 - r1) * t)
        sub_pts.append(clean_pts[-1])
        sub_rads.append(clean_rads[-1])

        # Smooth radius profile
        m = len(sub_pts)
        for _ in range(2):
            new_r = list(sub_rads)
            for k in range(1, m - 1):
                new_r[k] = 0.25 * sub_rads[k - 1] + 0.50 * sub_rads[k] + 0.25 * sub_rads[k + 1]
            sub_rads = new_r

        # Transform centerline points and half-widths to SVG pixel space
        svg_pts = [self._transform(pt) for pt in sub_pts]
        svg_ws = [max(0.24, min(0.75, 0.19 * r)) * self.ZOOM for r in sub_rads]

        # Compute smoothed unit tangents and normals in SVG pixel space
        tangents = []
        for k in range(m):
            p_prev = svg_pts[max(0, k - 2)]
            p_next = svg_pts[min(m - 1, k + 2)]
            tx = p_next[0] - p_prev[0]
            ty = p_next[1] - p_prev[1]
            tnorm = math.hypot(tx, ty)
            if tnorm < 1e-4:
                tangents.append(tangents[-1] if tangents else (1.0, 0.0))
            else:
                cand_t = (tx / tnorm, ty / tnorm)
                if tangents and (cand_t[0] * tangents[-1][0] + cand_t[1] * tangents[-1][1]) < -0.2:
                    cand_t = tangents[-1]
                tangents.append(cand_t)

        # Smooth tangents so normals never kink
        for _ in range(3):
            sm_t = list(tangents)
            for k in range(1, m - 1):
                sx = 0.25 * tangents[k - 1][0] + 0.50 * tangents[k][0] + 0.25 * tangents[k + 1][0]
                sy = 0.25 * tangents[k - 1][1] + 0.50 * tangents[k][1] + 0.25 * tangents[k + 1][1]
                sn = math.hypot(sx, sy)
                if sn > 1e-4:
                    sm_t[k] = (sx / sn, sy / sn)
            tangents = sm_t

        normals = [(-ty, tx) for (tx, ty) in tangents]

        # Build left boundary, end cap, right boundary, and start cap
        outline = []
        for k in range(m):
            px, py = svg_pts[k]
            nx, ny = normals[k]
            w = svg_ws[k]
            outline.append((round(px + nx * w, 2), round(py + ny * w, 2)))

        # Semicircular cap at end (m - 1)
        px_e, py_e = svg_pts[-1]
        tx_e, ty_e = tangents[-1]
        nx_e, ny_e = normals[-1]
        w_e = svg_ws[-1]
        cap_steps = 6
        for s in range(1, cap_steps):
            ang = (math.pi * s) / cap_steps
            ca, sa = math.cos(ang), math.sin(ang)
            cx = px_e + w_e * (ca * nx_e + sa * tx_e)
            cy = py_e + w_e * (ca * ny_e + sa * ty_e)
            outline.append((round(cx, 2), round(cy, 2)))

        for k in range(m - 1, -1, -1):
            px, py = svg_pts[k]
            nx, ny = normals[k]
            w = svg_ws[k]
            outline.append((round(px - nx * w, 2), round(py - ny * w, 2)))

        # Semicircular cap at start (0)
        px_0, py_0 = svg_pts[0]
        tx_0, ty_0 = tangents[0]
        nx_0, ny_0 = normals[0]
        w_0 = svg_ws[0]
        for s in range(1, cap_steps):
            ang = (math.pi * s) / cap_steps
            ca, sa = math.cos(ang), math.sin(ang)
            cx = px_0 + w_0 * (-ca * nx_0 - sa * tx_0)
            cy = py_0 + w_0 * (-ca * ny_0 - sa * ty_0)
            outline.append((round(cx, 2), round(cy, 2)))

        d_env = "M %s,%s " % outline[0] + " ".join("L %s,%s" % pt for pt in outline[1:]) + " Z"
        d_center = "M %s,%s " % (svg_pts[0][0], svg_pts[0][1]) + " ".join(
            "L %s,%s" % (pt[0], pt[1]) for pt in svg_pts[1:]
        )

        base_color = color if (color and str(color).lower() not in ("none", "grey", "gray", "#665646")) else "#00838f"
        fill_col = self._tint_color(base_color, blend=0.45)
        stroke_col = self._darken_color(base_color, factor=0.88)
        target_container = channel_group if channel_group else self.dwg

        # 1. Soft white halo around tunnel envelope
        halo_path = self.dwg.path(
            d=d_env,
            style="fill:none;stroke:#ffffff;stroke-width:2.8;stroke-opacity:%.2f;stroke-linejoin:round;stroke-linecap:round;"
            % (0.65 * opacity),
        )
        target_container.add(halo_path)

        # 2. Translucent variable-width tunnel envelope + crisp contour
        body_path = self.dwg.path(
            d=d_env,
            style="fill:%s;fill-opacity:%.2f;stroke:%s;stroke-width:1.45;stroke-opacity:%.2f;stroke-linejoin:round;stroke-linecap:round;"
            % (fill_col, 0.34 * opacity, stroke_col, 0.82 * opacity),
        )
        target_container.add(body_path)

        # 3. Inner dashed centerline showing tunnel trajectory
        center_path = self.dwg.path(
            d=d_center,
            style="fill:none;stroke:%s;stroke-width:1.1;stroke-dasharray:3.5,3.0;stroke-opacity:%.2f;stroke-linecap:round;"
            % (stroke_col, 0.58 * opacity),
        )
        target_container.add(center_path)


    def tunnel_label(self, text_name, coord, color="#00838f", group=False):
        """Render a crisp, high-contrast label badge at the mouth of a tunnel using two-pass halo text."""
        if not text_name:
            return
        base_color = color if (color and str(color).lower() not in ("none", "grey", "gray", "#665646")) else "#00838f"
        label_col = self._darken_color(base_color, factor=0.72)
        # 1. White halo pass underneath (compatible with all SVG/PDF renderers including PyMuPDF)
        self.text(
            str(text_name),
            coord,
            size=0.66,
            group=group,
            style="fill:#ffffff;stroke:#ffffff;stroke-width:2.8px;stroke-linejoin:round;stroke-linecap:round;font-family:Verdana, sans-serif;font-weight:bold;text-anchor:middle;dominant-baseline:central;",
        )
        # 2. Crisp colored text pass on top
        self.text(
            str(text_name),
            coord,
            size=0.66,
            group=group,
            style="fill:%s;stroke:none;font-family:Verdana, sans-serif;font-weight:bold;text-anchor:middle;dominant-baseline:central;"
            % label_col,
        )


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
        conn_opacity = round(opacity * 0.65, 3) if opacity == 1.0 else opacity
        self.line(angle_start, angle_target, group=sse_group, style="stroke: %s; stroke-width: 1.3; stroke-linecap: round; opacity: %s" % (color, conn_opacity))


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
