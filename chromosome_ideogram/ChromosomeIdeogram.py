from typing import Union, Dict, List, Optional
from collections import defaultdict
from os import makedirs, getcwd
from os.path import basename
import numpy as np
import matplotlib
import matplotlib.ticker
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.markers import MarkerStyle
from matplotlib.collections import PatchCollection
from matplotlib.transforms import Affine2D
import builtins
import pyBigWig
from click import echo


def __data_preprocessing(
    reversed_chr: dict,
    bed_file: str,
    sample_name: Optional[str],
    out_path: str
) -> None:
    with open(bed_file) as f, open(f'{out_path}/{basename(bed_file)}.processed', 'w') as o:
        for line in f:
            if not line.startswith('#') and line.strip():
                split = line.strip().split('\t')
                chr_name, start, end = split[0], int(split[1]), int(split[2])
                try:
                    split[1] = str(reversed_chr[chr_name] - end)
                    split[2] = str(reversed_chr[chr_name] - start)
                except KeyError:
                    pass
                if sample_name is not None:
                    split[0] = f'{sample_name}_{split[0]}'
                echo('\t'.join(split), o)


def data_preprocessing(
    karyotype: str,
    annotation: Optional[str] = None,
    synteny: Optional[str] = None,
    mark_file: Optional[str] = None,
    sample_name: Optional[str] = None,
    exclude_chr: Optional[Union[str, List[str]]] = None,
    force_revers_chr: Optional[Union[str, List[str]]] = None,
    out_path: Optional[str] = None,
) -> None:
    if exclude_chr is None:
        exclude_chr = []
    if force_revers_chr is None:
        force_revers_chr = []
    if out_path is None:
        out_path = getcwd()
    makedirs(out_path, exist_ok=True)

    reversed_chr = {}
    with open(karyotype) as f, open(f'{out_path}/{basename(karyotype)}.processed', 'w') as o:
        for line in f:
            if not line.startswith('#') and line.strip():
                split = line.strip().split('\t')
                chr_name, chr_len, cent_start, cent_end = split[0], int(split[1]), int(split[2]), int(split[3])
                if chr_len - cent_end > cent_start and chr_name not in exclude_chr:  # short arm towards the right
                    split[2] = str(chr_len - cent_end)
                    split[3] = str(chr_len - cent_start)
                    reversed_chr[split[0]] = chr_len
                elif chr_name in force_revers_chr:
                    split[2] = str(chr_len - cent_end)
                    split[3] = str(chr_len - cent_start)
                    reversed_chr[split[0]] = chr_len
                if sample_name is not None:
                    split[0] = f'{sample_name}_{chr_name}'
                echo('\t'.join(split), o)
    echo(f'Reversed chromosome: {" ".join(reversed_chr.keys())}', err=True)

    if annotation is not None:
        __data_preprocessing(
            reversed_chr=reversed_chr,
            bed_file=annotation,
            sample_name=sample_name,
            out_path=out_path
        )

    if mark_file is not None:
        __data_preprocessing(
            reversed_chr=reversed_chr,
            bed_file=mark_file,
            sample_name=sample_name,
            out_path=out_path
        )

    if synteny is not None:
        with open(synteny) as f, open(f'{out_path}/{basename(synteny)}.processed', 'w') as o:
            for line in f:
                if not line.startswith('#') and line.strip():
                    split = line.strip().split('\t')
                    chr1, start1, end1 = split[0], int(split[1]), int(split[2])
                    chr2, start2, end2 = split[3], int(split[4]), int(split[5])
                try:
                    split[1] = str(reversed_chr[chr1] - end1)
                    split[2] = str(reversed_chr[chr1] - start1)
                except KeyError:
                    pass
                try:
                    split[4] = str(reversed_chr[chr2] - end2)
                    split[5] = str(reversed_chr[chr2] - start2)
                except KeyError:
                    pass
                if sample_name is not None:
                    split[0] = f'{sample_name}_{split[0]}'
                    split[3] = f'{sample_name}_{split[3]}'
                echo('\t'.join(split), o)


class ChromosomeIdeogram:
    def __init__(
        self,
        karyotype: Union[str, Dict[str, List[int]]],
        radius: Optional[Union[int, float]] = 0.4,
        space: int = 4
    ):
        if isinstance(karyotype, dict):
            raw_karyotype = karyotype  # {chromosome: [length, centromere_start, centromere_end]}
        else:
            raw_karyotype = {}
            with open(karyotype) as f:
                for line in f:
                    split = line.strip().split('\t')
                    if len(split) < 4:
                        continue
                    chr_name, chr_len, cent_start, cent_end = split[0], split[1], split[2], split[3]
                    chr_len, cent_start, cent_end = int(chr_len), int(cent_start), int(cent_end)
                    raw_karyotype[chr_name] = [chr_len, cent_start, cent_end]
        if not raw_karyotype:
            raise ValueError('karyotype is empty')
        longest_chr = max([i[0] for i in raw_karyotype.values()])
        if longest_chr / 10 ** 3 < 1000:
            scale = 10 ** 3
            length_unit = 'Kb'
        elif longest_chr / 10 ** 6 < 1000:
            scale = 10 ** 6
            length_unit = 'Mb'
        else:
            scale = 10 ** 9
            length_unit = 'Gb'
        scaled_karyotype = {}
        for k, v in raw_karyotype.items():
            scaled_karyotype[k] = [i / scale for i in v]
        coordinate = [radius * (space * n - 2) for n in range(1, len(scaled_karyotype) + 1)]
        chr_coordinate = {j: i + radius for i, j in zip(coordinate, scaled_karyotype.keys())}

        self.karyotype = scaled_karyotype  # {chr: [length_scaled, cent_start_scaled, cent_end_scaled]}
        self.raw_karyotype = raw_karyotype  # {chr: [length, cent_start, cent_end]}
        self.scale = scale
        self.length_unit = length_unit
        self.chr_coordinate = chr_coordinate  # horizontal layout: y centers
        self.chr_ycoordinate = None  # vertical layout: x centers
        self.chromosome_patches = {}
        self.radius = radius
        self.handles = []
        self.orientation = 'horizontal'  # 'horizontal' or 'vertical'
        self.space = space

    def __chromosome_shape(
        self,
        chromosome_len: Union[int, float],
        centromere_start: Union[int, float],
        centromere_end: Union[int, float],
        y_center: Optional[Union[int, float]] = 0,
        centromere_angle: Optional[Union[int, float]] = 60.0
    ) -> tuple:
        angle_rad = np.radians(centromere_angle)
        slope = np.tan(angle_rad)
        r = self.radius

        left_line_x = r
        right_line_x = chromosome_len - r

        mid_x = (centromere_start + centromere_end) / 2.0
        dx = mid_x - centromere_start
        h_top = max(r - dx * slope, 0.0)
        h_bottom = -h_top

        theta_right = np.linspace(np.pi / 2, -np.pi / 2, 30)
        right_arc_x = right_line_x + r * np.cos(theta_right)
        right_arc_y = y_center + r * np.sin(theta_right)

        theta_left = np.linspace(3 * np.pi / 2, np.pi / 2, 30)
        left_arc_x = left_line_x + r * np.cos(theta_left)
        left_arc_y = y_center + r * np.sin(theta_left)

        vertices = []
        vertices.append((left_line_x, y_center + r))
        vertices.append((centromere_start, y_center + r))
        vertices.append((mid_x, y_center + h_top))
        vertices.append((centromere_end, y_center + r))
        vertices.append((right_line_x, y_center + r))

        for x, y in zip(right_arc_x, right_arc_y):
            vertices.append((x, y))

        vertices.append((right_line_x, y_center - r))
        vertices.append((centromere_end, y_center - r))
        vertices.append((mid_x, y_center + h_bottom))
        vertices.append((centromere_start, y_center - r))
        vertices.append((left_line_x, y_center - r))

        for x, y in zip(left_arc_x, left_arc_y):
            vertices.append((x, y))

        vertices.append(vertices[0])

        codes = [Path.MOVETO] + [Path.LINETO] * (len(vertices) - 2) + [Path.CLOSEPOLY]
        return np.array(vertices), codes

    def __chromosome_shape_v(
        self,
        chromosome_len: Union[int, float],
        centromere_start: Union[int, float],
        centromere_end: Union[int, float],
        x_center: Optional[Union[int, float]] = 0,
        centromere_angle: Optional[Union[int, float]] = 60.0
    ) -> tuple:
        angle_rad = np.radians(centromere_angle)
        slope = np.tan(angle_rad)
        r = self.radius

        bottom_line_y = r
        top_line_y = chromosome_len - r

        mid_y = (centromere_start + centromere_end) / 2.0
        dy = mid_y - centromere_start
        h_right = max(r - dy * slope, 0.0)
        h_left = -h_right

        theta_top = np.linspace(0, np.pi, 30)
        top_arc_x = x_center + r * np.cos(theta_top)
        top_arc_y = top_line_y + r * np.sin(theta_top)

        theta_bottom = np.linspace(np.pi, 2 * np.pi, 30)
        bottom_arc_x = x_center + r * np.cos(theta_bottom)
        bottom_arc_y = bottom_line_y + r * np.sin(theta_bottom)

        vertices = []
        vertices.append((x_center + r, bottom_line_y))
        vertices.append((x_center + r, centromere_start))
        vertices.append((x_center + h_right, mid_y))
        vertices.append((x_center + r, centromere_end))
        vertices.append((x_center + r, top_line_y))

        for x, y in zip(top_arc_x, top_arc_y):
            vertices.append((x, y))

        vertices.append((x_center - r, top_line_y))
        vertices.append((x_center - r, centromere_end))
        vertices.append((x_center + h_left, mid_y))
        vertices.append((x_center - r, centromere_start))
        vertices.append((x_center - r, bottom_line_y))

        for x, y in zip(bottom_arc_x, bottom_arc_y):
            vertices.append((x, y))

        vertices.append(vertices[0])

        codes = [Path.MOVETO] + [Path.LINETO] * (len(vertices) - 2) + [Path.CLOSEPOLY]
        return np.array(vertices), codes

    def draw_chromh(
        self,
        axes,
        centromere_angle: Optional[Union[int, float]] = 60.0,
        facecolor: Optional[str] = 'lightgrey',
        edgecolor: Optional[str] = None,
        linewidth: Optional[Union[int, float]] = 1.0
    ):
        self.orientation = 'horizontal'
        self.chr_ycoordinate = None
        for chr_name, l in self.karyotype.items():
            chromosome_len, centromere_start, centromere_end = l
            verts, codes = self.__chromosome_shape(
                chromosome_len=chromosome_len,
                centromere_start=centromere_start,
                centromere_end=centromere_end,
                y_center=self.chr_coordinate[chr_name],
                centromere_angle=centromere_angle
            )
            patch = PathPatch(Path(verts, codes), facecolor=facecolor, edgecolor=edgecolor, linewidth=linewidth)
            axes.add_patch(patch)
            self.chromosome_patches[chr_name] = patch
        axes.set_yticks(list(self.chr_coordinate.values()), self.chr_coordinate.keys())
        axes.set_xlim(-1, max([i[0] for i in self.karyotype.values()]) + 1)
        axes.set_ylim(0, max(self.chr_coordinate.values()) + self.radius)
        axes.set_xlabel(self.length_unit, loc='right')
        axes.tick_params(axis='y', length=0)
        axes.spines['top'].set_color(None)
        axes.spines['right'].set_color(None)
        axes.spines['left'].set_color(None)
        return axes

    def draw_chrom(
        self,
        axes,
        centromere_angle: Optional[Union[int, float]] = 60.0,
        facecolor: Optional[str] = 'lightgrey',
        edgecolor: Optional[str] = None,
        linewidth: Optional[Union[int, float]] = 1.0
    ):
        # Vertical chromosome ideogram with a bp-scale y axis (first base = 1)
        self.orientation = 'vertical'
        r = self.radius
        space = self.space
        chr_list = list(self.karyotype.keys())
        coordinate_x = [r * (space * n - 2) for n in range(1, len(chr_list) + 1)]
        chr_xcoord = {j: i + r for i, j in zip(coordinate_x, chr_list)}
        self.chr_ycoordinate = chr_xcoord

        max_len = max([i[0] for i in self.karyotype.values()])
        for chr_name, l in self.karyotype.items():
            chromosome_len, centromere_start, centromere_end = l
            verts, codes = self.__chromosome_shape_v(
                chromosome_len=chromosome_len,
                centromere_start=centromere_start,
                centromere_end=centromere_end,
                x_center=chr_xcoord[chr_name],
                centromere_angle=centromere_angle
            )
            patch = PathPatch(Path(verts, codes), facecolor=facecolor, edgecolor=edgecolor, linewidth=linewidth)
            axes.add_patch(patch)
            self.chromosome_patches[chr_name] = patch

        # Left y-axis: chromosome length scale in bp, first base = 1
        raw_max = max(v[0] for v in self.raw_karyotype.values())
        tick_vals = matplotlib.ticker.MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10]).tick_values(1, raw_max)
        tick_vals = [t for t in tick_vals if 1 <= t <= raw_max]
        if not tick_vals or tick_vals[0] != 1:
            tick_vals = [1] + [t for t in tick_vals if t > 1]

        axes.set_yticks([t / self.scale for t in tick_vals])
        axes.set_ylabel('chromosome length (Mb)')
        axes.set_xticks(list(chr_xcoord.values()), chr_xcoord.keys(), rotation=90)
        axes.set_ylim(-max_len * 0.01, max_len * 1.01)
        axes.set_xlim(0, max(chr_xcoord.values()) + 2 * r)
        axes.tick_params(axis='x', length=0)
        axes.spines['top'].set_color(None)
        axes.spines['right'].set_color(None)
        axes.spines['bottom'].set_color(None)
        return axes

    def annotations(
        self,
        axes,
        annotations: str,
        height: Optional[Union[int, float]] = None
    ):
        if height is None:
            height = self.radius * 1.5

        groups = defaultdict(list)
        feature_color = {}
        legend_handles = []

        with open(annotations) as f:
            for line in f:
                split = line.strip().split('\t')
                if len(split) < 5:
                    continue
                chr_name, start, end, feature_type, color = split[0], split[1], split[2], split[3], split[4]
                if self.orientation == 'horizontal':
                    if chr_name not in self.chr_coordinate:
                        continue
                else:
                    if self.chr_ycoordinate is None or chr_name not in self.chr_ycoordinate:
                        continue
                start = int(start) / self.scale
                end = int(end) / self.scale
                width = end - start
                if self.orientation == 'horizontal':
                    y_center = self.chr_coordinate[chr_name]
                    y_bottom = y_center - height / 2
                    groups[(chr_name, feature_type)].append({
                        'x': start,
                        'y': y_bottom,
                        'width': width,
                        'height': height,
                        'color': color
                    })
                else:
                    x_center = self.chr_ycoordinate[chr_name]
                    x_left = x_center - height / 2
                    groups[(chr_name, feature_type)].append({
                        'x': x_left,
                        'y': start,
                        'width': height,
                        'height': width,
                        'color': color
                    })
                if feature_type not in feature_color:
                    feature_color[feature_type] = color

        for (chr_name, feature_type), rects_data in groups.items():
            rects = []
            for data in rects_data:
                rect = matplotlib.patches.Rectangle(
                    (data['x'], data['y']),
                    data['width'],
                    data['height'],
                    facecolor=data['color'],
                    edgecolor='none'
                )
                rects.append(rect)

            collection = PatchCollection(rects, match_original=True)
            chrom_patch = self.chromosome_patches.get(chr_name)
            if chrom_patch is not None:
                collection.set_clip_path(chrom_patch)

            axes.add_collection(collection)

        for feature_type, color in feature_color.items():
            proxy = matplotlib.patches.Rectangle(
                (0, 0),
                width=1,
                height=1,
                facecolor=color,
                edgecolor='none',
                label=feature_type
            )
            legend_handles.append(proxy)
        self.handles.extend(legend_handles)
        return axes

    def draw_synteny(
        self,
        axes,
        synteny_file: str,
        bezier_scale: Optional[Union[int, float]] = 0.5,
        alpha: Optional[Union[int, float]] = 0.6,
        linewidth: Optional[Union[int, float]] = 0
    ):
        # Load data
        records = []
        with open(synteny_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 8:
                    parts = line.split(',')
                if len(parts) < 8:
                    continue
                chr1, s1, e1, chr2, s2, e2, label, color = parts[0], parts[1], parts[2], parts[3], parts[4], parts[5], parts[6], parts[7]
                try:
                    s1, e1 = int(s1), int(e1)
                    s2, e2 = int(s2), int(e2)
                except ValueError:
                    continue
                records.append((chr1, s1, e1, chr2, s2, e2, label.strip(), color.strip()))

        # Define the auxiliary function of the Bézier curve
        def _bezier_curve(P0, P1, P2, P3, num=100):
            t = np.linspace(0, 1, num)
            t_ = 1 - t
            x = t_ ** 3 * P0[0] + 3 * t_ ** 2 * t * P1[0] + 3 * t_ * t ** 2 * P2[0] + t ** 3 * P3[0]
            y = t_ ** 3 * P0[1] + 3 * t_ ** 2 * t * P1[1] + 3 * t_ * t ** 2 * P2[1] + t ** 3 * P3[1]
            return np.column_stack((x, y))

        legend_handles = []
        legend_added = set()
        r = self.radius
        horizontal = self.orientation == 'horizontal'

        for chr1, start1, end1, chr2, start2, end2, label, color in records:
            if horizontal:
                if chr1 not in self.chr_coordinate or chr2 not in self.chr_coordinate:
                    continue
                c1, c2 = self.chr_coordinate[chr1], self.chr_coordinate[chr2]
                s1, e1 = start1 / self.scale, end1 / self.scale
                s2, e2 = start2 / self.scale, end2 / self.scale
            else:
                if self.chr_ycoordinate is None or chr1 not in self.chr_ycoordinate or chr2 not in self.chr_ycoordinate:
                    continue
                c1, c2 = self.chr_ycoordinate[chr1], self.chr_ycoordinate[chr2]
                s1, e1 = start1 / self.scale, end1 / self.scale
                s2, e2 = start2 / self.scale, end2 / self.scale

            if c1 > c2:
                edge1, edge2 = c1 - r, c2 + r
            else:
                edge1, edge2 = c1 + r, c2 - r

            # horizontal: segments run along x, connect on y edges; vertical: rotated
            if horizontal:
                A = np.array([s1, edge1])
                B = np.array([e1, edge1])
                C = np.array([e2, edge2])
                D = np.array([s2, edge2])
                offset = np.array([0.0, (edge2 - edge1) * bezier_scale])
            else:
                A = np.array([edge1, s1])
                B = np.array([edge1, e1])
                C = np.array([edge2, e2])
                D = np.array([edge2, s2])
                offset = np.array([(edge2 - edge1) * bezier_scale, 0.0])

            ctrl_A, ctrl_B = A + offset, B + offset
            ctrl_D, ctrl_C = D - offset, C - offset

            curve_AD = _bezier_curve(A, ctrl_A, ctrl_D, D)
            curve_BC = _bezier_curve(B, ctrl_B, ctrl_C, C)

            poly_verts = np.vstack([A, B, curve_BC, C, D, curve_AD[::-1]])
            patch = matplotlib.patches.Polygon(
                poly_verts,
                closed=True,
                facecolor=color,
                edgecolor='none' if linewidth == 0 else color,
                linewidth=linewidth,
                alpha=alpha,
                label=None,
                zorder=0
            )
            axes.add_patch(patch)

            if label not in legend_added:
                legend_added.add(label)
                proxy = matplotlib.patches.Rectangle(
                    (0, 0), width=1, height=1, facecolor=color, edgecolor='none', label=label)
                legend_handles.append(proxy)
        self.handles.extend(legend_handles)
        return axes

    def mark_loci(
        self,
        axes,
        markers: str,
        pos: Optional[str] = 'top',
        size: Optional[Union[int, float]] = None
    ):
        MARKER_SHAPES = {
            'triangle': '^',
            'circle': 'o',
            'square': 's',
            'diamond': 'D',
            'star': '*',
            'pentagon': 'p',
            'hexagon': 'h',
            'octagon': '8',
            'plus': '+',
            'cross': 'x',
            'x': 'x',
            'triangle_down': 'v',
            'triangle_left': '<',
            'triangle_right': '>',
            'thin_diamond': 'd',
            'pentagram': '*',
        }

        valid_pos_h = ('top', 'bottom')
        valid_pos_v = ('left', 'right')
        if self.orientation == 'horizontal':
            if pos not in valid_pos_h:
                raise ValueError(f"Invalid pos: {pos!r}, expected 'top' or 'bottom'")
        else:
            if pos not in valid_pos_v:
                raise ValueError(f"Invalid pos: {pos!r}, expected 'left' or 'right'")
        if size is None:
            size = self.radius * 0.5

        legend_handles = []
        legend_added = set()
        lowest = float('inf')
        highest = float('-inf')
        leftmost = float('inf')
        rightmost = float('-inf')
        with open(markers) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 6:
                    parts = line.split(',')
                if len(parts) < 6:
                    continue
                chr_name, start, end, name, color, shape = [p.strip() for p in parts[:6]]
                if self.orientation == 'horizontal':
                    if chr_name not in self.chr_coordinate:
                        continue
                else:
                    if self.chr_ycoordinate is None or chr_name not in self.chr_ycoordinate:
                        continue
                position = (int(start) + int(end)) / 2 / self.scale

                try:
                    marker_style = MarkerStyle(MARKER_SHAPES.get(shape, shape))
                    marker_path = marker_style.get_path().transformed(marker_style.get_transform())
                except ValueError:
                    raise ValueError(f"Invalid shape: {shape!r}")

                raw_vertices = marker_path.vertices
                center = (raw_vertices.min(axis=0) + raw_vertices.max(axis=0)) / 2
                extent = raw_vertices.max(axis=0) - raw_vertices.min(axis=0)
                scale_factor = size / extent.max()

                if self.orientation == 'horizontal':
                    y_center = self.chr_coordinate[chr_name]
                    if pos == 'top':
                        y_edge = y_center + self.radius
                        y = y_edge + size / 2
                        highest = max(highest, y_edge + size)
                    else:
                        y_edge = y_center - self.radius
                        y = y_edge - size / 2
                        lowest = min(lowest, y_edge - size)
                    axes.plot(
                        (position, position),
                        (y_edge, y),
                        color='#4d4d4d',
                        linewidth=1,
                        alpha=0.8
                    )
                    verts = (raw_vertices - center) * scale_factor + np.array([position, y])
                else:
                    x_center = self.chr_ycoordinate[chr_name]
                    if pos == 'left':
                        x_edge = x_center - self.radius
                        x = x_edge - size / 2
                        leftmost = min(leftmost, x_edge - size)
                    else:
                        x_edge = x_center + self.radius
                        x = x_edge + size / 2
                        rightmost = max(rightmost, x_edge + size)
                    axes.plot(
                        (x_edge, x),
                        (position, position),
                        color='#4d4d4d',
                        linewidth=1,
                        alpha=0.8
                    )
                    verts = (raw_vertices - center) * scale_factor + np.array([x, position])
                patch = matplotlib.patches.PathPatch(
                    matplotlib.path.Path(verts, marker_path.codes),
                    facecolor=color,
                    edgecolor='none',
                    zorder=3
                )
                axes.add_patch(patch)

                key = (name, color, MARKER_SHAPES.get(shape, shape))
                if key not in legend_added:
                    legend_added.add(key)
                    proxy = matplotlib.lines.Line2D(
                        [], [],
                        marker=MARKER_SHAPES.get(shape, shape),
                        linestyle='None',
                        markerfacecolor=color,
                        markeredgecolor=color,
                        markersize=12,
                        label=name
                    )
                    legend_handles.append(proxy)

        if self.orientation == 'horizontal':
            if lowest < float('inf') or highest > float('-inf'):
                ymin, ymax = axes.get_ylim()
                axes.set_ylim(min(ymin, lowest), max(ymax, highest))
        else:
            if leftmost < float('inf') or rightmost > float('-inf'):
                xmin, xmax = axes.get_xlim()
                axes.set_xlim(min(xmin, leftmost), max(xmax, rightmost))

        self.handles.extend(legend_handles)

        return axes

    def show_tracks(
        self,
        axes,
        bigwig: str,
        scale: Optional[str] = 'auto',
        max: Optional[Union[int, float]] = None,
        pos: Optional[str] = None,
        height: Optional[Union[int, float]] = None,
        color: Optional[str] = 'steelblue',
        bin_size: Optional[int] = 10000,
        reversed_chr: Optional[Dict[str, int]] = None,
        sample_name: Optional[str] = None,
        label: Optional[str] = None,
        alpha: Optional[Union[int, float]] = 0.6,
        linewidth: Optional[Union[int, float]] = 0.8,
        smooth: Optional[bool] = True,
        zorder: Optional[Union[int, float]] = 2
    ):
        if height is None:
            height = self.radius * 3 / 4
        if height <= 0:
            raise ValueError(f'Invalid height: {height!r}, expected a positive number')
        if bin_size is None:
            bin_size = 10000
        if bin_size <= 0:
            raise ValueError(f'Invalid bin_size: {bin_size!r}, expected a positive integer')
        if reversed_chr is None:
            reversed_chr = {}
        if color is None:
            color = 'steelblue'
        if label is None:
            label = sample_name if sample_name is not None else basename(bigwig)
        if max is not None:
            if max <= 0 or not np.isfinite(max):
                raise ValueError(f'Invalid max: {max!r}, expected a positive number')
            depth_max = float(max)
        else:
            depth_max = None
            if scale is None:
                scale = 'linear'
            if scale not in ('linear', 'sqrt', 'log', 'auto'):
                raise ValueError(
                    f"Invalid scale: {scale!r}, expected 'linear', 'sqrt', 'log' or 'auto'"
                )

        horizontal = self.orientation == 'horizontal'
        if pos is None:
            pos = 'top' if horizontal else 'right'
        if horizontal:
            if pos not in ('top', 'bottom'):
                raise ValueError(f"Invalid pos: {pos!r}, expected 'top' or 'bottom'")
        else:
            if pos not in ('left', 'right'):
                raise ValueError(f"Invalid pos: {pos!r}, expected 'left' or 'right'")
            if self.chr_ycoordinate is None:
                raise ValueError('draw_chrom must be called before show_tracks for vertical chromosomes')

        bw = pyBigWig.open(bigwig)
        if bw is None:
            raise ValueError(f'Unable to open bigwig file: {bigwig!r}')
        try:
            bw_chroms = bw.chroms()
            profiles = {}
            reversed_names = set()
            bw_names = {}
            for chr_name in self.karyotype:
                chr_len = self.raw_karyotype[chr_name][0]
                candidates = []
                if sample_name is not None:
                    prefix = f'{sample_name}_'
                    if chr_name.startswith(prefix):
                        candidates.append(chr_name[len(prefix):])
                    candidates.append(chr_name)
                    candidates.append(f'{sample_name}_{chr_name}')
                else:
                    candidates.append(chr_name)
                bw_name = next((name for name in candidates if name in bw_chroms), None)
                if bw_name is None:
                    continue
                bw_names[chr_name] = bw_name

                n_bins = builtins.max(1, int(np.ceil(chr_len / bin_size)))
                region_len = min(chr_len, int(bw_chroms[bw_name]))
                n_data = builtins.max(1, min(n_bins, int(np.ceil(region_len / bin_size))))
                stats = bw.stats(bw_name, 0, region_len, nBins=n_data, type='mean')
                if stats is None or len(stats) == 0:
                    continue

                data = np.array([0.0 if v is None else float(v) for v in stats])
                values = np.zeros(n_bins)
                values[:len(data)] = data[:n_bins]
                positions = (np.arange(n_bins) + 0.5) * (chr_len / n_bins)

                if smooth and n_bins >= 5:
                    first_value, last_value = values[0], values[-1]
                    padded = np.pad(values, 2, mode='edge')
                    values = np.convolve(padded, np.ones(5) / 5.0, mode='valid')
                    values[0] = first_value
                    values[-1] = last_value
                if depth_max is not None:
                    values = np.minimum(values, depth_max)

                rev_key = None
                for key in (chr_name, bw_name):
                    if key in reversed_chr:
                        rev_key = key
                        break
                if rev_key is None and sample_name is not None:
                    prefix = f'{sample_name}_'
                    if chr_name.startswith(prefix) and chr_name[len(prefix):] in reversed_chr:
                        rev_key = chr_name[len(prefix):]
                if rev_key is not None:
                    values = values[::-1]
                    positions = chr_len - positions[::-1]
                    reversed_names.add(chr_name)

                if positions[0] > 0:
                    positions = np.concatenate(([0.0], positions))
                    values = np.concatenate(([values[0]], values))
                if positions[-1] < chr_len:
                    positions = np.concatenate((positions, [chr_len]))
                    values = np.concatenate((values, [values[-1]]))

                profiles[chr_name] = (positions, values)

            if not profiles:
                raise ValueError(f'No chromosome in {bigwig!r} matches the karyotype')

            for chr_name in self.karyotype:
                candidates = [chr_name, bw_names.get(chr_name)]
                if sample_name is not None:
                    prefix = f'{sample_name}_'
                    if chr_name.startswith(prefix):
                        candidates.append(chr_name[len(prefix):])
                    candidates.append(f'{sample_name}_{chr_name}')
                if any(name in reversed_chr for name in candidates if name is not None):
                    reversed_names.add(chr_name)

            for chr_name in reversed_names:
                patch = self.chromosome_patches.get(chr_name)
                if patch is None or getattr(patch, '_reversed', False):
                    continue
                chr_len, cent_start, cent_end = self.raw_karyotype[chr_name]
                if chr_len - cent_end <= cent_start:
                    continue
                chr_len = chr_len / self.scale
                if horizontal:
                    center = (chr_len / 2.0, self.chr_coordinate[chr_name])
                else:
                    center = (self.chr_ycoordinate[chr_name], chr_len / 2.0)
                patch.set_transform(
                    Affine2D().rotate_deg_around(center[0], center[1], 180.0) + axes.transData
                )
                patch._reversed = True

            r = self.radius
            labeled = False
            track_ymin = float('inf')
            track_ymax = float('-inf')
            track_xmin = float('inf')
            track_xmax = float('-inf')
            for chr_name, (positions, values) in profiles.items():
                coords = positions / self.scale
                if depth_max is not None:
                    mode = 'linear'
                else:
                    mode = scale
                    if mode == 'auto':
                        positive = np.abs(values[values > 0])
                        if positive.size:
                            median = float(np.median(positive))
                            dynamic = float(positive.max()) / median if median > 0 else float('inf')
                        else:
                            dynamic = 1.0
                        if dynamic >= 100:
                            mode = 'log'
                        elif dynamic >= 10:
                            mode = 'sqrt'
                        else:
                            mode = 'linear'
                display_values = values
                if mode == 'log':
                    display_values = np.sign(values) * np.log1p(np.abs(values))
                elif mode == 'sqrt':
                    display_values = np.sign(values) * np.sqrt(np.abs(values))
                chr_absmax = float(np.abs(display_values).max())
                depth_scale = height / chr_absmax if chr_absmax > 0 else 0.0
                curve = np.clip(display_values * depth_scale, -height, height)
                if horizontal:
                    if pos == 'top':
                        base = self.chr_coordinate[chr_name] + r
                        sign = 1.0
                    else:
                        base = self.chr_coordinate[chr_name] - r
                        sign = -1.0
                    y = base + sign * curve
                    axes.fill_between(
                        coords,
                        base,
                        y,
                        facecolor=color,
                        edgecolor='none',
                        alpha=alpha,
                        zorder=zorder,
                        label=None if labeled else label
                    )
                    if linewidth:
                        axes.plot(coords, y, color=color, linewidth=linewidth, zorder=zorder)
                    track_ymin = min(track_ymin, float(y.min()))
                    track_ymax = builtins.max(track_ymax, float(y.max()))
                else:
                    if pos == 'right':
                        base = self.chr_ycoordinate[chr_name] + r
                        sign = 1.0
                    else:
                        base = self.chr_ycoordinate[chr_name] - r
                        sign = -1.0
                    x = base + sign * curve
                    axes.fill_betweenx(
                        coords,
                        base,
                        x,
                        facecolor=color,
                        edgecolor='none',
                        alpha=alpha,
                        zorder=zorder,
                        label=None if labeled else label
                    )
                    if linewidth:
                        axes.plot(x, coords, color=color, linewidth=linewidth, zorder=zorder)
                    track_xmin = min(track_xmin, float(x.min()))
                    track_xmax = builtins.max(track_xmax, float(x.max()))
                if not labeled:
                    proxy = matplotlib.patches.Rectangle(
                        (0, 0),
                        width=1,
                        height=1,
                        facecolor=color,
                        edgecolor='none',
                        alpha=alpha,
                        label=label
                    )
                    self.handles.append(proxy)
                    labeled = True

            if horizontal:
                ymin, ymax = axes.get_ylim()
                if track_ymin < float('inf'):
                    ymin = min(ymin, track_ymin)
                if track_ymax > float('-inf'):
                    ymax = builtins.max(ymax, track_ymax)
                axes.set_ylim(ymin, ymax)
            else:
                xmin, xmax = axes.get_xlim()
                if track_xmin < float('inf'):
                    xmin = min(xmin, track_xmin)
                if track_xmax > float('-inf'):
                    xmax = builtins.max(xmax, track_xmax)
                axes.set_xlim(xmin, xmax)
        finally:
            bw.close()

        return axes
