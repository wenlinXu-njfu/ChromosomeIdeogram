from typing import Union, Dict, List
from collections import defaultdict
from os import makedirs, getcwd
from os.path import basename
import numpy as np
import matplotlib
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.collections import PatchCollection
from click import echo


def data_preprocessing(
    karyotype: str,
    annotation: Union[str, None] = None,
    synteny: Union[str, None] = None,
    sample_name: Union[str, None] = None,
    exclude_chr: Union[str, List[str], None] = None,
    force_revers_chr: Union[str, List[str], None] = None,
    out_path: Union[str, None] = None,
) -> None:
    if exclude_chr is None:
        exclude_chr = []
    if force_revers_chr is None:
        force_revers_chr = []
    if out_path is None:
        out_path = getcwd()
    makedirs(out_path, exist_ok=True)

    d = {}
    with open(karyotype) as f, open(f'{out_path}/{basename(karyotype)}.processed', 'w') as o:
        for line in f:
            if not line.startswith('#') and line.strip():
                split = line.strip().split('\t')
                chr_name, chr_len, cent_start, cent_end = split[0], int(split[1]), int(split[2]), int(split[3])
                if chr_len - cent_end > cent_start and chr_name not in exclude_chr:  # short arm towards the right
                    split[2] = str(chr_len - cent_end)
                    split[3] = str(chr_len - cent_start)
                    d[split[0]] = chr_len
                elif chr_name in force_revers_chr:
                    split[2] = str(chr_len - cent_end)
                    split[3] = str(chr_len - cent_start)
                    d[split[0]] = chr_len
                if sample_name is not None:
                    split[0] = f'{sample_name}_{chr_name}'
                echo('\t'.join(split), o)
    echo(f'Reversed chromosome: {" ".join(d.keys())}', err=True)

    if annotation is not None:
        with open(annotation) as f, open(f'{out_path}/{basename(annotation)}.processed', 'w') as o:
            for line in f:
                if not line.startswith('#') and line.strip():
                    split = line.strip().split('\t')
                    chr_name, start, end = split[0], int(split[1]), int(split[2])
                    try:
                        split[1] = str(d[chr_name] - end)
                        split[2] = str(d[chr_name] - start)
                    except KeyError:
                        pass
                    if sample_name is not None:
                        split[0] = f'{sample_name}_{split[0]}'
                    echo('\t'.join(split), o)

    if synteny is not None:
        with open(synteny) as f, open(f'{out_path}/{basename(synteny)}.processed', 'w') as o:
            for line in f:
                if not line.startswith('#') and line.strip():
                    split = line.strip().split('\t')
                    chr1, start1, end1 = split[0], int(split[1]), int(split[2])
                    chr2, start2, end2 = split[3], int(split[4]), int(split[5])
                try:
                    split[1] = str(d[chr1] - end1)
                    split[2] = str(d[chr1] - start1)
                except KeyError:
                    pass
                try:
                    split[4] = str(d[chr2] - end2)
                    split[5] = str(d[chr2] - start2)
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
        radius=0.4,
        space: int = 4
    ):
        if isinstance(karyotype, dict):
            raw_karyotype = karyotype  # {chromosome: [length, centromere_start, centromere_end]}
        else:
            raw_karyotype = {}
            with open(karyotype) as f:
                for line in f:
                    split = line.strip().split('\t')
                    chr_name, chr_len, cent_start, cent_end = split
                    chr_len, cent_start, cent_end = int(chr_len), int(cent_start), int(cent_end)
                    raw_karyotype[chr_name] = [chr_len, cent_start, cent_end]
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

        self.karyotype = scaled_karyotype
        self.scale = scale
        self.length_unit = length_unit
        self.chr_coordinate = chr_coordinate
        self.chromosome_patches = {}
        self.radius = radius

    def __chromosome_shape(
        self,
        chromosome_len: int,
        centromere_start: int,
        centromere_end: int,
        y_center: Union[int, float] = 0,
        centromere_angle: Union[int, float] = 60
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

    def draw_chromh(
        self,
        axes: matplotlib.axes.Axes,
        centromere_angle: Union[int, float] = 60.0,
        facecolor: str = 'lightgrey',
        edgecolor: str = None,
        linewidth: float = 1.0
    ) -> matplotlib.axes.Axes:
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

    def annotations(
        self,
        axes: matplotlib.axes.Axes,
        annotations: str
    ):
        groups = defaultdict(list)
        feature_color = {}

        with open(annotations) as f:
            for line in f:
                split = line.strip().split('\t')
                if len(split) < 5:
                    continue
                chr_name, start, end, feature_type, color = split
                start = int(start) / self.scale
                end = int(end) / self.scale
                width = end - start
                y_center = self.chr_coordinate[chr_name]
                height = 0.7
                y_bottom = y_center - height / 2

                groups[(chr_name, feature_type)].append({
                    'x': start,
                    'y': y_bottom,
                    'width': width,
                    'height': height,
                    'color': color
                })
                if feature_type not in feature_color:
                    feature_color[feature_type] = color

        legend_added = set()
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

            if feature_type not in legend_added:
                legend_added.add(feature_type)
                proxy = matplotlib.patches.Rectangle(
                    (max([i[0] for i in self.karyotype.values()]) * 10000, 0), 1, 1,
                    facecolor=feature_color[feature_type],
                    edgecolor='none',
                    label=feature_type
                )
                axes.add_patch(proxy)

    def draw_synteny(
        self,
        axes: matplotlib.axes.Axes,
        synteny_file: str,
        bezier_scale: float = 0.5,
        alpha: float = 0.6,
        linewidth: float = 0
    ) -> matplotlib.axes.Axes:
        # load data
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
                chr1, s1, e1, chr2, s2, e2, label, color = parts[:8]
                try:
                    s1, e1 = int(s1), int(e1)
                    s2, e2 = int(s2), int(e2)
                except ValueError:
                    continue
                records.append((chr1, s1, e1, chr2, s2, e2, label.strip(), color.strip()))

        # Cubic Bezier curve auxiliary function
        def _bezier_curve(P0, P1, P2, P3, num=100):
            t = np.linspace(0, 1, num)
            t_ = 1 - t
            x = t_ ** 3 * P0[0] + 3 * t_ ** 2 * t * P1[0] + 3 * t_ * t ** 2 * P2[0] + t ** 3 * P3[0]
            y = t_ ** 3 * P0[1] + 3 * t_ ** 2 * t * P1[1] + 3 * t_ * t ** 2 * P2[1] + t ** 3 * P3[1]
            return np.column_stack((x, y))

        legend_added = set()
        proxy_x = max([i[0] for i in self.karyotype.values()]) * 10000

        for chr1, start1, end1, chr2, start2, end2, label, color in records:
            if chr1 not in self.chr_coordinate or chr2 not in self.chr_coordinate:
                continue

            x1s, x1e = start1 / self.scale, end1 / self.scale
            x2s, x2e = start2 / self.scale, end2 / self.scale
            y1 = self.chr_coordinate[chr1]
            y2 = self.chr_coordinate[chr2]

            if y1 > y2:
                y1_edge, y2_edge = y1 - self.radius, y2 + self.radius
            else:
                y1_edge, y2_edge = y1 + self.radius, y2 - self.radius

            dy = (y2_edge - y1_edge) * bezier_scale
            A = np.array([x1s, y1_edge])
            B = np.array([x1e, y1_edge])
            C = np.array([x2e, y2_edge])
            D = np.array([x2s, y2_edge])

            ctrl_A, ctrl_B = A + [0, dy], B + [0, dy]
            ctrl_D, ctrl_C = D - [0, dy], C - [0, dy]

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
                    xy=(proxy_x, 0),
                    width=1,
                    height=1,
                    facecolor=color,
                    edgecolor='none',
                    label=label
                )
                axes.add_patch(proxy)

        return axes