from typing import Union, Dict, List
from collections import defaultdict
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.collections import PatchCollection


class ChromosomeIdeogram:
    def __init__(
        self,
        karyotype: Union[str, Dict[str, List[int]]],
        radius=0.4,
        space: int = 3
    ):
        if isinstance(karyotype, dict):
            raw_karyotype = karyotype
        else:
            raw_karyotype = {}
            with open(karyotype) as f:
                for line in f:
                    split = line.strip().split('\t')
                    chr_name, chr_len, cent_start, cent_end = split
                    chr_len, cent_start, cent_end = int(chr_len), int(cent_start), int(cent_end)
                    raw_karyotype[chr_name] = [chr_len, cent_start, cent_end]
        longest_chr = max([i[0] for i in raw_karyotype.values()])
        if longest_chr / 10 ** 3 <= 100:
            scale = 10 ** 3
            length_unit = 'Kb'
        elif longest_chr / 10 ** 6 <= 100:
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
        color: str = 'lightgrey'
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
            patch = PathPatch(Path(verts, codes), facecolor=color, edgecolor='w', linewidth=0)
            axes.add_patch(patch)
            self.chromosome_patches[chr_name] = patch
        axes.set_yticks(list(self.chr_coordinate.values()), self.chr_coordinate.keys())
        axes.set_xlim(-1, max([i[0] for i in self.karyotype.values()]) + 1)
        axes.set_ylim(0, max(self.chr_coordinate.values()) + self.radius)
        axes.set_xlabel(self.length_unit, loc='right')
        plt.tick_params(axis='y', length=0)
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
                plt.legend()
