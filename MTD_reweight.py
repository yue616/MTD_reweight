import numpy as np
import argparse
from scipy.ndimage import gaussian_filter


class Reweight:

    def __init__(self, rcv: list, bias: list, colvar="COLVAR",
                 temp=298.15, nbin=100, min=None, max=None):
        """
        :param rcv: list. Column(s) in COLVAR file to be reweighted, note that the first column starts from 1 rather than 0.
        :param bias: list. Column(s) in COLVAR file contains bias potential, e.g., MTD bias, wall bias, etc., note
                           that the first column starts from 1.
        :param colvar: str. Filename of colvar file, COLVAR is used by default.
        :param temp: float. System temperature.
        :param nbin: int. The number of bins for the outputted FES.
        :param min: float. the minimum of Grids boundaries, program will detect it if omitted.
        :param max: float. the maximum of Grids boundaries, program will detect it if omitted.
        """
        self.rcv = rcv
        self.bias = bias
        try:
            # read colvar file
            self.colvar = np.loadtxt(colvar)
        except FileNotFoundError:
            raise

        self.kbt = temp * 8.31446261815324e-3  # Boltzmann constant 8.31446261815324e-3 kj/mol/k
        self.nbin = nbin
        self.min = []
        self.max = []

        # find the border of CVs if min and max are not specified.
        if min is None and max is None:
            for cvi in self.rcv:
                col_idx = cvi - 1
                min, max = np.min(self.colvar[:, col_idx]), np.max(self.colvar[:, col_idx])
                self.min.append(min)
                self.max.append(max)
        elif min is None and max is not None:
            for cvi in self.rcv:
                col_idx = cvi - 1
                min = np.min(self.colvar[:, col_idx])
                self.min.append(min)
            self.max.extend(max)
        elif max is None and min is not None:
            for cvi in self.rcv:
                col_idx = cvi - 1
                max = np.max(self.colvar[:, col_idx])
                self.max.append(max)
            self.min.extend(min)
        else:
            self.max.extend(max)
            self.min.extend(min)

        assert len(self.rcv) == len(self.max) == len(self.min)

    def reweight(self) -> tuple:

        # initialize grid for storing bias value.
        grids = []
        for i in range(len(self.rcv)):
            grid = np.linspace(self.min[i], self.max[i], num=self.nbin, endpoint=True)
            grids.append(grid)

        # initialize FES grids
        fes = np.zeros([self.nbin] * len(self.rcv))

        # loop row in COLVAR and accumulate the calculated weights
        for row in self.colvar:
            # allocate the CVs into grids and find CVs indices
            diff = [np.abs(row[j - 1] - grids[i]) for i, j in enumerate(self.rcv)]
            cvs_loc = [np.argmin(i) for i in diff]
            cvs_loc_idx = tuple(cvs_loc)

            # calculate the weights through the formulae w ~ exp(V(s)/kbt), where the V(s) is bias value
            bias = 0.0
            for i in self.bias:
                col_idx = i - 1
                bias += row[col_idx]

            w = np.exp(bias / self.kbt)

            # accumulate weights into fes grids
            fes[cvs_loc_idx] += w

        # FES is calculated by -KbT*ln(p), where the kj/mol is used by default
        # ignore the divide zero error
        np.seterr(divide="ignore")
        fes = - self.kbt * np.log(fes / np.sum(fes))

        return grids, fes

    def output_fes(self, grids, fes):
        """
        for three CVs, output FES in gaussian-type cube format and also a plain text file.
        for two CVs, output format is suitable for gnuplot.
        for one CV, first column corresponding to grid value, and second column for free energy.
        :return: output grids and FES when CVs higher than three
        """

        if len(grids) > 3:
            print("warning! Can not output fes when CVs number higher than 3, return FES and corresponding grids")
            return grids, fes
        elif len(grids) == 3:
            # firstly, output fes into a plain text file
            out = []
            out.extend(["#Fields", "\t", "CV1", "\t", "CV2", "\t", "CV3", "\t", "free energy", "\n"])
            for i, iv in enumerate(grids[0]):
                for j, jv in enumerate(grids[1]):
                    for k, kv in enumerate(grids[2]):
                        out.extend(
                            [f"{iv:.7f}", "\t", f"{jv:.7f}", "\t", f"{kv:.7f}", "\t", f"{fes[i, j, k]:.7f}", "\n"])

            with open("fes_3d.txt", "w") as f:
                f.write("".join(out))

            # secondly, output fes into a gaussian-type cube format.
            grids_num_tot = self.nbin ** 3
            grids_vec1 = (grids[0].max() - grids[0].min()) / len(grids[0])
            grids_vec2 = (grids[1].max() - grids[1].min()) / len(grids[1])
            grids_vec3 = (grids[2].max() - grids[2].min()) / len(grids[2])

            cube = []
            cube.extend(["FES in Gaussian-type cube format", "\n"])
            cube.extend([f"Total {grids_num_tot} grids", "\n"])
            cube.extend(["  -1", "\t", f"{0:.6f}", "\t", f"{0:.6f}", "\t", f"{0:.6f}", "\n"])
            cube.extend([f"  {self.nbin}", "\t", f"{grids_vec1:.6f}", "\t", f"{0:.6f}", "\t", f"{0:.6f}", "\n"])
            cube.extend([f"  {self.nbin}", "\t", f"{0:.6f}", "\t", f"{grids_vec2:.6f}", "\t", f"{0:.6f}", "\n"])
            cube.extend([f"  {self.nbin}", "\t", f"{0:.6f}", "\t", f"{0:.6f}", "\t", f"{grids_vec3:.6f}", "\n"])
            cube.extend(["  1", "\t", "1", "\n"])

            for i, _ in enumerate(grids[0]):
                for j, _ in enumerate(grids[1]):
                    for k, _ in enumerate(grids[2]):
                        if (k + 1) % 6 == 0 or k == len(grids[2]) - 1:
                            cube.extend([f"{fes[i, j, k]:.5e}", "\n"])
                        else:
                            cube.extend([f"{fes[i, j, k]:.5e}", "\t"])

            with open("fes.cube", "w") as f:
                f.write("".join(cube))

        elif len(grids) == 2:
            out = []
            out.extend(["#Fields", "\t", "CV1", "\t", "CV2", "\t", "free energy", "\n"])
            for i, iv in enumerate(grids[0]):
                for j, jv in enumerate(grids[1]):
                    if j != len(grids[1]) - 1:
                        out.extend([f"{iv:.7f}", "\t", f"{jv:.7f}", "\t", f"{fes[i, j]:.7f}", "\n"])
                    elif j == len(grids[1]) - 1:
                        out.extend([f"{iv:.7f}", "\t", f"{jv:.7f}", "\t", f"{fes[i, j]:.7f}", "\n", "\n"])

            with open("fes_2d.txt", "w") as f:
                f.write("".join(out))

        elif len(grids) == 1:
            out = []
            out.extend(["#Fields", "\t", "CV1", "\t", "free energy", "\n"])
            for i, iv in enumerate(grids[0]):
                out.extend([f"{iv:.7f}", "\t", f"{fes[i]:.7f}", "\n"])

            with open("fes_1d.txt", "w") as f:
                f.write("".join(out))

    def smoothing(self, fes, sigma):

        smoothed_fes = gaussian_filter(fes, sigma)

        return smoothed_fes


def parse_args():
    parser = argparse.ArgumentParser(description='Free energy reconstruction (reweight) based on the arthrogram '
                                                 'proposed by Pratyush Tiwary and Michele Parrinello '
                                                 '(J. Phys. Chem. B 2015, 119, 736−742)')
    parser.add_argument('--rcv', dest='rcv', type=int, nargs='+', help='Column(s) in COLVAR file to be '
                                                'reweighted, note that the first column starts from 1 rather than 0.')
    parser.add_argument('--bias', dest='bias', type=int, nargs='+', help='Column(s) of bias potential in COLVAR file. '
                                                                'These biases can be, e.g., MTD bias, wall bias, etc.')
    parser.add_argument('--colvar', dest='colvar', type=str, help='Filename of colvar file, '
                                                                  'COLVAR is used by default.', default="COLVAR")
    parser.add_argument('--temp', dest='temp', type=float, help='System temperature.', default=298.15)
    parser.add_argument('--nbin', dest='nbin', type=int, help='The number of bins for the outputted FES', default=100)
    parser.add_argument('--min', dest='min', type=float, nargs='+', help='The minimum of Grids boundaries, program will '
                                                              'detect it if omitted.', default=None)
    parser.add_argument('--max', dest='max', type=float, nargs='+', help='The maximum of Grids boundaries, program will '
                                                              'detect it if omitted.', default=None)
    parser.add_argument('--smooth', dest='smooth', action='store_true', help='Performing Gaussian smoothing for the '
                                                                            'calculated Free Energy Surface (FES).')
    parser.add_argument('--sigma', dest='sigma', type=float, nargs='+', help='The sigma parameter used in a Gaussian '
                       'function determines the level of smoothing applied to the fes. A larger sigma value results in '
                       'a smoother fes. Singe sigma value means the level of smoothing is equalvent on all grid dimention, '
                       'Note that if multiple sigma are provided, each sigma will determine the level of smoothing on '
                       'corresponding grid dimention, total sigma should equal to the number of reweighted CVs.', default=0.1)

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    r = Reweight(rcv=args.rcv, bias=args.bias, colvar=args.colvar, temp=args.temp, nbin=args.nbin, min=args.min,
                 max=args.max)
    grids, fes = r.reweight()
    if args.smooth:
        fes = r.smoothing(fes, args.sigma)
    r.output_fes(grids=grids, fes=fes)


