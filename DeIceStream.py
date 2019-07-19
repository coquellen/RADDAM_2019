import os
import sys
if sys.version_info < (2,7):
        raise RuntimeError("minimum Python 2.7.0 is required!")
import numpy as np
import argparse

def cos2(a):
    return np.cos(a) * np.cos(a)

def sin2(a):
    return np.sin(a) * np.sin(a)


def resolution(scell, hkl):
    """

    :param scell: 6-tuple - will be internally converted to float
    :param hkl: numpy array of shape (N,3)
    :return:  numpy array of shape (N,)
    """

    a = float(scell[0])
    b = float(scell[1])
    c = float(scell[2])
    al = np.deg2rad(float(scell[3]))
    be = np.deg2rad(float(scell[4]))
    ga = np.deg2rad(float(scell[5]))  # in degrees
    h = hkl[:, 0]
    k = hkl[:, 1]
    l = hkl[:, 2]
    pf = 1.0 - cos2(al) - cos2(be) - cos2(ga) + 2.0 * np.cos(al) * np.cos(be) * np.cos(ga)
    n1 = h * h * sin2(al) / (a * a) + k * k * sin2(be) / (b * b) + l * l * sin2(ga) / (c * c)
    n2a = 2.0 * k * l * (np.cos(be) * np.cos(ga) - np.cos(al)) / (b * c)
    n2b = 2.0 * l * h * (np.cos(ga) * np.cos(al) - np.cos(be)) / (c * a)
    n2c = 2.0 * h * k * (np.cos(al) * np.cos(be) - np.cos(ga)) / (a * b)

    return np.sqrt((n1 + n2a + n2b + n2c) / pf)





class Stream(object):

    def __init__(self, streamfile):
        self.streamfile = streamfile

        self.header, self.frames, self.geom, self.all_reflections = self.parse_stream(self.streamfile)
        self.all_reflections = np.array(self.all_reflections)

    def parse_stream(self,streamfile):
        append = 0
        crystal = 0
        count_shots = 0
        count_crystals = 0

        frame_stream = []
        frames = []
        self.cells = []
        head_check = 0
        header = ''
        geom = Geom()
        stream = open(streamfile)
        all_reflections = []
        REFLECTIONS = 0
        GEOM = 0

        for index, line in enumerate(stream):

            # Get header from stream
            while (head_check == 0):

                header += line
                if '----- Begin geometry file -----' in line: GEOM = 1

                if GEOM:
                    if 'clen' in line:
                        geom.distance = float(line.split('=')[1].split(';')[0].strip())
                    if 'corner_x' in line:
                        geom.bx = -1 * float(line.split('=')[1].split(';')[0].strip())
                    if 'corner_y' in line:
                        geom.by = -1 * float(line.split('=')[1].split(';')[0].strip())
                    if 'res' in line:
                        geom.ps = 1 / float(line.split('=')[1].split(';')[0].strip())
                    if 'photon_energy' in line:
                        geom.energy = float(line.split('=')[1].split(';')[0].strip())
                        geom.wl = 12398.425 / geom.energy
                    if 'max_fs' in line:
                        geom.max_fs = float(line.split('=')[1].split(';')[0].strip())
                    if 'max_ss' in line:
                        geom.max_ss = float(line.split('=')[1].split(';')[0].strip())

                break

            if 'Begin chunk' in line:
                count_shots += 1
                append = 1

            elif 'Begin crystal' in line:
                crystal = 1
                count_crystals += 1
                frame = Frame()
            elif "Cell parameters" in line:
                ls = line.split()
                a, b, c = map(float, ls[2:5])
                al, be, ga = map(float, ls[6:9])
                self.cells.append([a, b, c, al, be, ga])

            elif "----- End unit cell -----" in line:
                head_check = 1

            elif 'End of reflections' in line:
                REFLECTIONS = 0

            elif "End chunk" in line:
                if crystal == 1:
                    frame.all = frame_stream
                    frames.append(frame)
                append = 0
                frame_stream = []
                crystal = 0

            if REFLECTIONS:
                ls = line.split()
                h, k, l = map(int, ls[0:3])
                I, SIG, P, BKG, X, Y = map(float, ls[3:9])
                all_reflections.append([h, k, l, I, SIG, BKG, X, Y])
                frame.hkl_stream.append(line)

            if append == 1: frame_stream.append(line)

            if '  h    k    l          I   sigma(I)' in line:
                REFLECTIONS = 1
                append = 0

            if count_shots % 100 == 0: print '%7i frames parsed, %7i crystals found\r' % (count_shots, count_crystals),
            sys.stdout.flush()

        print '%7i frames parsed, %7i crystals found\r' % (count_shots, count_crystals),
        sys.stdout.flush()
        print

        return header, frames, geom, all_reflections

    def bkg_stats(self,maxres, bins, sig, outputfn):
        
        self.bins = 1./ (1. /maxres * np.arange(0.0001, 1.05, 1. / bins) ** (1. / 3.))[::-1]

        hkl = self.all_reflections[:,0:3]
        cell = np.mean(np.array(self.cells), axis=0)
        RES = resolution(cell, hkl)
        BKG = self.all_reflections[:,5]
        #RES = self.geom.get_resolution(self.all_reflections[:, 6], self.all_reflections[:, 7])

        # Calculate the mean and std bkg values for each resolution bin
        TOFIX = []
        mean = []
        std = []
        for i in range(self.bins.size):
            # Get resolution range limits
            if i == self.bins.size - 1:
                rmax = RES.max()
            else:
                rmax = self.bins[i + 1]
            rmin = self.bins[i]

            #Get bkg values for the considered bin
            idx = np.where( (RES < rmax) & ( RES >= rmin))
            IR0 = BKG[idx]
            m = IR0.mean()
            s = IR0.std()
            mean.append(m)
            std.append(s)
            message=""
            #Worrisome bins are flagged if the mean bkg level is if std > mean
            if s > m:
                message = " -- WORRISOME !! -- "
                TOFIX.append(True)
            else:
                TOFIX.append(False)
            print("Shell %4.2f - %4.2f : %4.2f +/- %4.2f %s" % (rmax, rmin, m, s, message))
        print("/n/n")

        #Now flagging the outliers

        #The threshold to reject outliers reflections is computed using thefor worrisome resolutions bins from the closest regular bin and is equal to mean(bkg) + 2 * std(bkg)
        threshold_idx = np.where(np.array(TOFIX) == True)[0][0]
        threshold = mean[threshold_idx] + sig * std[threshold_idx]

        IDX = []
        Nremoved = 0
        NWorrBins = 0
        for i in range(len(TOFIX)):
            if i == self.bins.size - 1:
                rmax = 50
            else:
                rmax = self.bins[i + 1]
            rmin = self.bins[i]

            # If the bin is OK, just updating the threshold value
            if not TOFIX[i]:
                threshold = mean[i] + sig * std[i]

            # If the bin is worrisome
            if TOFIX[i]:

                #getting the reflections in the bins below the threshold (to compute the new mean and sd after outliers rejection)
                idx = np.where(((RES <= rmax) & (RES > rmin)) & (BKG <= threshold))

                m = BKG[idx[0]].mean()
                s = BKG[idx[0]].std()
                #threshold = mean[i] + sig * std[i]
                
                #getting the outliers
                tmp = list(np.where(((RES <= rmax) & (RES > rmin)) & (BKG >= threshold))[0])
                Nremoved += len(tmp)
                NWorrBins += len(tmp) + idx[0].size
                IDX.extend(tmp)
                print("Shell %4.2f - %4.2f : %4.2f +/- %4.2f -- %8i outliers removed out of %8i reflections, i.e (%4.2f %%)" % (rmax, rmin, m, s, len(tmp), len(tmp) + idx[0].size, float(len(tmp)) / (len(tmp) + idx[0].size) * 100))
            else:
                print("Shell %4.2f - %4.2f : %4.2f +/- %4.2f" % (rmax, rmin, mean[i], std[i]))
        print("In the worrisome resolution bins, %8i outliers out of %8i reflections will be deleted from the stream (%4.2f%%)"%(Nremoved, NWorrBins, float(Nremoved) / NWorrBins * 100))

        IDX = sorted(list(set(IDX)))
        self.count = 0
        self.idx = 0


        out = open(outputfn, 'w')
        print >> out, self.header,
        j = 0
        for frame in self.frames:
            toberemoved = []
            while j < len(IDX) and IDX[j] < len(frame.hkl_stream) + self.count:
                    toberemoved.append(IDX[j] - self.count)
                    j += 1
            self.count += len(frame.hkl_stream)
            for i in sorted(toberemoved, reverse=True):
                    del frame.hkl_stream[i]
            print >> out, ''.join(frame.all),
            print >> out, ''.join(frame.hkl_stream),
            print >> out, "End of reflections\n--- End crystal\n----- End chunk -----"
        out.close()

class Frame(object):
    def __init__(self):
        self.a = 0
        self.b = 0
        self.c = 0
        self.alpha = 90.
        self.beta = 90.
        self.gamma = 90.
        self.res = 5.
        self.index = 'none'
        self.filename = 'example.h5'
        self.timeline = 0
        self.indexing = ''
        self.astar = np.zeros((3,))
        self.bstar = np.zeros((3,))
        self.cstar = np.zeros((3,))
        self.reflections = []
        self.hkl_stream = []


class Geom(object):
    def __init__(self):
        self.distance = 0.1376
        self.energy = 13450
        self.wl = 0.92
        self.bx = 1131.08
        self.by = 1192.0
        self.ps = 7.5e-5
        self.RADIUS = 0
        self.RES = 0
        self.max_fs = 2069
        self.max_ss = 2166

    def get_resolution(self, x, y):
        dx = (x - self.bx) * self.ps
        dy = (y - self.by) * self.ps
        radius = np.sqrt(dx ** 2 + dy ** 2)
        theta = 0.5 * np.arctan(radius / self.distance)
        return self.wl / (2. * np.sin(theta))

    def get_maxres(self):
        #TODO: needs to be adapted for CSPAD or MPCCD
        x1, y1 = 0, 0
        d1 = np.sqrt( (x1 - self.bx ) ** 2 + (y1 - self.by) ** 2)
        x2, y2 = self.max_fs, self.max_ss
        d2 = np.sqrt((x2 - self.bx) ** 2 + (y2 - self.by) ** 2)
        x3, y3 = self.max_fs, 0
        d3 = np.sqrt((x3 - self.bx) ** 2 + (y3 - self.by) ** 2)
        x4, y4 = 0, self.max_ss
        d4 = np.sqrt((x4 - self.bx) ** 2 + (y4 - self.by) ** 2)

        d = max(d1, d2, d3, d4)
        theta = 0.5 * np.arctan( d * self.ps / self.distance)
        return (d, self.wl / (2. * np.sin(theta)))

    def get_bins(self, maxres, nbins):
        RADIUS, RES = self.get_maxres()
        if maxres < RES:
            print("Sorry. Given the geometry, highest resolution is %4.2f"%RES)
            maxres = RES
            RMAX = RADIUS
        else:
            theta =  np.arcsin(self.wl / (2 * maxres))
            RMAX = np.tan(2 * theta) * self.distance / self.ps
        S = np.pi * RMAX ** 2 / float(nbins)
        radius = np.sqrt(np.arange(1, nbins+1, 1) * S / np.pi) * self.ps
        theta = 0.5 * np.arctan(radius / self.distance)
        self.bins = self.wl / (2. * np.sin(theta))
        return self.bins

    
    
class ArgumentParser(argparse.ArgumentParser):
    def __init__(self):
        desc = """
        """
        argparse.ArgumentParser.__init__(self, description=desc)
        self.add_argument("-i","--input", nargs=1, type=str, default= None,
                       help="Input stream")
        self.add_argument("-r", "--resolution","--res", nargs=1, type=float, default= None,
                          help="High resolution limit (default: High resolution at detector edge)")
        self.add_argument("-n", "--bins", nargs=1, default= [50], type=int,
                       help="Number of bins to be used (default: 50)")

        self.add_argument("-s","--sigma", nargs= 1, type=int, default= [2],
                       help="Sigma outlier rejection (default:2)")

        self.add_argument("-o","--output", nargs= 1, type=str, default= None,
                       help="Optional - Output stream")

        self.args = self.parse_args()
        self.checkInput()
        
    def checkInput(self):
        
        if self.args.resolution == None:
            return False
        else:
            self.res = self.args.resolution[0]
        
        if self.args.input == None:
            return False
        
        if not os.path.isfile(self.args.input[0]):
            print("%s: No such crystFEL stream\n" % self.args.input[0])
            return False
        else: 
            self.stream = self.args.input[0]
            if self.args.output == None:
                filename, file_extension = os.path.splitext(self.stream)
                self.output = '%s_OutliersRemoved.%s' % (filename, file_extension)
            else:
                self.output = self.args.output[0]
            
        self.nbins = self.args.bins[0]
        self.sig = self.args.sigma[0]
        
        return True
        

if __name__ == '__main__':
    
    
    parser = ArgumentParser()

    if parser.checkInput():
        print("Input CrystFEL stream   : %s" % parser.stream)
        print("Output CrystFEL stream  : %s" % parser.output)
        print("N resolutions bins      : %i" % parser.nbins)
        print("High resolution cutoff  : %4.2f" % parser.res)
        print("Sigma outlier rejection : %i" % parser.sig)

        s = Stream(parser.stream)
        s.bkg_stats(parser.res, parser.nbins, parser.sig, parser.output)

    else:
        parser.print_help()


