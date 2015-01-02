# -*- coding: utf-8; -*-

"""
Originally created on May 4, 2014. Enhanced on July 23, 2014.
Part of the script was adapted from
http://verahill.blogspot.com/2013/09/514-extracting-data-form-pes-scan-with.html
Uses Python 2.7 and libraries as implemented in Anaconda from Contibuum Analytics

Run from the terminal window (cmd) or shell as:
>> python pes_nbo3.py output_file.out

Requires Gaussian PES output file (output_file.out) from the Gaussian PES job.
Examples of such files are part of the download in the GitHub repo.
"""

# author:   'Marcel Patek'
# filename: 'test.py'
# date:      7/23/2014
# version:  '1.1'
# email:    'chemgplus@gmail.com'
# license:  'GNU3'
# usage:     python pes_nbo3.py output_file.out

'''
 * Copyright (C) 2014 Marcel Patek
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * For a copy of the GNU General Public License,
 * see <http://www.gnu.org/licenses/>.
'''

import sys
import os
import re


def print_frame_top(m, c):
    print(c * m)


def print_frame_bot(m, c):
    print(c * m)


def main(argv):
    if len(argv) < 2:
        print_frame_top(60, '+')
        sys.stderr.write("\n   Usage: >>[python] %s gau_output.out\n\n" % (argv[0],))
        sys.stderr.write(
            "\n   or  : >>[python] %sy gau_output.out\n\n" % sys.argv[0].split("\\")[len(sys.argv[0].split("\\")) - 1][
                0:-1])
        print_frame_bot(60, '+')
        return 1

    if not os.path.exists(argv[1]):
        print_frame_top(60, '+')
        sys.stderr.write("\n   ERROR: *.out file %r was not found!\n\n" % (argv[1],))
        print_frame_bot(60, '+')
        return 1

    if len(getscan(sys.argv[1])) < 1:
        print_frame_top(60, '+')
        sys.stderr.write("\n ERROR: This does not seem to be the right file. Scan coordinate is missing.\n\n")
        print_frame_bot(60, '+')
        return 1


def rundif(it, zero):
    """
    Create values for relative energy in kcal/mol
    :param it: list of energies
    :param zero: energy to which other values will be referenced to
    """
    for x in it:
        ener = x - zero
        yield ener * 627.51


def getscan(infile):
    """
    Find the 'Scan' keyword in gau.out file
    :param infile: Gaussian PES output file
    :return: line as string containing word Scan
    """
    try:
        f = open(infile, 'r')
        getcoord = ''
        fi = f.readlines()
        for line in fi:
            if '!' in line and "Scan" in line:
                getcoord = line.split()  # splits by words
                f.close()
        return getcoord
    except IOError:
        print "This does not seem to be the right file. Scan coordinate is missing."


def getrawdata(infile):
    f = open(infile, 'r')
    opt = 0
    geo = 0
    optpar = 0
    coords = []
    struct = []
    structure = []
    energies = []
    energy = []
    for line in f:

        if opt == 1 and geo == 1 and not ("---" in line):  # grab the lines of XYZ coordinates
            structure += [line.rstrip()]

        if 'Optimized Parameters' in line:  # Set flags to grab the right strings
            optpar = 1

        if 'Coordinates (Angstroms)' in line:
            if opt == 0:
                opt = 1
                structure = []

        if opt == 1 and "--------------------------" in line:
            if geo == 0:
                geo = 1
            elif geo == 1:
                geo = 0
                opt = 0
        if 'SCF Done' in line:
            energy = filter(None, line.rstrip('\n').split(' '))

        if 'Optimization completed' in line and (opt == 0 and geo == 0):
            energies += [float(energy[4])]
            opt = 0
            geo = 0
            struct += [structure]
            structure = []

        if optpar == 1 and '! ' + scanned in line:
            coord = filter(None, line.rstrip('\n').split(' '))
            coords += [coord[3]]
            optpar = 0

    return struct, energies, coords


def periodictable(elementnumber):
    ptable = {1: 'H', 2: 'He',
              3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
              11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar',
              19: 'K', 20: 'Ca',
              21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
              31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr',
              37: 'Rb', 38: 'Sr',
              39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd',
              49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe',
              55: 'Cs', 56: 'Ba',
              57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy',
              67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu',
              72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
              81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn',
              87: 'Fr', 88: 'Ra',
              89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf',
              99: 'Es', 100: 'Fm', 101: 'Md',
              102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds',
              111: 'Rg', 112: 'Cn',
              113: 'Uut', 114: 'Fl', 115: 'Uup', 116: 'Lv', 117: 'Uus', 118: 'Uuo'}
    element = ptable[elementnumber]
    return element


def genxyzstring(coor, elementnumber):
    x_str = '%10.5f' % coor[0]
    y_str = '%10.5f' % coor[1]
    z_str = '%10.5f' % coor[2]
    element = periodictable(int(elementnumber))
    xyz_string = element + (3 - len(element)) * ' ' + 10 * ' ' + \
        (8 - len(x_str)) * ' ' + x_str + 10 * ' ' + (8 - len(y_str)) * ' ' + y_str + 10 * ' ' + \
        (8 - len(z_str)) * ' ' + z_str + '\n'
    return xyz_string


def getstructures(rawdata, coords, nbo):
    for structure, deg in zip(rawdata, coords):

        g = open(geo_fn + '_' + deg + '_' + '.gjf', 'w')
        cartesian = []
        chk = "%chk=" + geo_fn + '_' + deg + '_' + ".chk" + '\n'
        g.write(chk)
        g.write(card)
        note = geo_fn + '_' + deg + "_" + ", sp at " + met + '/' + bas + '\n\n'
        g.write(note)
        for item in structure:
            coords = filter(None, item.split(' '))
            coordinates = [float(coords[3]), float(coords[4]), float(coords[5])]
            element = coords[1]
            cartesian += [genxyzstring(coordinates, element)]
        g.write('0 1' + '\n')
        for line in cartesian:
            g.write(line)
        nbo = re.sub("_\d\.?\d?\.?\d?\.?", "_" + deg, nbo)
        g.write(nbo)
        g.close()
    return 0


def getbatch(directory, coords):
    """
    Prepare for creating Gaussian batch (bcf) file
    :param directory: path to the destination directory
    :param coords: suffix (deg, ditances) to distiguish files by the coordinate step
    :return: string of input and destination files
    """
    batch = []
    for deg in coords:
        gjf = directory + '\\' + geo_fn + '_' + deg + '_' + '.gjf'
        out = geo_fn + '_' + deg + '_' + '.out'
        batch.append(gjf + ', ' + out)
    return batch


# Round number for pretty output; coordinates A(int), R(1), D(int), F(1)
def list_get(l, coordinate, v = "R"):
    coords_rnd = []
    if coordinate == "A" or coordinate == "D":
        for value in l:
            coords_rnd += [str(int(round(float(value))))]
        return coords_rnd
    else:
        for value in l:
            coords_rnd += [str(round(float(value), 2))]
        return coords_rnd


if __name__ == "__main__":
    # Errors in the input - get the usage and errors
    if main(sys.argv) == 1:
        raise SystemExit

    # Read data in from terminal
    infile = sys.argv[1]

    # Extract file name for later file naming
    geo_fn = infile.split('.')[0]

    # ######## Menu entries ##

    # THEORY
    print '\n'
    print (15 * '-')
    print ("  SET THEORY   ")
    print (15 * '-')
    print (" 1. HF")
    print (" 2. B3LYP")
    print (" 3. M06-2X")
    print (" 4. wB97XD")
    print (" 5. MP2")
    print (" 6. other ..")
    print (15 * '-')
    is_valid = 0
    while not is_valid:
        try:
            metraw = int(raw_input('Enter your choice [1-6] : '))
            is_valid = 1  # set it to 1 to validate input and to terminate the while..not loop
        except ValueError, e:
            print ("'%s' is not a valid entry." % e.args[0].split(": ")[1])
    if metraw == 1:
        met = 'HF'
    elif metraw == 2:
        met = 'B3LYP'
    elif metraw == 3:
        met = 'M06-2X'
    elif metraw == 4:
        met = 'wB97XD'
    elif metraw == 5:
        met = 'MP2'
    elif metraw == 6:
        met = raw_input('\nType the theory level: ')
        if len(met) < 2:
            print "\n ---> Wrong entry. B3LYP will be used."
            met = 'B3LYP'
    else:
        met = 'B3LYP'

    # BASIS SET
    print '\n'
    print (15 * '-')
    print ("  BASIS SET  ")
    print (15 * '-')
    print (" 1. 6-31+G(d)")
    print (" 2. 6-311++G(d,p)")
    print (" 3. Type other ..")
    print (15 * '-')

    is_valid = 0
    while not is_valid:
        try:
            basraw = int(raw_input('Enter your choice [1-3] : '))
            is_valid = 1  # set it to 1 to validate input and to terminate the while..not loop
        except ValueError, e:
            print ("'%s' is not a valid entry." % e.args[0].split(": ")[1])
    if basraw == 1:
        bas = '6-31+G(d)'
    elif basraw == 2:
        bas = '6-311++G(d,p)'
    elif basraw == 3:
        bas = raw_input('\nType the basis set: ')
        if len(bas) < 2:
            print "\n ---> Wrong entry. 6-311++G(d,p) will be used."
            bas = '6-311++G(d,p)'
    else:
        bas = '6-311++G(d,p)'

    # How to run NBO
    print '\n'
    print (15 * '-')
    print ("  NBO OPTIONS  ")
    print (15 * '-')
    print (" 1. create .47 file only (archive)")
    print (" 2. run linked G09-NBO")
    print (" 3. run compiled G09-NBO binaries")
    print (15 * '-')

    is_valid = 0
    while not is_valid:
        try:
            nboraw = int(raw_input('Enter your choice [1-3] : '))
            is_valid = 1  # set it to 1 to validate input and to terminate the while..not loop
        except ValueError, e:
            print ("'%s' is not a valid entry." % e.args[0].split(": ")[1])

    if nboraw == 2:
        option = 'run linked G09-NBO'
    if nboraw == 3:
        option = 'run compiled G09-NBO binaries'

    # NBO keywords
    if nboraw == 1:
        deg = ''
        nbo = '\n' + "$NBO archive FILE=" + geo_fn + '_0' + "_ $END" + '\n\n'
        option = 'create .47 file only (archive)'

    # NBO keywords for option 2,3
    else:
        print '\n'
        print '\n' + "    You will need to choose NBO keywords."
        print "    Use GennboHelper to copy/paste keywords to the input 3."
        print (15 * '-')
        print ("  NBO KEYWORDS  ")
        print (15 * '-')
        print ("    1. Keyword set 1 (NBOSUM DIST BNDIDX DIPOLE=0.02 E2PERT=5 PRINT=2)")
        print ("    2. Keyword set 2 (STERIC=0.5 DIST E2PERT=5 PRINT=2)")
        print ("    3. other ..")
        print (15 * '-')

        is_valid = 0
        while not is_valid:
            try:
                nbokey = int(raw_input('Enter NBO keywords : '))
                is_valid = 1  # set it to 1 to validate input and to terminate the while..not loop
            except ValueError, e:
                print ("'%s' is not a valid entry." % e.args[0].split(": ")[1])

        if nbokey == 1:
            keywords = 'NBOSUM DIST BNDIDX DIPOLE=0.02 E2PERT=5 PRINT=2'
            nbo = '\n' + "$NBO " + keywords + " $END" + '\n\n'
        elif nbokey == 2:
            keywords = 'STERIC=0.5 DIST E2PERT=5 PRINT=2'
            nbo = '\n' + "$NBO " + keywords + " $END" + '\n\n'
        elif nbokey == 3:
            keywords = raw_input('\nType/paste the keywords (space separated): ')
            nbo = '\n' + "$NBO " + keywords + " $END" + '\n\n'
            if len(keywords) < 3:
                print "\n ---> Wrong entry. 'DEFAULT' will be used."
                keywords = 'NBOSUM NRT STERIC=0.5 DIST BNDIDX DIPOLE=0.02 E2PERT=5 PRINT=2'
        else:
            keywords = 'NBOSUM NRT STERIC=0.5 DIST BNDIDX DIPOLE=0.02 E2PERT=5 PRINT=2'
            nbo = '\n' + "$NBO archive FILE=" + geo_fn + '_0' + "_ $END" + '\n\n'
            option = 'create .47 file only (archive)'

            # ####### Menu ends ########
    print "\n\n  Theory/Basis:   " + met + '/' + bas + '\n'
    print "  NBO options:   " + option + '\n'

    scanned = getscan(infile)[1]
    sccoord = getscan(infile)[2]
    print "  Scanned coordinate is: " + sccoord + " with label: " + scanned + '\n'

    # Route card
    if nboraw == 1:
        card = "# " + met + '/' + bas + " pop=nboread sp nosymm" + '\n\n'
    elif nboraw == 2:
        card = "# " + met + '/' + bas + " external=C:\G09W\gaunbo6.bat POP=NBO6Read sp nosymm" + '\n\n'
    elif nboraw == 3:
        card = "# " + met + '/' + bas + " pop=nbo6read sp nosymm" + '\n\n'
    else:
        card = "# " + met + '/' + bas + " pop=nboread sp nosymm" + '\n\n'

    rawdata, energies, coords = getrawdata(infile)

    # Format coords steps
    regexR = re.compile("^R.+")
    regexA = re.compile("^A.+")
    regexD = re.compile("^D.+")

    if regexA.match(scanned):
        coordinate = "A"
    elif regexR.match(scanned):
        coordinate = "R"
    elif regexD.match(scanned):
        coordinate = "D"
    else:
        coordinate = "F"

    # Call rounding function
    coords_round = list_get(coords, coordinate)

    # Print string list
    # for value in coords_round:
        # print value

    structures = getstructures(rawdata, coords_round, nbo)

    # Write results to a file
    g = open(geo_fn + '_energies.dat', 'w')  # get energies for graph
    for n in range(0, len(coords_round)):
        g.write(coords_round[n] + '\t' + str(energies[n]) + '\n')
    g.close()
    decor = len(os.path.dirname(os.path.realpath(__file__))) + 31
    print_frame_top(decor, '*')
    directory = os.path.dirname(os.path.realpath(__file__))
    print str(
        len(energies)) + " files and Gaussian " + geo_fn + "_batch.bcf batch file are in directory: " + '\n' + directory
    print_frame_bot(decor, '*')

    # write gaussian batch file .bcf
    batchf = getbatch(directory, coords_round)
    b = open(geo_fn + '_batch.bcf', 'w')  # get files into batch file
    b.write("!" + '\n'
                  "!User created batch file" + '\n'
                                               "!start=1" + '\n'
                                                            "!" + '\n')
    for n in range(0, len(batchf)):
        b.write(str(batchf[n]) + '\n')
    b.close()

    # Prepare for plots and prints
    maxmindif = (max(energies) - min(energies)) * 627.51
    coords = map(float, coords)  # list of strings to floats
    rangeX = abs(max(coords) - min(coords))  # find the range
    firstEne = min(energies)  # reference energy

    # reformat energy list
    ene = list(rundif(energies, firstEne))  # subtract current energy from reference*627.51
    ene = ["%.2f" % member for member in ene]  # format numbers
    ene = map(float, ene)

    # ploting
    if coordinate == "A":
        plotcoord = "Angle, deg"
    elif coordinate == "R":
        plotcoord = "Distance, Angstrom"
    elif coordinate == "D":
        plotcoord = "Dihedral, deg"
    else:
        plotcoord = ""

    try:
        import pylab as pl
        from pylab import *

        pylab_available = True
    except ImportError:
        pylab_available = False
        print "Pylab and matplotlib modules were not imported. Use the .dat file to print"

    if pylab_available:

        X = coords
        y = ene
        pl.ylim(min(y) - 0.1 * maxmindif, max(y) + 0.1 * maxmindif)
        pl.xlim(min(X) - 0.1 * rangeX, max(X) + 0.1 * rangeX)
        pl.xlabel('Coordinate (' + plotcoord + ')')
        pl.ylabel('rel Energy (kcal/mol)')
        pl.plot(X, y, 'bo', label='(max-min)dE=%5.2f' % maxmindif + ' kcal/mol')
        pl.plot(X, y, ':k')
        plt.axhline(y=0, xmin=0, xmax=1, linewidth=1, color='b')
        pl.legend(loc='upper right')
        locs, labels = yticks()
        yticks(locs, map(lambda x: "%.1f" % x, locs * 1e0))
        if coordinate == "A" or coordinate == "D":
            pl.xticks(np.arange(min(X), max(X)+1, max(X)/(len(coords)-1)))
        text(0.0, 1.01, '', fontsize=10, transform=gca().transAxes)
        pl.show()
        print
    else:
        exit(1)
