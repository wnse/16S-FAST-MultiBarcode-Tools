import argparse
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Barcode Analysis')
# parser.add_argument('-FB', '--Fbarcode', help='F barcode list, eg. F1,F2,F3...', required=True, nargs="+")
parser.add_argument('-RB', '--Rbarcode', help='R barcode list, eg. R1,R2,R3...', required=True, nargs="+")
parser.add_argument('-N', '--Name', help='Sample names for each barcode', required=True, nargs="+")
parser.add_argument('-d', '--debug', help='for debug', action='store_true')
args = parser.parse_args()

for i in  list(zip(args.Name, args.Rbarcode)):
    print(i[0],i[1])

print(args.debug)
