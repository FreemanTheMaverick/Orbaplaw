import re
import argparse
from Orbaplaw import WaveFunction as wfn
from Orbaplaw import Population as pop
from Orbaplaw import Localization as loc
from Orbaplaw import OrbitalAlignment as oa
from Orbaplaw import NaturalBondOrbitalMethods as nbo

def parse_range(value):
	"""Parse a string like '3-10,12,17-20' into a list of integers."""
	ranges = []
	parts = value.split(',')
	
	for part in parts:
		if '-' in part:
			# Handle ranges like "3-10"
			subparts = part.split('-')
			if len(subparts) != 2:
				raise argparse.ArgumentTypeError(f"Invalid range: '{part}'")
			try:
				start = int(subparts[0])
				end = int(subparts[1])
			except ValueError:
				raise argparse.ArgumentTypeError(f"Non-integer in range: '{part}'")
			if start > end:
				raise argparse.ArgumentTypeError(f"Invalid range order: '{part}'")
			ranges.extend(range(start, end + 1))  # Inclusive of end
		else:
			# Handle single numbers like "12"
			try:
				num = int(part)
				ranges.append(num)
			except ValueError:
				raise argparse.ArgumentTypeError(f"Invalid number: '{part}'")
	return ranges

def parse_alpha_range(value):
	if bool(re.search(r'[a-zA-Z]', value)):
		return value
	else:
		return parse_range(value)

def main():
	parser = argparse.ArgumentParser(description = "Command-line tool for Orbaplaw")
	subparsers = parser.add_subparsers(
			title = "Job types",
			dest = "job_type",
			required = True
	)

	# Parser for population analysis
	parser_pop = subparsers.add_parser(
			"pop",
			help = "Population Analysis"
	)
	parser_pop.add_argument(
			"-i", "--input",
			required = True,
			help = "Mwfn file (Required)"
	)
	parser_pop.add_argument(
			"--charge",
			default = "Lowdin",
			help = "Charge type (Default: %(default)s)"
	)
	parser_pop.add_argument(
			"--space",
			default = "occ",
			type = parse_alpha_range,
			help = "Orbital space to analyze ('occ', 'vir', '0-2,4,6,8-10') (Default: %(default)s)"
	)

	# Parser for localization
	parser_loc = subparsers.add_parser(
			"loc",
			help = "Orbital localization"
	)
	parser_loc.add_argument(
			"-i", "--input",
			required = True,
			help = "Original mwfn file (Required)"
	)
	parser_loc.add_argument(
			"-o", "--output",
			required = True,
			help = "Mwfn file for exported localized orbitals (Required)"
	)
	parser_loc.add_argument(
			"--method",
			default = "PipekMezey-Lowdin",
			help = "Localization method ('PipekMezey[-Lowdin/-Mulliken]', 'Foster-Boys', 'Orbitalet[ 0.7959/ gamma_e]') (Default: %(default)s)"
	)
	parser_loc.add_argument(
			"--space",
			default = "occ",
			type = parse_alpha_range,
			help = "Orbital space to localize ('occ', 'vir', '0-2,4,6,8-10') (Default: %(default)s)"
	)

	# Parser for fragment alignment
	parser_famo = subparsers.add_parser(
			"famo",
			help = "Fragment aligned molecular orbital"
	)
	parser_famo.add_argument(
			"-i", "--input",
			required = True,
			help = "Original mwfn file of the whole molecule (Required)"
	)
	parser_famo.add_argument(
			"-f", "--fragments",
			required = True,
			nargs = '+',
			help = "Fragment mwfn file(s)"
	)
	parser_famo.add_argument(
			"-o", "--output",
			required = True,
			help = "Mwfn file for exported FAMOs (Required)"
	)
	parser_famo.add_argument(
			"--diagmat",
			default = False,
			type = bool,
			help = "Whether to diagonalize the matched space (Default: %(default)s)"
	)
	parser_famo.add_argument(
			"--diagmis",
			default = True,
			type = bool,
			help = "Whether to diagonalize the mismatched space (Default: %(default)s)"
	)

	# Parser for spin alignment
	parser_sno = subparsers.add_parser(
			"sno",
			help = "Spin natural orbital"
	)
	parser_sno.add_argument(
			"-i", "--input",
			required = True,
			help = "Original mwfn file (Required)"
	)
	parser_sno.add_argument(
			"-o", "--output",
			required = True,
			help = "Mwfn file for exported SNOs (Required)"
	)
	parser_sno.add_argument(
			"--diagmat",
			default = False,
			type = bool,
			help = "Whether to diagonalize the matched space (Default: %(default)s)"
	)
	parser_sno.add_argument(
			"--diagmis",
			default = True,
			type = bool,
			help = "Whether to diagonalize the mismatched space (Default: %(default)s)"
	)

	# Parser for natural atomic orbitals
	parser_nao = subparsers.add_parser(
			"nao",
			help = "Natural atomic orbital"
	)
	parser_nao.add_argument(
			"-i", "--input",
			required = True,
			help = "Original mwfn file (Required)"
	)
	parser_nao.add_argument(
			"-o", "--output",
			required = True,
			help = "Mwfn file for exported NAOs (Required)"
	)

	# Parser for principal interacting orbitals
	parser_pio = subparsers.add_parser(
			"pio",
			help = "Principal interacting orbital"
	)
	parser_pio.add_argument(
			"-i", "--input",
			required = True,
			help = "Original mwfn file (Required)"
	)
	parser_pio.add_argument(
			"--fragments",
			required = True,
			nargs = 2,
			help = "Fragmentation scheme (e.g., '0-2,4 6,8-10') (Required)"
	)
	parser_pio.add_argument(
			"-o", "--output",
			required = True,
			nargs = 2,
			help = "Mwfn file for exported PIOs and PIMOs, respectively (Required)"
	)

	# Parser for natural bond orbitals
	parser_nbo = subparsers.add_parser(
			"nbo",
			help = "Natural (fragment) bond orbital"
	)
	parser_nbo.add_argument(
			"-i", "--input",
			required = True,
			help = "Original mwfn file (Required)"
	)
	parser_nbo.add_argument(
			"-o", "--output",
			required = True,
			nargs = 2,
			help = "Mwfn file for exported NHOs and NBOs, respectively (Required)"
	)
	parser_nbo.add_argument(
			"--fragments",
			nargs = '*',
			default = [],
			help = "Fragmentation scheme (e.g., '0-2,4 6,8-10 ...') (Default: _EMPTY_ -> each atom is a fragment)"
	)
	parser_nbo.add_argument(
			"--maxnfrags",
			default = -1,
			type = int,
			help = "Maximum number of fragments involved in a bonding scheme (Default: %(default)s -> no limit)"
	)
	parser_nbo.add_argument(
			"--maxnnbos",
			default = -1,
			type = int,
			help = "Maximum number of bonding schemes to find (Default: %(default)s -> no limit)"
	)
	parser_nbo.add_argument(
			"--occ_thres",
			default = 0.95,
			type = float,
			help = "Occupation threshold for recognizing pNFBOs (ranging from 0 to 1) (Default: %(default)s)"
	)
	parser_nbo.add_argument(
			"--multi_thres",
			default = 1.,
			type = float,
			help = "Occupation threshold for recognizing multi-electron bonds (ranging from 0 to 1) (Default: %(default)s)"
	)
	parser_nbo.add_argument(
			"--pdeg_thres",
			default = 1e-5,
			type = float,
			help = "Occupation difference threshold for recognizing degenerate pNFBOs (Default: %(default)s)"
	)
	parser_nbo.add_argument(
			"--deg_thres",
			default = 0.,
			type = float,
			help = "Occupation difference threshold for recognizing degenerate NFBOs (Default: %(default)s)"
	)


	args = parser.parse_args()


	# Running population analysis
	if args.job_type == "pop":
		mwfn = wfn.MultiWaveFunction(args.input)
		pop.PopulationAnalyzer(mwfn, args.charge, args.space)

	# Running localization
	if args.job_type == "loc":
		orig_mwfn = wfn.MultiWaveFunction(args.input)
		loc_mwfn = loc.Localizer(orig_mwfn, args.method, args.space)
		loc_mwfn.Export(args.output)

	# Running fragment alignment
	if args.job_type == "famo":
		whole_mwfn = wfn.MultiWaveFunction(args.input)
		frag_mwfns = []
		for frag_mwfn in args.fragments:
			frag_mwfns.append(mwfn.MultiWaveFunction(frag_mwfn))
		famo_mwfn = oa.FragmentAlignment(whole_mwfn, frag_mwfns, args.diagmat, args.diagmis)
		famo_mwfn.Export(args.output)

	# Running spin alignment
	if args.job_type == "sno":
		orig_mwfn = wfn.MultiWaveFunction(args.input)
		sno_mwfn = oa.SpinAlignment(orig_mwfn, args.diagmat, args.diagmis)
		sno_mwfn.Export(args.output)

	# Running natural atomic orbitals
	if args.job_type == "nao":
		orig_mwfn = wfn.MultiWaveFunction(args.input)
		nao_mwfn = nbo.NaturalAtomicOrbital(orig_mwfn)
		nao_mwfn.Export(args.output)

	# Running principal interacting orbitals
	if args.job_type == "pio":
		orig_mwfn = wfn.MultiWaveFunction(args.input)
		nao_mwfn = nbo.NaturalAtomicOrbital(orig_mwfn)
		frags = []
		for frag in args.fragments:
			frags.append(parse_range(frag))
		pio_mwfn, pimo_mwfn = nbo.PrincipalInteractingOrbital(nao_mwfn, frags)
		pio_mwfn.Export(args.output[0])
		pimo_mwfn.Export(args.output[1])

	# Running natural bond orbitals
	if args.job_type == "nbo":
		orig_mwfn = wfn.MultiWaveFunction(args.input)
		nao_mwfn = nbo.NaturalAtomicOrbital(orig_mwfn)
		frags = []
		for frag in args.fragments:
			frags.append(parse_range(frag))
		nho_mwfn, nbo_mwfn = nbo.NaturalBondOrbital(nao_mwfn, frags, args.maxnfrags, args.maxnnbos, args.occ_thres, args.multi_thres, args.pdeg_thres, args.deg_thres)
		nho_mwfn.Export(args.output[0])
		nbo_mwfn.Export(args.output[1])
