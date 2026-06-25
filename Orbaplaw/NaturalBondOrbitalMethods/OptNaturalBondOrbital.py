import numpy as np
import scipy.linalg as sl
import pandas as pd
import itertools as it
import copy as cp
import Maniverse as mv
from Orbaplaw import Integrals as eint
from Orbaplaw import Localization as loc


def RecursiveBridge(pNBO, this_comb, this_i, is_meb, code):
	this_pos = list(pNBO.columns.get_level_values(0).unique()).index(this_comb)
	for comb, subdf in list(pNBO.groupby(level=0, axis=1))[this_pos:]:
		if comb is this_comb or len(set(comb) & set(this_comb)) == 0:
			continue
		for i in subdf.columns:
			if is_meb[i] > 0:
				continue
			for frag in comb:
				if frag in this_comb:
					if abs(np.dot(pNBO[(this_comb, this_i), (frag, slice(None))], pNBO[(comb, i), (frag, slice(None))])) > multi_thres:
						is_meb[i] = code
						RecursiveBridge(pNBO, comb, i, is_meb, code)


class Obj(mv.Objective):
	def __init__(self, P, basis_list):
		super().__init__()
		self.P = P
		self.this_P = np.zeros_like(P)
		self.basis_list = basis_list

	def Calculate(self, X, derivatives):
		nmats = len(X) // 2
		t = X[:nmats]
		C = X[nmats:]
		P = self.P
		this_P = self.this_P
		basis_list = self.basis_list
		if 0 in derivatives:
			this_P.fill(0)
			for imat in range(nmats):
				ti = t[imat][:, 0]
				ni = ( np.tanh(ti) + 1 ) / 2
				print(ni)
				Ci = C[imat]
				bases = self.basis_list[imat]
				this_P[np.ix_(bases, bases)] += Ci * ni @ Ci.T
			self.Value = 0.5 * np.linalg.norm(this_P - P) ** 2
		if 1 in derivatives:
			Gradient = [None] * len(basis_list) * 2
			E = this_P - P
			for imat in range(nmats):
				ti = t[imat][:, 0]
				ni = ( np.tanh(ti) + 1 ) / 2
				Ci = C[imat]
				bases = self.basis_list[imat]
				Gradient[imat] = np.diag(Ci.T @ E[np.ix_(bases, bases)] @ Ci) / np.cosh(ti) ** 2 / 2
				Gradient[imat + nmats] = 2 * E[np.ix_(bases, bases)] @ Ci @ np.diag(ni)
			self.Gradient = Gradient


def group_close_indices(lst, threshold):
	"""
	Group indices of a sorted list where consecutive element differences
	are <= threshold.

	Args:
		lst: sorted list of numbers (or any comparable type)
		threshold: non-negative number; maximum allowed gap to stay in same group

	Returns:
		List of lists, each containing indices that form a group.
	"""
	if not lst:
		return []
	
	groups = []
	current_group = [0] # start with index of first element
	
	for i in range(1, len(lst)):
		if abs(lst[i] - lst[i-1]) <= threshold:
			# close enough -> stay in same group
			current_group.append(i)
		else:
			# gap too large -> finish current group and start a new one
			groups.append(current_group)
			current_group = [i]
	
	# append the last group
	groups.append(current_group)
	return groups


def generateOptNaturalBondOrbital(basis_indices_by_frag, adjacency, P_, T, S, maxnfrags, maxnnbos, occ_thres, multi_thres, pdeg_thres, deg_thres):
	# P - The NAO-based density matrix.
	# T - The NAO coefficient matrix in AO basis set.
	# S - The AO-based overlap matrix.

	nbasis = P_.shape[0]
	nfrags_tot = len(basis_indices_by_frag)
	idx = pd.IndexSlice

	# Storing the density matrix as a pd.DataFrame
	frags_basis = []
	for ibasis in range(nbasis):
		for jfrag in range(nfrags_tot):
			if ibasis in basis_indices_by_frag[jfrag]:
				frags_basis.append((jfrag, ibasis))
				continue
	P = pd.DataFrame(
			P_,
			index = pd.MultiIndex.from_tuples(frags_basis),
			columns = pd.MultiIndex.from_tuples(frags_basis)
	)
	Ptot = P.copy()

	# Searching for pNBOs and pNHOs.
	pNHO = pd.DataFrame(
			index = pd.MultiIndex.from_tuples([("occ", -1, -1)] + [("coeff", frag, basis) for (frag, basis) in P.index.tolist()]),
			columns = pd.MultiIndex.from_tuples([], names=["comb", "pnbo", "frag"])
	)
	npnbos_tot = 0
	for nfrags in range(1, maxnfrags + 1):
		pNBO = pd.DataFrame(
				index = pd.MultiIndex.from_tuples([("occ", -1, -1)] + [("coeff", frag, basis) for (frag, basis) in P.index.tolist()]),
				columns = pd.MultiIndex.from_tuples([], names=["comb", "pnbo"])
		)
		combs = list(it.combinations(range(nfrags_tot), nfrags)) # Combinations of fragments.
		npnbos = 0
		for comb in combs:
			if len(comb)==2:
				continue
			bonded = True
			for ifrag in range(len(comb)):
				fragi = list(comb)[ifrag]
				for jfrag in range(ifrag):
					fragj = list(comb)[jfrag]
					if not adjacency[fragi][fragj]:
						bonded = False
			if not bonded:
				continue
			Pblock = np.asarray(P.loc[idx[list(comb), :], idx[list(comb), :]]) # The density matrix block of this combination
			Nblock, Hblock = np.linalg.eigh(Pblock) # Occupation and orbital vectors
			to_delete = np.where(Nblock < occ_thres)[0] # Ignoring all small occupancies.
			Nblock = np.delete(Nblock, to_delete)
			if len(Nblock) == 0:
				continue
			Hblock = np.delete(Hblock, to_delete, axis=1)
			frag_bases = [index for index in pNBO.index if index[1] in comb] # Basis indices of this combination [("coeff", frag, basis)]
			print("fuck", [(frozenset(comb), i) for i in range(npnbos, npnbos + len(Nblock))])
			pNBO = pNBO.join( pd.DataFrame(
					np.vstack((Nblock, Hblock)),
					index = pd.MultiIndex.from_tuples([("occ", -1, -1)] + frag_bases),
					columns = pd.MultiIndex.from_tuples(
							[(frozenset(comb), i) for i in range(npnbos, npnbos + len(Nblock))]
					)
			))
			npnbos += len(Nblock)

		# Checking and excluding bridge bonds between combinations
		bridge = [0] * len(pNBO.columns)
		for combi, ipnbo in pNBO.columns:
			for combj, jpnbo in pNBO.columns:
				if combi == combj or ipnbo >= jpnbo:
					continue
				for frag in set(combi) & set(combj):
					vi = pNBO.loc[idx["coeff", frag, :], idx[:, ipnbo]].to_numpy()[:, 0]
					vj = pNBO.loc[idx["coeff", frag, :], idx[:, jpnbo]].to_numpy()[:, 0]
					if abs(np.dot(vi, vj)) > multi_thres ** 2 * np.linalg.norm(vi) * np.linalg.norm(vj):
						bridge[ipnbo] = bridge[jpnbo] = 1
		for combi, ipnbo in pNBO.columns:
			if bridge[ipnbo] == 1:
				pNBO.drop((frozenset(combi), ipnbo), axis=1)
		pNBO.columns = pNBO.columns.set_levels(range(len(pNBO.columns)), level=1)
		if len(pNBO.columns) == 0:
			continue

		while True:
			t = []
			C = []
			basis_list = []
			for comb, subdf in pNBO.groupby(level=0, axis=1):
				count = 0
				this_t = []
				this_C = []
				basis_list.append([index[1] for index in P.columns if index[0] in comb]) # Basis indices of this combination [(frag, basis)]
				for i in subdf.columns:
					count += 1
					n = min(pNBO.at[("occ", -1, -1), i], 0.99999)
					this_t.append([ np.arctanh(2 * n - 1) ])
					this_C.append(pNBO.loc[idx["coeff", list(comb), :], idx[:, i[1]]].T.values[0])
				this_t = np.array(this_t)
				this_C = np.array(this_C).T
				t.append(mv.Euclidean(this_t))
				C.append(mv.Stiefel(this_C))
			obj = Obj(P.to_numpy(), basis_list)
			M = mv.Iterate(obj, t + C)
			converged = mv.LBFGS(
					M, (1e-5, 1e-4, 1e-3),
					100, 10000, 0.1, 0.75, 10, 1
			)
			if not converged:
				raise RuntimeError("Not converged!")
			kpnbo = 0;
			nmats = len(t)
			print(pNBO.columns)
			for imat in range(nmats):
				for j in range(M.Ms[imat].P.shape[0]):
					print(kpnbo)
					pNBO.loc[idx["occ", -1, -1], idx[:, kpnbo]] = ( np.tanh(M.Ms[imat].P[j, 0]) + 1 ) / 2
					pNBO.loc[idx["coeff", :, basis_list[imat]], idx[:, kpnbo]] = M.Ms[imat + nmats].P[:, j]
					kpnbo += 1
			for combi, ipnbo in pNBO.columns:
				if pNBO.at[idx["occ", -1, -1], idx[combi, ipnbo]] < occ_thres:
					pNBO.drop((combi, ipnbo), axis=1, inplace=True)
			if len(pNBO.columns) == 0:
				break
			else:
				level0_old = pNBO.columns.get_level_values(0)
				level1_new = range(npnbos_tot, npnbos_tot + len(pNBO.columns)) if kpnbo == len(pNBO.columns) else range(len(pNBO.columns))
				pNBO.columns = pd.MultiIndex.from_arrays([level0_old, level1_new])
			if kpnbo == len(pNBO.columns):
				npnbos_tot += len(pNBO.columns)
				break

		if len(pNBO.columns) == 0:
			continue

		# Identifying and partially localizing the degenerate orbitals
		skip_list = []
		for comb, ipnbo in pNBO.columns:
			deg_list = [ipnbo]
			for combj, jpnbo in pNBO.columns:
				if comb != combj or ipnbo >= jpnbo:
					continue
				if jpnbo in skip_list:
					continue
				if abs(
						pNBO.at[idx["occ", -1, -1], idx[comb, ipnbo]]
						- pNBO.at[idx["occ", -1, -1], idx[comb, jpnbo]]
				) < pdeg_thres:
					deg_list.append(jpnbo)
			if len(deg_list) == 1:
				continue
			skip_list += deg_list
			pNBO.fillna(0, inplace=True)
			if len(comb) == 1:
				continue
			H = pNBO.loc[idx["coeff", :, :], idx[:, deg_list]].to_numpy()
			U = loc.oldPipekMezey(T @ H, S, basis_indices_by_frag, "Lowdin", None) # T - NAO in AO basis; T@H - degenerate pNBO in AO basis
			pNBO.loc[idx["coeff", :, :], idx[:, deg_list]] = H @ U

		# Subtracting the density of the pNBOs
		for comb, ipnbo in pNBO.columns:
			occ = pNBO.at[idx["occ", -1, -1], idx[comb, ipnbo]]
			vec = pNBO.loc[idx["coeff", list(comb), :], idx[:, ipnbo]].fillna(0).to_numpy()[:, 0]
			P.loc[idx[list(comb), :], idx[list(comb), :]] -= occ * np.outer(vec, vec)

		# Decomposing pNBOs into pNHOs in average over bridge
		pNHO_dict = dict()
		for combi, ipnbo in pNBO.columns:
			for frag in combi:
				vec = np.zeros(P.shape[0])
				vec[basis_indices_by_frag[frag]] = pNBO.loc[idx["coeff", frag, :], idx[:, ipnbo]].to_numpy()[:, 0]
				vec = np.insert(vec, 0, [0])
				pNHO_dict[(combi, ipnbo, frag)] = vec
		this_pNHO = pd.DataFrame(pNHO_dict, index = pNHO.index)
		this_pNHO.columns = pd.MultiIndex.from_tuples(this_pNHO.columns)
		replicate = []
		for combi, ipnbo, fragi in this_pNHO.columns:
			veci = this_pNHO.loc[idx["coeff", fragi, :], idx[:, ipnbo, fragi]].to_numpy()[:, 0]
			if np.linalg.norm(veci) < 1e-12:
				continue
			for combj, jpnbo, fragj in this_pNHO.columns:
				if combi != combj or ipnbo >= jpnbo or fragi != fragj:
					continue
				if jpnbo in replicate:
					continue
				vecj = this_pNHO.loc[idx["coeff", fragj, :], idx[:, jpnbo, fragj]].to_numpy()[:, 0]
				if np.linalg.norm(vecj) < 1e-12:
					continue
				overlap = np.dot(veci, vecj)
				if abs(overlap) > multi_thres * np.linalg.norm(veci) * np.linalg.norm(vecj):
					this_pNHO.loc[idx["coeff", fragi, :], idx[:, ipnbo, fragi]] += vecj.reshape([-1, 1]) * np.sign(overlap)
					replicate.append((combj, jpnbo, fragj))
		for combi, ipnbo, fragi in this_pNHO.columns:
			if (combi, ipnbo, fragi) in replicate:
				this_pNHO.drop((combi, ipnbo, fragi), axis=1, inplace=True)
			else:
				this_pNHO.loc[idx["coeff", fragi, :], idx[:, ipnbo, fragi]] /= np.linalg.norm(this_pNHO.loc[idx["coeff", fragi, :], idx[:, ipnbo, fragi]])

		# Merging into the pNHO list
		pNHO = pNHO.join(this_pNHO)

	# Orthogonalizing the pNHOs to NHOs
	for frag in range(nfrags):
		I = pNHO.loc[idx["coeff", frag, :], idx[:, :, frag]].to_numpy()
		J = np.zeros([I.shape[1], I.shape[1]])
		for i in range(I.shape[1]):
			J[i, i] = np.linalg.norm(I[:, i])
		Sblock = I.T @ I
		eigval, eigvec = np.linalg.eigh(J @ Sblock @ J)
		O = J @ eigvec @ np.diag( np.sqrt( 1 / np.abs(eigval) ) ) @  eigvec.T
		I = I @ O
		pNHO.loc[idx["coeff", frag, :], idx[:, :, frag]] = I
		Pblock = Ptot.loc[idx[frag, :], idx[frag, :]].to_numpy()
		J = I.T @ Pblock @ I
		pNHO.loc[idx["occ", -1, -1], idx[:, :, frag]] = np.diag(J)
	pNHO.fillna(0, inplace=True)

	# NHO-based density matrix and occupation numbers of NHOs
	NHOmat = pNHO.loc[idx["coeff", :, :]].to_numpy()
	P_NHO = NHOmat.T @ Ptot.to_numpy() @ NHOmat
	for inho, (combi, ipnbo, fragi) in enumerate(pNHO.columns):
		pNHO.loc[idx["occ", -1, -1], idx[:, ipnbo, fragi]] = np.diag(P_NHO)[inho]
	P_NHO = pd.DataFrame(P_NHO, index=pNHO.columns, columns=pNHO.columns)

	# Diagonalization of NHO-based density matrix
	NBO = pd.DataFrame(index = [("occ", -1)] + [("coeff", i) for i in range(pNHO.shape[1])])
	current_pnbo = -1
	jnbo = 0
	for combi, ipnbo, fragi in pNHO.columns:
		if ipnbo == current_pnbo:
			continue
		current_pnbo = ipnbo
		occs, vecs = np.linalg.eigh(P_NHO.loc[idx[:, ipnbo, :], idx[:, ipnbo, :]])
		occs = occs[::-1]
		vecs = vecs[:, ::-1]
		deg_lists = group_close_indices(list(occs), deg_thres)
		for deg_list in deg_lists:
			if len(deg_list) == 1:
				continue
			I = pNHO.loc[idx["coeff", :, :], idx[:, ipnbo, :]].to_numpy()
			H = vecs[:, deg_list]
			U = loc.oldPipekMezey(T @ I @ H, S, basis_indices_by_frag, "Lowdin", None) # T - NAO in AO basis
			vecs[:, deg_list] = H @ U
		nho_positions = P_NHO.index.get_locs(idx[:, ipnbo, :]).tolist()
		for k in range(len(occs)):
			NBO = NBO.join(pd.DataFrame(
				[occs[k]] + vecs[:, k].tolist(),
				index = [("occ", -1)] + [("coeff", pos) for pos in nho_positions],
				columns=[(combi, ipnbo, jnbo)]
			))
			jnbo += 1
	NBO.fillna(0, inplace=True)
	NBO.index = pd.MultiIndex.from_tuples(NBO.index)
	NBO.columns = pd.MultiIndex.from_tuples(NBO.columns)

	return pNHO, NBO

def OptNaturalBondOrbital(
		nao_mwfn,
		nao_info,
		frags = [],
		adjacency = None,
		maxnfrags = -1,
		maxnnbos = -1,
		occ_thres = 0.95,
		multi_thres = 1,
		pdeg_thres = 1e-5,
		deg_thres = 1e-5): # By default, every atom is a fragment, which is the case of NBO. By combining atoms into fragments one extends NBO to natural fragment bond orbital (NFBO).
	idx = pd.IndexSlice
	if maxnfrags == -1:
		if frags == []:
			maxnfrags = nao_mwfn.getNumCenters()
		else:
			maxnfrags = len(frags)
	frags=[[i] for i in range(nao_mwfn.getNumCenters())] if frags == [] else frags
	if adjacency is None: # All connected
		adjacency = np.zeros([len(frags), len(frags)]) + 1
	basis_indices_by_center = nao_mwfn.Atom2BasisList()
	basis_indices_by_frag = []
	for frag in frags:
		basis_indices_this_fragment = []
		for icenter in frag:
			if icenter >= nao_mwfn.getNumCenters():
				raise RuntimeError("Atom index out of range!")
			basis_indices_this_fragment.extend(basis_indices_by_center[icenter])
		basis_indices_by_frag.append(basis_indices_this_fragment)
	nbasis = nao_mwfn.getNumBasis()
	if nao_mwfn.Overlap.shape != tuple([nbasis] * 2):
		nao_mwfn.Overlap = eint.PyscfOverlap(nao_mwfn, nao_mwfn)
	S = nao_mwfn.Overlap
	nho_mwfn = nao_mwfn.Clone()
	nbo_mwfn = nao_mwfn.Clone()
	print("Natural (fragment) bond orbitals:")
	for spin in nao_mwfn.getSpins():
		C = nao_mwfn.getCoefficientMatrix(spin)
		if maxnnbos == -1:
			maxnnbos = round(nao_mwfn.getNumElec(spin))
			if spin == 0:
				maxnnbos //= 2
		P = None
		match spin:
			case 0:
				P = nao_info["NAO_density_matrix"]
			case 1:
				P = nao_info["NAO_alpha_density_matrix"]
			case 2:
				P = nao_info["NAO_beta_density_matrix"]
		if spin == 0:
			occ_thres *= 2
			deg_thres *= 2
		nhos, nbos = generateOptNaturalBondOrbital(
				basis_indices_by_frag,
				adjacency,
				P, C, S,
				maxnfrags, maxnnbos,
				occ_thres, multi_thres,
				pdeg_thres, deg_thres
		)
		print("Spin", spin) # Printing NBO and NHO information
		this_comb = None
		this_ipnbo = -1
		for combi, ipnbo, inbo in nbos.columns:
			if combi != this_comb:
				print("\nFragment combination", combi)
				this_comb = combi
			if ipnbo != this_ipnbo:
				print("Interaction", ipnbo)
				this_ipnbo = ipnbo
			occ = nbos.at[idx["occ", -1], idx[combi, ipnbo, inbo]]
			print(f"NBO_{inbo + (nbasis if spin == 2 else 0)} ({occ:.3f}) =", end='')
			for index in nbos.index:
				if index[0] == "occ":
					continue
				jnho = index[1]
				coeff = nbos.at[idx["coeff", jnho], idx[combi, ipnbo, inbo]]
				if coeff == 0:
					continue
				occ = nhos.loc[idx["occ", -1, -1]].iloc[jnho]
				frag = nhos.columns.tolist()[jnho][2]
				print(f" {coeff: .3f}*NHO_{jnho + (nbasis if spin==2 else 0)}({occ:.3f}, F_{frag})", end='')
			print()
		I = np.zeros([nbasis, nbasis])
		I[:P.shape[0], :nhos.shape[1]] = nhos.loc[idx["coeff", :, :]]
		H = np.zeros([nbasis, nbasis])
		H[:nhos.shape[1], :nbos.shape[1]] = nbos.loc["coeff"]
		J = nhos.loc[idx["occ", -1, -1]].to_numpy()
		N = nbos.loc[idx["occ", -1]].to_numpy()
		nho_mwfn.setEnergy([0 for i in range(nbasis)], spin)
		nho_mwfn.setOccupation(J, spin)
		nho_mwfn.setCoefficientMatrix(C @ I, spin)
		nbo_mwfn.setEnergy([0 for i in range(nbasis)], spin)
		nbo_mwfn.setOccupation(N, spin)
		nbo_mwfn.setCoefficientMatrix(C @ I @ H, spin)
	nbo_mwfn.Orthogonalize("GramSchmidt")
	print("Warning: The indeces of orbitals and fragments above start from 0!")
	return nho_mwfn, nbo_mwfn
