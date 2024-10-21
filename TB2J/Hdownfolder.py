import numpy as np

from .io_exchange import ExchangeIO
from .mathutils import get_rotation_arrays
from .io_exchange.structure import get_attribute_array

from lawaf.utils.kpoints import monkhorst_pack, build_Rgrid
from lawaf.mathutils.kR_convert import k_to_R, R_to_k
from lawaf.interfaces.magnon.magnon_downfolder import MagnonDownfolder, MagnonWrapper

def combine_arrays(u, v):

    return np.concatenate([u*v, np.roll(u, -1, axis=-1)*v, np.roll(v, -1, axis=-1)*u], axis=-1)


class ExchangeDownfolder(ExchangeIO):

    def __init__(self, **kwargs):

        reference_axes = kwargs.pop('reference_axes', None) 
        kwargs['kmesh'] = kwargs.pop('kmesh', [7, 7, 7])
        super().__init__(**kwargs)
        self._old_values = None

        if reference_axes is None:
            reference_axes = np.zeros((6, 3))
            reference_axes[:3] = np.eye(3)
            reference_axes[3:] = np.array([
                [0.0, 1.0, 1.0],
                [1.0, 0.0, 1.0],
                [1.0, 1.0, 0.0]
            ]) / np.sqrt(2)
        self.reference_axes = reference_axes

    @property
    def reference_axes(self):
        return self._axes

    @reference_axes.setter
    def reference_axes(self, values):

        axes = get_attribute_array(values, 'reference_axes')
        if axes.ndim != 2 or axes.shape[-1] != 3:
            raise ValueError("The reference axes must be an array of shape (n, 3)")
        self._axes = axes

    @property
    def kpoints(self):
        return self._kpoints

    @ExchangeIO.kmesh.setter
    def kmesh(self, value):

        ExchangeIO.kmesh.fset(self, value)
        self._kpoints = monkhorst_pack(self._kmesh)

    def set_downfolded_magnetic_sites(self, metals):

        vectors = self.vectors
        old_pairs = self.interacting_pairs
        self._old_magnetic_elements = self.magnetic_elements

        self.magnetic_elements = metals
        indices = [old_pairs.index(pair) for pair in self.interacting_pairs]
        self._old_values = self._exchange_values.copy()
        self.set_vectors(values=self.vectors[indices])

    def _generate_u_matrix(self):
        '''
        Constructs the matrix with the coefficients that relate the exchange tensor J_{ij} to the matrix A_{ij}.
        These coefficients only depend on the vectors u_i which depend on the orientations of the magnetic moments.

        '''
        i, j = self.i, self.j
        magmoms = self.magmoms[self._index_spin]
        magmoms /= np.linalg.norm(magmoms, axis=-1)[:, None]

        U, _ = zip(*[get_rotation_arrays(magmoms, u=u[None, :]) for u in self.reference_axes])
        U = np.stack(U).swapaxes(0, 1)

        u1 = combine_arrays(U[i], U[j].conj())
        ur = combine_arrays(U[i], U[j])
        ui = combine_arrays(U[i].conj(), U[j].conj())
        u2 = combine_arrays(U[i].conj(), U[j])
        u = np.concatenate([u1, ur, ui, u2], axis=1)
        
        self.u_matrix = u

    def _compute_AB_coefficients(self, Hk):
        '''
        Computes the coefficients corresponding to A_{ij}(d), B_{ij}(d) from the dynamical matrix h(k). They are
        required to reconstruct the exchange tensor J_{ij}(d).
        These 18 coefficients are stored in an array with the shape (npairs, 18, nvectors)

        '''
        i, j = self.i, self.j
        n = i.max() + 1

        ABk = Hk[:, :, [i, i, i+n, i+n], [j, j+n, j, j+n]]
        ABk = np.moveaxis(ABk, [2, 3], [1, 0]).reshape(len(i), 24, -1)

        exp_summand = np.exp( 2j*np.pi*self.vectors @ self.kpoints.T )
        AB = np.einsum('nik,ndk->nid', ABk[..., ::-1], exp_summand) / len(self.kpoints)

        self.AB_coefficients = AB

    def compute_exchange_tensor(self, Hk):
        '''
        Computes the exchange tensor that best fits the dynamical matrix Hk

        Parameters
        ----------
        Hk : ndarray
            Dynamical matrix corresponding to the points on a k-grid, constructed for different reference axes.
            It has the shape (naxes, nkpoints, 2*natoms, 2*natoms)

        Returns
        -------
        output : ndarray
            An exchange tensor with shape (npairs, ninteractions, 3, 3)
        '''
        n = Hk.shape[-1] // 2
        self.i, self.j = np.triu_indices(n)

        self._generate_u_matrix()
        self._compute_AB_coefficients(Hk)

        ii = np.where(self.i == self.j)
        i0 = np.where( (self.vectors == 0).all(axis=-1) )
        J = np.stack([
            np.linalg.lstsq(mat, coeffs, rcond=None)[0] 
            for mat, coeffs in zip(self.u_matrix, self.AB_coefficients)
        ]).swapaxes(1, 2)
        J = J[:, :, [0, 6, 5, 3, 1, 7, 8, 4, 2]].reshape(len(self.i), -1, 3, 3)
        J *= -1
        J[ii] *= 2
        J[i0] *= 0

        del self.i, self.j, self.u_matrix, self.AB_coefficients

        return J

    def set_exchange_tensor(self, J):

        idig = np.diag_indices(3)
        Jiso = J[:, :, *idig].mean(axis=-1)

        idmi = ([1, 2, 0], [2, 0, 1])
        DMI = (J - J.swapaxes(-1, -2)) / 2

        Jani = (J + J.swapaxes(-1, -2)) / 2
        Jani[:, :, *idig] -= Jiso[:, :, None]

        self._exchange_values[:, :, 3] = Jiso
        self._exchange_values[:, :, 6:9] = DMI[:, :, *idmi]
        self._exchange_values[:, :, 9:] = Jani.reshape(Jani.shape[:2] + (9,))

    def _downfold_matrix(self, Mk, **params):

        MR = k_to_R(kpts=self.kpoints, Rlist=self.Rlist, Mk=Mk, kweights=None)
        model = MagnonWrapper(HR=MR, Rlist=self.Rlist, atoms=self.atoms)

        wann = MagnonDownfolder(model)
        wann.set_parameters(kmesh=self.kmesh, **params)
        ewf = wann.downfold()
        M = np.stack([ewf.get_wann_Hk(k) for k in self.kpoints])

        return M

    def _compute_downfolded_H(self, u, **params):

        Hq = self.Hq(kpoints=self.kpoints, anisotropic=True, u=u)

        n = Hq.shape[-1] // 2
        A = self._downfold_matrix(Hq[:, :n, :n], **params)
        B = self._downfold_matrix(Hq[:, :n, n:], **params)

        ik = self.kpoints.shape[0] // 2
        B[:ik] = B[-1:ik:-1].swapaxes(-1, -2)

        H = np.block([
            [A, B],
            [B.conj().swapaxes(-1, -2), A[::-1].conj()]
        ])

        return H

    def downfold(self, metals, **params):

        try:
            metals = metals.split()
        except AttributeError:
            try:
                metals = list(metals)
            except (ValueError, TypeError):
                raise TypeError("argument must be a list of element symbols.")

        if any(metal not in self.magnetic_elements for metal in metals):
            wrong_symbols = [metal for metal in metals if metal not in self.magnetic_elements]
            raise ValueError(f"The metal symbols '{wrong_symbols}' are not magnetic elements.")
        else:
            magnetic_sites = [symbol for symbol in self.elements if symbol in self.magnetic_elements]
            nsites = len(magnetic_sites)
            metal_indices = [i for i, element in enumerate(magnetic_sites) if element in metals]
            selected_basis = metal_indices #+ [x+nsites for x in metal_indices]
            params |= {'nwann': len(selected_basis), 'selected_basis': selected_basis}

        self.atoms = self.to_ase()
        self.metal_indices = metal_indices
        self.Rlist = build_Rgrid(self.kmesh)

        Hk = np.stack([
            self._compute_downfolded_H(u[None, :], **params) for u in self.reference_axes
        ])

        self.set_downfolded_magnetic_sites(metals)
        J = self.compute_exchange_tensor(Hk)
        self.set_exchange_tensor(J.real)

        del self.atoms, self.metal_indices, self.Rlist

        return J

    def reset(self):

        if self.old_values is not None:
            self.magnetic_elements = self._old_magnetic_elements
            self._exchange_values = self._old_values
            self.old_values = None
            self.old_magnetic_elements = None