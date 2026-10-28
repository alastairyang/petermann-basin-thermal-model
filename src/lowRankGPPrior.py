import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as spla

class lowRankGPPrior:
    """
    A class representing a low-rank Gaussian Process (GP) model
    built on graph laplacian. The basis are formed in the null space
    of the simulation PCA basis.
    """

    def __init__(
        self,
        domain_mask, # Boolean mask indicating the active domain region
        x,
        y,
        n_mode, # Number of low-rank modes to retain
        sigma_amp, # Constant, discrepancy amplitude (in normalized space)
        nu, # Matern smoothness
        rho, # Approximate practical range (in grid cells)
    ):
        self.domain_mask = domain_mask
        self.x  = x
        self.y  = y
        self.n_mode = n_mode
        self.sigma_amp = sigma_amp 
        self.nu = nu 
        self.rho = rho 
        self.eigenvec = None    

        self.nx = domain_mask.shape[1]
        self.ny = domain_mask.shape[0]
        self.hx = x[1] - x[0]
        self.hy = y[1] - y[0]
        self.extent = [self.x[0]/1e3, 
                       self.x[-1]/1e3, 
                       self.y[0]/1e3, 
                       self.y[-1]/1e3] # in km
    def grid_laplacian(self, hx=1.0, hy=1.0):
        mask = np.asarray(self.domain_mask, dtype=bool)
        ny, nx = mask.shape

        ids = np.arange(ny * nx).reshape(ny, nx)

        # Horizontal edges
        horizontal = mask[:, :-1] & mask[:, 1:]
        left = ids[:, :-1][horizontal]
        right = ids[:, 1:][horizontal]

        # Vertical edges
        vertical = mask[:-1, :] & mask[1:, :]
        lower = ids[:-1, :][vertical]
        upper = ids[1:, :][vertical]

        rows = np.concatenate([
            left, right,
            lower, upper
        ])

        cols = np.concatenate([
            right, left,
            upper, lower
        ])

        weights = np.concatenate([
            np.full(left.size, 1.0 / hx**2),
            np.full(right.size, 1.0 / hx**2),
            np.full(lower.size, 1.0 / hy**2),
            np.full(upper.size, 1.0 / hy**2),
        ])

        A_full = sparse.coo_matrix(
            (weights, (rows, cols)),
            shape=(nx * ny, nx * ny)
        ).tocsr()

        active = mask.ravel()
        full_ids = np.arange(nx * ny)

        A = A_full[active][:, active]

        degree = np.asarray(A.sum(axis=1)).ravel()
        L = sparse.diags(degree) - A

        return L, np.flatnonzero(active), full_ids

    def eigen_laplacian(self, k=50, hx=1.0, hy=1.0):
        """
        Compute the smallest k eigenvalues and 
        corresponding eigenvectors of the grid Laplacian.
        """

        L, active_ids, full_ids = self.grid_laplacian(hx=hx, hy=hy)

        # Smallest eigenvalues/eigenvectors of a symmetric Laplacian
        gl_eigenvalues, gl_eigenvectors = spla.eigsh(
            L,
            k=k,
            which="SM",
            tol=1e-6,
            maxiter=10000
        )

        return gl_eigenvalues, gl_eigenvectors, active_ids, full_ids

    def sample_graph_matern(
        self,
        gl_eigenvalues,
        gl_eigenvectors,
        dimension=2,
        n_samples=1,
        remove_constant=False,
        rng=None,
    ):
        """
        Generate truncated graph-Matern samples using Laplacian eigenpairs.

        Parameters
        ----------
        eigenvalues : (k,) ndarray
            Laplacian eigenvalues.
        eigenvectors : (n_active, k) ndarray
            Corresponding orthonormal eigenvectors.
        sigma : float
            Target average marginal standard deviation.
        nu : float
            Matern smoothness.
        rho : float
            Approximate practical range, in units consistent with the Laplacian.
            With an unscaled grid Laplacian, this is approximately in grid cells.
        dimension : int
            Spatial dimension. Use 2 for a two-dimensional grid.
        n_samples : int
            Number of independent fields.
        remove_constant : bool
            If True, remove zero/near-zero Laplacian modes.
        rng : np.random.Generator or None
            Random-number generator.

        Returns
        -------
        samples : (n_active, n_samples) ndarray
            Random-field samples on active nodes.
        spectral_variance : (k_used,) ndarray
            Variance assigned to each retained Laplacian mode.
        retained : (k,) boolean ndarray
            Indicator for retained modes.
        raw_variance : (k_used,) ndarray
            Unnormalized Matern spectral variances for each retained mode.
        """
        eigenvalues = np.asarray(gl_eigenvalues, dtype=float)
        eigenvectors = np.asarray(gl_eigenvectors, dtype=float)

        # eigsh does not always return eigenpairs in sorted order
        order = np.argsort(eigenvalues)
        eigenvalues = eigenvalues[order]
        eigenvectors = eigenvectors[:, order]

        # Remove tiny negative values caused by numerical roundoff
        eigenvalues = np.maximum(eigenvalues, 0.0)

        retained = np.ones(eigenvalues.size, dtype=bool)

        if remove_constant:
            tolerance = 1e-10 * max(1.0, eigenvalues.max())
            retained = eigenvalues > tolerance

        lam = eigenvalues[retained]
        U = eigenvectors[:, retained]

        if lam.size == 0:
            raise ValueError("No Laplacian modes remain after filtering.")

        alpha = self.nu + dimension / 2.0
        kappa = np.sqrt(8.0 * self.nu) / self.rho

        # Unnormalized Matern spectral variances
        raw_variance = (kappa**2 + lam) ** (-alpha)

        # For orthonormal U:
        # trace(U diag(w) U.T) = sum(w).
        # Normalize so average marginal variance is sigma**2.
        n_active = U.shape[0]
        normalization = n_active / raw_variance.sum()

        spectral_variance = (
            self.sigma_amp**2 * normalization * raw_variance
        )

        if rng is None:
            rng = np.random.default_rng()

        z = rng.standard_normal((lam.size, n_samples))

        samples = U @ (
            np.sqrt(spectral_variance)[:, None] * z
        )

        return samples, spectral_variance, retained, raw_variance
    def basis(self, plotting=False):
        """
        Create basis function built on Matern kernel Gaussian Process 
        Parameters
        ----------
        hx : float
            Grid spacing in the x direction.
        hy : float
            Grid spacing in the y direction.
        n_mode : int
            Number of discrepancy basis functions to generate.

        Returns
        -------
        F : (n, n_mode) ndarray
            Basis functions generated from the Matern kernel Gaussian Process.

        """
        print("Step 1: Eigen-decomposition on graph Laplacian...")
        gl_eigenval, gl_eigenvec, active_ids, full_ids = self.eigen_laplacian(
            k = self.n_mode,
            hx = self.hx, 
            hy = self.hy
        )
        print("...Done.")


        print("Step 2: Get mode covariance under Matern kernel...")
        samples_active, mode_variance, retained, raw_variance = self.sample_graph_matern(
            gl_eigenvalues  = gl_eigenval,
            gl_eigenvectors = gl_eigenvec,
            n_samples=5,
            remove_constant=False,
            rng=np.random.default_rng(42),
        )
        print("...Done.")

        print("Step 3: scale the basis with variance...")
        # construct the covariance-scaled graph basis
        eigenvec_full = np.zeros((self.nx * self.ny, gl_eigenvec.shape[1]))
        eigenvec_full[active_ids, :] = gl_eigenvec
        eigenvec_full = eigenvec_full.reshape(self.ny, self.nx, gl_eigenvec.shape[1])
        eigenvec_flat = eigenvec_full.reshape(self.ny * self.nx, eigenvec_full.shape[2])
        F = eigenvec_flat @ np.diag(np.sqrt(mode_variance))
        self.matern_basis = F
        return
    def plot_random_fields(
            self,
            fields,
            vmin=0,
            vmax=1,
            save_path=None
        ):
        import matplotlib.pyplot as plt

        n_fields = fields.shape[1]

        fig = plt.figure(figsize=(25, 5))
        gs = fig.add_gridspec(
            1, n_fields + 1,
            width_ratios=[1] * n_fields + [0.05],
            wspace=0.25
        )

        for i in range(n_fields):
            ax = fig.add_subplot(gs[0, i])

            im = ax.imshow(
                fields[:, i].reshape(self.ny, self.nx),
                extent=self.extent,
                cmap='RdBu_r',
                vmin=vmin,
                vmax=vmax,
                origin='lower'
            )
            ax.set_title(f"Field {i}")
            # tilt the ticks
            plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

        # Separate axis: no space is taken from any field plot.
        cax = fig.add_subplot(gs[0, -1])
        fig.colorbar(im, cax=cax)

        if save_path is not None:
            fig.savefig(save_path, bbox_inches='tight')

        plt.show()

    def generate_random_field(self, n):
        """
        Generate random fields from eigenvectors and singular values

        Parameters
        ----------
        n : int
            Number of random fields to generate.

        Returns
        -------
        fields : (n, self.matern_basis.shape[1]) ndarray
            Generated random fields.

        """
        return self.sigma_amp * (self.matern_basis @ np.random.randn(self.matern_basis.shape[1], n))
