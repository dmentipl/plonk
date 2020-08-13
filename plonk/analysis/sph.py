"""SPH sums over neighbours."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable, Dict, Tuple

import numba
import numpy as np
from numba.typed import List

from .._logging import logger
from .._units import Quantity
from ..utils.kernels import (
    kernel_cubic,
    kernel_gradient_cubic,
    kernel_gradient_quintic,
    kernel_gradient_wendland_c4,
    kernel_quintic,
    kernel_wendland_c4,
)

if TYPE_CHECKING:
    from ..snap import Snap


def derivative(
    snap: Snap,
    derivative: str,
    quantity: str,
    kernel: str = 'cubic',
    chunk_size: int = None,
    verbose: bool = False,
) -> Quantity:
    """Calculate derivatives of quantities.

    WARNING: This function is experimental.

    Calculate derivatives such as grad, div, curl on any particle
    quantity using the SPH kernel gradient.

    Parameters
    ----------
    snap
        The Snap object.
    derivative
        Options are 'grad', 'div', 'curl'.
    quantity
        The quantity to take the derivative of. Must be a string to pass
        to Snap object which returns a scalar.
    kernel
        Kernel to compute density. E.g. 'cubic', 'quintic', or
        'Wendland C4'.
    chunk_size : optional
        The size of chunks, in terms of particle number, for neighbour
        finding. If the chunk size is too large then the neighbour
        finding algorithm (scipy.spatial.cKDTree.query_ball_point) runs
        out of memory. Default is None.
    verbose : optional
        If True, print progress. Default is False.

    Returns
    -------
    Quantity
        The derivative of the quantity.
    """
    logger.warning('sph.derivative is experimental')
    logger.info(f'Calculating {derivative}... may take some time...')

    if derivative not in ('grad', 'div', 'curl'):
        raise ValueError('derivative must be in ("grad", "div", "curl")')

    density: Quantity = snap['density']
    quantity_array: Quantity = snap[quantity]

    compute_function = _compute_derivative
    compute_function_kwargs = {
        'density': density.magnitude,
        'quantity_array': quantity_array.magnitude,
    }

    if derivative == 'grad':
        if quantity_array.ndim > 1:
            raise ValueError('Quantity must be scalar for grad')
        result_shape = (len(snap), 3)
        compute_function_kwargs[
            'compute_derivative_over_neighbours'
        ] = _compute_grad_over_neighbours

    elif derivative == 'div':
        if quantity_array.ndim != 2:
            raise ValueError('Quantity must be vector for div')
        result_shape = (len(snap), 1)
        compute_function_kwargs[
            'compute_derivative_over_neighbours'
        ] = _compute_div_over_neighbours

    elif derivative == 'curl':
        if quantity_array.ndim != 2:
            raise ValueError('Quantity must be vector for curl')
        result_shape = (len(snap), 3)
        compute_function_kwargs[
            'compute_derivative_over_neighbours'
        ] = _compute_curl_over_neighbours

    compute_function_kwargs['result_axis_1_size'] = result_shape[1]

    result = summation(
        snap=snap,
        result_shape=result_shape,
        compute_function=compute_function,
        compute_function_kwargs=compute_function_kwargs,
        kernel=kernel,
        chunk_size=chunk_size,
        verbose=verbose,
    )

    logger.info(f'Calculating {derivative}... Done!')

    if derivative == 'div':
        return np.squeeze(result)
    return result


def summation(
    snap: Snap,
    result_shape: Tuple[int, ...],
    compute_function: Callable,
    compute_function_kwargs: Dict[str, Any],
    kernel: str = 'cubic',
    chunk_size: int = None,
    verbose: bool = False,
) -> Quantity:
    """Calculate SPH sums.

    WARNING: This function is experimental.

    Parameters
    ----------
    snap
        The Snap object.
    result_shape
        The shape, as a tuple, of the returned sum.
    compute_function
        The function that computes the sums.
    compute_function_kwargs
        The keyword arguments to the function that computes the sums.
    kernel
        Kernel to compute density. E.g. 'cubic', 'quintic', or
        'Wendland C4'.
    chunk_size : optional
        The size of chunks, in terms of particle number, for neighbour
        finding. If the chunk size is too large then the neighbour
        finding algorithm (scipy.spatial.cKDTree.query_ball_point) runs
        out of memory. Default is None.
    verbose : optional
        If True, print progress. Default is False.

    Returns
    -------
    Quantity
        The result of the SPH summation.
    """
    logger.warning('sph.summation is experimental')
    if kernel not in ('cubic', 'quintic', 'Wendland C4'):
        raise ValueError('Kernel must be in ("cubic", "quintic", "Wendland C4")')
    if kernel == 'cubic':
        kernel_function = kernel_cubic
        kernel_gradient_function = kernel_gradient_cubic
    elif kernel == 'quintic':
        kernel_function = kernel_quintic
        kernel_gradient_function = kernel_gradient_quintic
    elif kernel == 'Wendland C4':
        kernel_function = kernel_wendland_c4
        kernel_gradient_function = kernel_gradient_wendland_c4

    snap.set_kernel(kernel)

    position: Quantity = snap['position']
    smoothing_length: Quantity = snap['smoothing_length']
    mass: Quantity = snap['mass']

    _position = position.magnitude
    _smoothing_length = smoothing_length.magnitude
    _mass = mass.magnitude

    result = np.zeros(result_shape)

    for particle_type in snap.particle_type:
        if verbose:
            logger.info(f'Summing over {particle_type} particles...')

        ind = snap.particle_indices(particle_type)
        if isinstance(ind, list):
            _indices = ind
        else:
            _indices = [ind]

        for type_indices in _indices:
            if len(type_indices) == 0:
                continue
            if chunk_size is not None:
                n_chunks = int(len(type_indices) / chunk_size)
                array_chunks = np.array_split(type_indices, n_chunks)
            else:
                n_chunks = 1
                array_chunks = [type_indices]

            if verbose and n_chunks > 1:
                logger.info(f'Number of chunks: {n_chunks}')
                logger.info(f'Chunk size: {chunk_size}')

            for idx, indices in enumerate(array_chunks):
                if verbose:
                    if n_chunks > 1:
                        logger.info(f'Finding neighbours for chunk: {idx}...')
                    else:
                        logger.info('Finding neighbours...')
                _neighbours = snap[indices].neighbours(np.arange(len(indices)))
                neighbours = List()
                for neigh in _neighbours:
                    neighbours.append(np.array(neigh))

                if verbose:
                    if n_chunks > 1:
                        logger.info(f'Summing over neighbours for chunk: {idx}...')
                    else:
                        logger.info('Summing over neighbours...')

                result[indices] = compute_function(
                    indices=indices,
                    neighbours=neighbours,
                    type_indices=type_indices,
                    position=_position,
                    smoothing_length=_smoothing_length,
                    mass=_mass,
                    kernel_function=kernel_function,
                    kernel_gradient_function=kernel_gradient_function,
                    verbose=verbose,
                    **compute_function_kwargs,
                )

    return result


@numba.njit
def _compute_derivative(
    indices,
    neighbours,
    type_indices,
    position,
    smoothing_length,
    mass,
    kernel_function,
    kernel_gradient_function,
    density,
    quantity_array,
    compute_derivative_over_neighbours,
    result_axis_1_size,
    verbose,
):
    num_particles = len(neighbours)
    result = np.zeros((num_particles, result_axis_1_size))

    for idxi in range(num_particles):
        index = indices[idxi]

        percent = idxi / num_particles * 100
        if np.mod(percent, 10) == 0 and verbose:
            print(percent, '%')

        hi = smoothing_length[index]
        if hi < 0:
            continue

        __neigh = neighbours[idxi]
        _neigh = type_indices[__neigh]
        neigh = _neigh[(_neigh != index) & (smoothing_length[_neigh] > 0)]

        posi = position[index]
        quani = quantity_array[index]

        result[idxi] = compute_derivative_over_neighbours(
            posi=posi,
            hi=hi,
            quani=quani,
            quanj=quantity_array[neigh],
            posj=position[neigh],
            mj=mass[neigh],
            rhoj=density[neigh],
            kernel_gradient_function=kernel_gradient_function,
        )

    return result


@numba.njit
def _compute_grad_over_neighbours(
    posi, hi, quani, quanj, posj, mj, rhoj, kernel_gradient_function,
):
    result = np.zeros(3)

    for idxj in range(len(quanj)):
        dr = posi - posj[idxj]
        rij = np.sqrt(dr[0] ** 2 + dr[1] ** 2 + dr[2] ** 2)
        dr = dr / rij
        qi = rij / hi
        grad_kern = kernel_gradient_function(qi)
        result += (
            mj[idxj] / rhoj[idxj] * (quanj[idxj] - quani) * dr * grad_kern / hi ** 4
        )

    return result


@numba.njit
def _compute_div_over_neighbours(
    posi, hi, quani, quanj, posj, mj, rhoj, kernel_gradient_function,
):
    result = 0.0

    for idxj in range(len(quanj)):
        dq = quanj[idxj] - quani
        dr = posi - posj[idxj]
        rij = np.sqrt(dr[0] ** 2 + dr[1] ** 2 + dr[2] ** 2)
        dr = dr / rij
        qi = rij / hi
        grad_kern = kernel_gradient_function(qi)
        dot = dq[0] * dr[0] + dq[1] * dr[1] + dq[2] * dr[2]
        result += mj[idxj] / rhoj[idxj] * dot * grad_kern / hi ** 4

    return result


@numba.njit
def _compute_curl_over_neighbours(
    posi, hi, quani, quanj, posj, mj, rhoj, kernel_gradient_function,
):
    result = np.zeros(3)
    cross = np.zeros(3)

    for idxj in range(len(quanj)):
        dq = quanj[idxj] - quani
        dr = posi - posj[idxj]
        rij = np.sqrt(dr[0] ** 2 + dr[1] ** 2 + dr[2] ** 2)
        dr = dr / rij
        qi = rij / hi
        grad_kern = kernel_gradient_function(qi)
        cross[0] = dq[1] * dr[2] - dq[2] * dr[1]
        cross[1] = dq[2] * dr[0] - dq[0] * dr[2]
        cross[2] = dq[0] * dr[1] - dq[1] * dr[0]
        result += mj[idxj] / rhoj[idxj] * cross * grad_kern / hi ** 4

    return result
