r"""Adaptive substepping of the forward's thickness/velocity split step.

The forward advances the thickness with the velocity solved at the previous
geometry, then solves the velocity at the new geometry. That lag makes the
coupled system explicit in the velocity's response to thickness, so it has a
stability limit on the step: where the velocity is hypersensitive to the
thickness (near-flotation grounding zones carrying cell-scale friction
contrast) a step above the limit amplifies a cell-scale perturbation each
step, which shows as a period-2 flip-flop of thickness and speed and then a
blow-up. The limit moves as the geometry evolves, so no fixed step is both
safe and affordable for a 286-year run.

The controller keeps the macro step ``dt`` (forcing, output, checkpoints and
budgets are unchanged) and splits it into ``m`` equal substeps, choosing
``m`` from the backward-Euler local error estimate on the DG0 thickness,

    err_j = tau_j / (tau_j + tau_{j-1}) * max | dh_j - (tau_j / tau_{j-1}) dh_{j-1} |,

the step-to-step change of the thickness increment, taken only over cells
whose increment reversed sign AND grew. The lagged step's instability is a
real negative mode stepped past -2/tau (amplification g = 1 - tau*lam < -1),
so it shows as a sign reversal every substep with a growing amplitude, and
there the estimate is about the size of the increment itself. A reversal
that shrinks (-1 < g < 0) is a damped oscillation the step already handles:
on the same control, a fixed dt of 0.00625 ran stably through years in which
a reversal-only test held the controller at 32 substeps of 0.05 (dt 0.0016)
at four times the cost. Cells whose increment keeps its sign are left out
because a jump
in the thinning rate (a cell going afloat, a neighbour calved) makes the
estimate tau times the jump at any step, and refining only spends substeps
on it without making the trajectory any more stable: on the 2 km SEP1
control one such jump of ~2500 m/yr needed 64 substeps to pass a 1 m
tolerance. A substep whose estimate exceeds the tolerance rejects the macro
step, which the caller rewinds and repeats with more substeps; a run of quiet
macro steps halves ``m`` again. Cells thinner than ``hmin`` in any of the
three states are left out too, so a calving or collapse removal cannot drive
``m`` to its ceiling.

Refining has to pay for itself. Some reversals do not shrink with the step
(a grounding-zone cell flickering across flotation): on the 2 km SEP1 CESM2
historical at 2004.5 the estimate went 1.37, 1.13, 1.12 m at 4, 8 and 16
substeps. When a retry's estimate is still above ``ineffective`` times the
estimate that rejected the previous attempt at this macro step, the retry is
accepted at its count instead of doubling again. The refinement bought
nothing, so the next macro step returns to the count of that rejected attempt
and accepts an estimate no larger than the one that rejected it; a larger one
rejects and tests refining again, and a macro step within the tolerance ends
the rule. A non-finite thickness always rejects, so a real blow-up still
climbs to the ceiling and stalls.

The increment is the transport's own: the thickness the implicit upwind
step leaves, before the positivity floor and the front's removals act on it
(the caller passes that state as ``h_after``, as the free-surface
stabilization takes its tendency). The removals are events, not a time
discretization. Under an advected front a cell the front retreats off keeps
the part of it still behind the front, a fraction set by where the front
lies in the cell and not by the step; a band cell hands its excess back to
the ice it came from; icebergs and collapsed cells are emptied. Read into the
increment, such a jump after a gain is a reversal that grows, tens to
hundreds of metres that no substep count reduces. A cell whose side of an
advected front changed in the substep changes the transport itself, which
exports from the cells inside the front only, so the caller names those cells
(``skip``) and they sit out that comparison and the next.

At the ceiling a finite estimate within ``ceiling`` times the tolerance is
accepted rather than stalling the run. On the 2 km screening of 11 October
every member that stalled did so at one Thwaites front cell, with 1 to 4 m
against a 1 m tolerance at 64 substeps, a time-discretization estimate on a
cell whose mass the transport conserves exactly either way; a run stopped
there has no better step to take. ``ceiling`` 0 restores the stall.
"""

import math

import numpy as np


class SubstepController:
    def __init__(self, *, tol, m_init=1, m_max=64, quiet_steps=20, hmin=10.0,
                 ineffective=0.7, ceiling=20.0, comm=None):
        if tol <= 0:
            raise ValueError("substep tolerance must be positive")
        if ceiling < 0:
            raise ValueError("the ceiling multiple must be nonnegative")
        if not 1 <= m_init <= m_max:
            raise ValueError("need 1 <= initial substeps <= maximum substeps")
        if quiet_steps < 1:
            raise ValueError("quiet window must be at least one macro step")
        self.tol = float(tol)
        self.m = int(m_init)
        self.m_max = int(m_max)
        self.quiet_steps = int(quiet_steps)
        self.hmin = float(hmin)
        self.ineffective = float(ineffective)
        self.comm = comm
        # last accepted substep: its increment, the cells it was taken over
        # and its length; None until the first substep of the run
        self._dh = None
        self._keep = None
        self._tau = None
        self._saved = None
        self._quiet = 0
        self.err_step = 0.0          # largest estimate in the current attempt
        self.err_xy = None           # where the latest estimate was largest
        self.last_reject = None      # estimate that rejected the last attempt
        self.persistent = None       # estimate refining failed to reduce, carried
        self.tolerated_against = None  # estimate the latest tolerated one was held to
        self._m_rejected = None      # count of the last rejected attempt
        self._tolerated_now = False
        self.rejections = 0
        self.tolerated = 0           # retries accepted above tol (ineffective)
        self.ceiling = float(ceiling)  # accept this multiple of tol at m_max
        self.ceiling_accepted = 0    # attempts accepted at the ceiling
        self.ceiling_now = False

    def begin_attempt(self):
        r"""Mark the start of an attempt at a macro step, so a rejection can
        rewind the substep history along with the model state."""
        self._saved = (self._dh, self._keep, self._tau)
        self.err_step = 0.0
        self._tolerated_now = False

    def observe(self, h_before, h_after, tau, xy=None, skip=None):
        r"""Record one substep's thickness change and return its error
        estimate (metres, the same on every rank). The first substep of a
        run has no predecessor and returns 0. ``xy`` (one row of coordinates
        per entry of ``h``) locates the largest estimate in ``err_xy``.
        ``skip`` (one boolean per entry of ``h``) leaves entries out of this
        comparison and the next, as a thin cell is left out: the cells an
        advected front crossed in the substep (see the module docstring)."""
        h_before = np.asarray(h_before)
        h_after = np.asarray(h_after)
        dh = h_after - h_before
        keep = (h_before >= self.hmin) & (h_after >= self.hmin)
        if skip is not None:
            keep &= ~np.asarray(skip, dtype=bool)
        where = None
        if not np.all(np.isfinite(h_after[h_before >= self.hmin])):
            err = math.inf                 # a blown-up state, whatever its history
        elif self._dh is None:
            err = 0.0
        else:
            # cells whose increment reversed sign and grew (per unit time)
            prev = (tau / self._tau) * self._dh
            both = (keep & self._keep & ~(dh * self._dh >= 0)
                    & ~(np.abs(dh) <= np.abs(prev)))
            diff = np.abs(dh[both] - prev[both])
            local = float(np.max(diff)) if diff.size else 0.0
            if not math.isfinite(local):
                local = math.inf
            elif diff.size and xy is not None:
                where = tuple(float(c) for c in np.asarray(xy)[both][int(np.argmax(diff))])
            err = tau / (tau + self._tau) * local
        if self.comm is not None:
            pairs = self.comm.allgather((err, where))
            err, where = max(pairs, key=lambda p: p[0])
        self.err_xy = where
        self._dh, self._keep, self._tau = dh.copy(), keep, float(tau)
        self.err_step = max(self.err_step, err)
        return err

    def exceeds(self, err):
        r"""Whether ``err`` rejects the attempt: above the tolerance, unless a
        refinement already failed to reduce the estimate (see the module
        docstring). NaN and inf always reject."""
        if not math.isfinite(err):
            return True
        if err <= self.tol:
            return False
        if self.last_reject is not None:
            tolerate = err > self.ineffective * self.last_reject
            ref = self.last_reject
        else:
            tolerate = self.persistent is not None and err <= self.persistent
            ref = self.persistent
        if tolerate:
            self.tolerated += 1
            self.tolerated_against = ref
            self._tolerated_now = True
            return False
        if self.m >= self.m_max and self.ceiling > 0 and err <= self.ceiling * self.tol:
            self.ceiling_accepted += 1
            self.ceiling_now = True
            return False
        return True

    def refine(self, err=None):
        r"""Rewind the substep history and return the substep count for the
        retry, or None when the ceiling is already reached. ``err`` is the
        estimate that rejected the attempt; None means the velocity solve
        failed, which doubles ``m``. An estimate jumps straight to the count
        that would meet the tolerance if the error scaled as tau^2, at most
        eightfold at once."""
        self._dh, self._keep, self._tau = self._saved
        self._quiet = 0
        self.rejections += 1
        self.last_reject = err
        self._m_rejected = self.m
        if self.m >= self.m_max:
            return None
        factor = 2
        if err is not None and math.isfinite(err) and err > 0:
            need = math.sqrt(err / self.tol) / 0.9
            factor = min(8, max(2, 2 ** math.ceil(math.log2(need))))
        self.m = min(self.m_max, self.m * factor)
        return self.m

    def accept(self):
        r"""Close an accepted macro step. After ``quiet_steps`` accepted macro
        steps in a row, none rejected and none with an estimate above a
        quarter of the tolerance (what halving m roughly multiplies it by),
        the substep count halves. An attempt accepted because refining did
        not reduce its estimate returns m to the count of the attempt it
        rejected and carries that estimate to the next macro steps. Returns
        True when the quiet steps halved m."""
        if self._tolerated_now and self.last_reject is not None:
            self.persistent = self.last_reject
            self.m = self._m_rejected
        elif not self._tolerated_now and self.err_step <= self.tol:
            self.persistent = None
        self.last_reject = None
        self._m_rejected = None
        self._tolerated_now = False
        self.ceiling_now = False
        if self.err_step < 0.25 * self.tol:
            self._quiet += 1
        else:
            self._quiet = 0
        if self._quiet >= self.quiet_steps and self.m > 1:
            self.m //= 2
            self._quiet = 0
            return True
        return False
