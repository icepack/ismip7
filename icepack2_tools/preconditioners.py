r"""Repository-specific Firedrake preconditioners.

Keep Firedrake imports out of :mod:`icepack2_tools.solverconfig`; PETSc loads
this module lazily only when the corresponding Python PC is selected.
"""

from firedrake.slate.slate import AssembledVector
from firedrake.slate.static_condensation.la_utils import (
    LAContext,
    SchurComplementBuilder,
)
from firedrake.slate.static_condensation.scpc import SCPC


class ISMIP7SCPC(SCPC):
    r"""SCPC for ``(u, M, tau)`` with the retained field ordered first.

    Firedrake's stock three-field helper correctly slices arbitrary eliminated
    fields, but passes their original indices into the newly sliced two-field
    tensor.  Eliminating fields ``1,2`` therefore asks that tensor for blocks
    ``(1,2)`` instead of its renumbered blocks ``(0,1)``.  The stock/default
    order (eliminate ``0,1`` and retain ``2``) does not expose the bug.

    Reordering this model's mixed state would touch checkpoint layout and every
    use of ``z.subfunctions``.  This narrow override keeps the established
    ``(u, M, tau)`` order and changes only the two local indices passed to the
    installed Slate ``SchurComplementBuilder``.  Reconstruction remains the
    upstream implementation and uses the original field numbers.
    """

    def condensed_system(self, A, rhs, elim_fields, prefix, pc):
        elim_fields = sorted(map(int, elim_fields))
        if elim_fields != [1, 2]:
            raise ValueError(
                "ISMIP7SCPC requires pc_sc_eliminate_fields=1,2; "
                f"got {elim_fields}"
            )

        blocks = A.blocks
        rhs_blocks = AssembledVector(rhs).blocks
        # Original field 0 (velocity) is retained. Original fields 1:3 become
        # local fields 0:2 after slicing into Aee/Afe/Aef.
        Aff = blocks[0:1, 0:1]
        Aef = blocks[1:3, 0:1]
        Afe = blocks[0:1, 1:3]
        Aee = blocks[1:3, 1:3]
        bf = rhs_blocks[0:1]
        be = rhs_blocks[1:3]

        builder = SchurComplementBuilder(
            prefix,
            Aee,
            Afe,
            Aef,
            pc,
            0,
            1,
            non_zero_saddle_mat=Aff,
        )
        reduced_rhs, reduced_operator = builder.build_schur(
            be, non_zero_saddle_rhs=bf
        )
        return (
            LAContext(
                lhs=reduced_operator,
                rhs=reduced_rhs,
                field_idx=[0],
            ),
            builder,
        )
